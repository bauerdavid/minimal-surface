#include "Transport.h"
#include "SimpleITK.h"
#include "Utils.h"
#include <sitkImageOperators.h>
namespace sitk = itk::simple;
using namespace std;

#ifdef ORDERED_TRAVERSAL

bool TransportFunction::Finished() {
	return mNarrowBand.empty();
}

double TransportFunction::ComputeValue(std::vector<int> p) {
	int xx(p[0]), yy(p[1]), zz(p[2]);
	double* grad_buffers[] = { mGradX.GetBufferAsDouble(), mGradY.GetBufferAsDouble() , mGradZ.GetBufferAsDouble() };
	double* transport_buffer = mTransportFunctionMap[0].GetBufferAsDouble();
	int32_t* frozen_buffer = mReadOnlyMap.GetBufferAsInt32();
	vector<unsigned> size = mGradX.GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	realnum val = 0;
	realnum grad_sum = -0.003;
	for (int i = 0; i < 3; i++) {
		int coord = p[i];
		std::vector<int> neighbor_pos = p;
		realnum grad = BUF_IDX3D(grad_buffers[i], xs, ys, zs, xx, yy, zz);
		int neigh_coord = grad > 0 ? coord + 1 : coord - 1;
		neighbor_pos[i] = neigh_coord;
		if (!BUF_IDX3D(frozen_buffer, xs, ys, zs, neighbor_pos[0], neighbor_pos[1], neighbor_pos[2])) continue;
		if (grad < 0) grad *= -1;
		realnum neighbor_value = BUF_IDX3D(transport_buffer, xs, ys, zs, neighbor_pos[0], neighbor_pos[1], neighbor_pos[2]);
		val += grad * neighbor_value;
		grad_sum += grad;
	}
	val /= grad_sum;
	return val;
}


void TransportFunction::Iterate() {
	_PROFILING;
	POINT3D p;
	double dist;
	tie(p, dist) = mNarrowBand.top();
	auto [x, y, z] = representation_to_point<int>(p);
	mNarrowBand.pop();
	int32_t* frozen_buffer = mReadOnlyMap.GetBufferAsInt32();
	double* transport_buffer = mTransportFunctionMap[0].GetBufferAsDouble();
	vector<unsigned> size = mTransportFunctionMap[0].GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	BUF_IDX3D(frozen_buffer, xs, ys, zs, x, y, z) = 1;
	BUF_IDX3D(transport_buffer, xs, ys, zs, x, y, z) = dist;
	//OMP_PARALLEL_FOR
	for (int i = 0; i < 6; i++) {
		const vector<int>& offset = NEIGH6_OFFSET[i];
		int xn = x + offset[0], yn = y + offset[1], zn = z + offset[2];
		if (xn < 0 || xn >= xs || yn < 0 || yn >= ys || zn < 0 || zn >= zs
			|| BUF_IDX3D(frozen_buffer, xs, ys, zs, xn, yn, zn)) continue;
		double dist_n = ComputeValue({ xn, yn, zn });
		mNarrowBand.push_or_promote(point_to_representation(xn, yn, zn), dist_n);
	}
}

void TransportFunction::Initialize(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax)
{
	_PROFILING;
	vector<unsigned> size = distanceMap.GetSize();
	vector<unsigned> inimap_size = initialSlice.GetSize();
	if (size[1] != inimap_size[0] && size[2] != inimap_size[1]) return;
	double* distmap_buffer = distanceMap.GetBufferAsDouble();
	int xs = size[0], ys = size[1], zs = size[2];

	mReadOnlyMap = sitk::Image(size, sitk::sitkInt32);
	mTransportFunctionMap[0] = sitk::Image(size, sitk::sitkFloat64);
	double min_val = GetImageMinimum(initialSlice);
	for (int zz = 0; zz < zs; ++zz)
		for (int yy = 0; yy < ys; ++yy)
			for (int xx = 0; xx < xs; ++xx) {
				if (distmap_buffer[xx + xs * (yy + ys * zz)] < -0.5) // not calculated
					distmap_buffer[xx + xs * (yy + ys * zz)] = 1.2 * distanceMax; // new
				if (!xx || xx == xs - 1 || !yy || yy == ys - 1 || !zz || zz == zs - 1) // external borders
					distmap_buffer[xx + xs * (yy + ys * zz)] = 1.2 * distanceMax; // new
			}
	mGradX = sitk::Derivative(distanceMap, 0, 1, false);
	mGradY = sitk::Derivative(distanceMap, 1, 1, false);
	mGradZ = sitk::Derivative(distanceMap, 2, 1, false);
	double* gx_buffer = mGradX.GetBufferAsDouble();
	double* gy_buffer = mGradY.GetBufferAsDouble();
	double* gz_buffer = mGradZ.GetBufferAsDouble();
	const double* inimap_buffer = initialSlice.GetBufferAsDouble();
	double* trf0_buffer = mTransportFunctionMap[0].GetBufferAsDouble();
	int* frozen_buffer = mReadOnlyMap.GetBufferAsInt32();
	vector<POINT3D> points;
	vector<double> values;
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				if (xx == initialSliceIdx) { // _NO_BOU_ = not boundary, < _NO_BOU_ = calculated boundary
					BUF_IDX3D(trf0_buffer, xs, ys, zs, xx, yy, zz) = BUF_IDX2D(inimap_buffer, ys, zs, yy, zz); // passed arrival plane
					BUF_IDX3D(frozen_buffer, xs, ys, zs, xx, yy, zz) = 1; // 1: read only boundary
					points.push_back(point_to_representation(xx, yy, zz));
					values.push_back(BUF_IDX2D(inimap_buffer, ys, zs, yy, zz));
				}
				else if (!xx || xx == xs - 1 || !yy || yy == ys - 1 || !zz || zz == zs - 1) { // external borders
					BUF_IDX3D(trf0_buffer, xs, ys, zs, xx, yy, zz) = min_val;
					BUF_IDX3D(frozen_buffer, xs, ys, zs, xx, yy, zz) = 1; // 1: read only boundary
					points.push_back(point_to_representation(xx, yy, zz));
					values.push_back(min_val);
				}

			}
		}
	}
	mNarrowBand = IndexedPriorityQueue<POINT3D_MAP(double), std::greater<double>>(points.begin(), points.end(), values.begin(), true);
	m_active = 1;
}

void TransportFunction::Calculate(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax) {
	Initialize(distanceMap, initialSlice, initialSliceIdx, distanceMax);
	while (!Finished()) {
		Iterate();
	}
}
#else
void TransportFunction::Initialize(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax)
{
	_PROFILING;
	mIterationCount = 0;
	vector<unsigned> size = distanceMap.GetSize();
	vector<unsigned> inimap_size = initialSlice.GetSize();
	if (size[1] != inimap_size[0] && size[2] != inimap_size[1]) return;
	double* distmap_buffer = distanceMap.GetBufferAsDouble();
	int xs = size[0], ys = size[1], zs = size[2];

	mReadOnlyMap = sitk::Image(size, sitk::sitkInt32);
	mTransportFunctionMap[0] = sitk::Image(size, sitk::sitkFloat64);
	mTransportFunctionMap[1] = sitk::Image(size, sitk::sitkFloat64);
	double min_val = GetImageMinimum(initialSlice);
	for (int zz = 0; zz < zs; ++zz)
		for (int yy = 0; yy < ys; ++yy)
			for (int xx = 0; xx < xs; ++xx) {
				if (distmap_buffer[xx + xs * (yy + ys * zz)] < -0.5) // not calculated
					distmap_buffer[xx + xs * (yy + ys * zz)] = 1.2 * distanceMax; // new
				if (!xx || xx == xs - 1 || !yy || yy == ys - 1 || !zz || zz == zs - 1) // external borders
					distmap_buffer[xx + xs * (yy + ys * zz)] = 1.2 * distanceMax; // new
			}
	mGradX = sitk::Derivative(distanceMap, 0, 1, false);
	mGradY = sitk::Derivative(distanceMap, 1, 1, false);
	mGradZ = sitk::Derivative(distanceMap, 2, 1, false);
	const double* inimap_buffer = initialSlice.GetBufferAsDouble();
	double* trf0_buffer = mTransportFunctionMap[0].GetBufferAsDouble();
	double* trf1_buffer = mTransportFunctionMap[1].GetBufferAsDouble();
	int* isboundary_buffer = mReadOnlyMap.GetBufferAsInt32();
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				if (xx == initialSliceIdx) { // _NO_BOU_ = not boundary, < _NO_BOU_ = calculated boundary
					BUF_IDX3D(trf1_buffer, xs, ys, zs, xx, yy, zz) = BUF_IDX3D(trf0_buffer, xs, ys, zs, xx, yy, zz) = BUF_IDX2D(inimap_buffer, ys, zs, yy, zz); // passed arrival plane
					BUF_IDX3D(isboundary_buffer, xs, ys, zs, xx, yy, zz) = 1; // 1: read only boundary
				}
				else if (!xx || xx == xs - 1 || !yy || yy == ys - 1 || !zz || zz == zs - 1) { // external borders
					BUF_IDX3D(trf1_buffer, xs, ys, zs, xx, yy, zz) = BUF_IDX3D(trf0_buffer, xs, ys, zs, xx, yy, zz) = min_val;
					BUF_IDX3D(isboundary_buffer, xs, ys, zs, xx, yy, zz) = 1; // 1: read only boundary
				}
			}
		}
	}
	m_active = 1;
}

bool TransportFunction::Iterate()
{
	bool updated = false;
	_PROFILING;
	int src_idx = ++mIterationCount % 2;
	int dst_idx = (src_idx + 1) % 2;
	double* dst_buffer = mTransportFunctionMap[dst_idx].GetBufferAsDouble();
	const double* src_buffer = mTransportFunctionMap[src_idx].GetBufferAsDouble();
	sitk::Image& bound = mReadOnlyMap;
	int* bound_buffer = bound.GetBufferAsInt32();
	vector<uint32_t> size = mTransportFunctionMap[src_idx].GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	static const realnum alpha = 0.003;// small value

	double* gx_buffer = mGradX.GetBufferAsDouble();
	double* gy_buffer = mGradY.GetBufferAsDouble();
	double* gz_buffer = mGradZ.GetBufferAsDouble();

	//bool decreasing_x = mIterationCount % 2;
	//bool decreasing_y = mIterationCount / 2 % 2;
	//bool decreasing_z = mIterationCount / 4 % 2;

OMP_PARALLEL_FOR_NUM_THREADS(4)
	for (int zyx = 0; zyx < zs * ys * xs; zyx++) {
		int zz = zyx / (ys * xs);
		if (zz > 0 && zz < zs - 1)
		{
			//if (decreasing_z) zz = zs - 1 - zz;
			int yy = (zyx / xs) % ys;
			if (yy > 0 && yy < ys - 1)
			{
				//if (decreasing_y) yy = ys - 1 - yy;
				int xx = zyx % xs;
				if (xx > 0 && xx < xs - 1)
				{
					//if (decreasing_x) xx = xs - 1 - xx;
					if (BUF_IDX3D(bound_buffer, xs, ys, zs, xx, yy, zz) != 0) {
						continue; // boundary to completely omit (no transport from uninitialized distance map)
					}

					realnum val;
					realnum dx, dy, dz; // reference values in the transport function along the given axis
					int x, y, z;
					realnum gx(BUF_IDX3D(gx_buffer, xs, ys, zs, xx, yy, zz)), gy(BUF_IDX3D(gy_buffer, xs, ys, zs, xx, yy, zz)), gz(BUF_IDX3D(gz_buffer, xs, ys, zs, xx, yy, zz));
					x = gx > 0 ? xx + 1 : xx - 1;
					y = gy > 0 ? yy + 1 : yy - 1;
					z = gz > 0 ? zz + 1 : zz - 1;
					if (gx < 0) gx *= -1;
					if (gy < 0) gy *= -1;
					if (gz < 0) gz *= -1;
					dx = BUF_IDX3D(src_buffer, xs, ys, zs, x, yy, zz);
					dy = BUF_IDX3D(src_buffer, xs, ys, zs, xx, y, zz);
					dz = BUF_IDX3D(src_buffer, xs, ys, zs, xx, yy, z);
					val = gx * dx + gy * dy + gz * dz;
					val /= gx + gy + gz + alpha;// zero: at points where calculation is needed, 1 read only boundary points
					if (abs(val - BUF_IDX3D(src_buffer, xs, ys, zs, xx, yy, zz)) > abs(BUF_IDX3D(src_buffer, xs, ys, zs, xx, yy, zz) * 0.001))
						updated = true;
					BUF_IDX3D(dst_buffer, xs, ys, zs, xx, yy, zz) = val;
				}
			}
		}
	}
	return updated;
}

void TransportFunction::Calculate(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax, int maxIterations)
{
	_PROFILING;
	Initialize(distanceMap, initialSlice, initialSliceIdx, distanceMax);
	if (!m_active) return;
	for (int ii = 0; ii < maxIterations; ++ii) {
		bool updated = Iterate();
		if (!updated)
			return;
	}
	cout << "transport function did not converge" << endl;
}

#endif
void TransportFunction::RotateTransportFunction(const std::vector<double>& rotation_matrix, vector<unsigned> sample_size) {
	resample_img(mTransportFunctionMap[0], mTransportFunctionMap[0], rotation_matrix, sample_size, true);
	resample_img(mTransportFunctionMap[1], mTransportFunctionMap[1], rotation_matrix, sample_size, true);
}

sitk::Image& TransportFunction::GetTransportFunction() {
	return mTransportFunctionMap[(mIterationCount + 1) % 2];
}

const sitk::Image& TransportFunction::GetTransportFunction() const {
	return mTransportFunctionMap[(mIterationCount + 1) % 2];
}

const sitk::Image& TransportFunction::GetTransportFunction(int i) const {
	return mTransportFunctionMap[i];
}

sitk::Image& TransportFunction::GetTransportFunction(int i) {
	return mTransportFunctionMap[i];
}

const sitk::Image& TransportFunction::GetReadOnlyMap() const {
	return mReadOnlyMap;
}

const sitk::Image& TransportFunction::GetGradX() const {
	return mGradX;
}

const sitk::Image& TransportFunction::GetGradY() const {
	return mGradY;
}

const sitk::Image& TransportFunction::GetGradZ() const {
	return mGradZ;
}
