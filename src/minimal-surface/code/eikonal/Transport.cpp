#include "stdafx.h"
#include "Transport.h"
#include "SimpleITK.h"
#include "Utils.h"
namespace sitk = itk::simple;
using namespace std;


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
	int src_idx = mIterationCount++ % 2;
	int dst_idx = (src_idx + 1) % 2;
	const sitk::Image& trf = mTransportFunctionMap[src_idx];// mTransportFunctionMap[bev];
	sitk::Image& trft = mTransportFunctionMap[dst_idx];// mTransportFunctionMap[targi];
	const double* trf_buffer = trf.GetBufferAsDouble();
	double* trft_buffer = trft.GetBufferAsDouble();
	sitk::Image& bound = mReadOnlyMap;
	int* bound_buffer = bound.GetBufferAsInt32();
	vector<uint32_t> size = trf.GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	static const realnum alpha = 0.003;// small value

	double* gx_buffer = mGradX.GetBufferAsDouble();
	double* gy_buffer = mGradY.GetBufferAsDouble();
	double* gz_buffer = mGradZ.GetBufferAsDouble();

OMP_PARALLEL_FOR
	for (int zyx = 0; zyx < zs * ys * xs; zyx++) {
		int zz = zyx / (ys * xs);
		if (zz > 0 && zz < zs - 1)
		{
			int yy = (zyx / xs) % ys;
			if (yy > 0 && yy < ys - 1)
			{
				int xx = zyx % xs;
				if (xx > 0 && xx < xs - 1)
				{
					if (BUF_IDX3D(bound_buffer, xs, ys, zs, xx, yy, zz) != 0) {
						continue; // boundary to completely omit (no transport from uninitialized distance map)
					}

					realnum val;
					int x, y, z;
					realnum gx(BUF_IDX3D(gx_buffer, xs, ys, zs, xx, yy, zz)), gy(BUF_IDX3D(gy_buffer, xs, ys, zs, xx, yy, zz)), gz(BUF_IDX3D(gz_buffer, xs, ys, zs, xx, yy, zz));
					x = gx > 0 ? xx + 1 : xx - 1;
					y = gy > 0 ? yy + 1 : yy - 1;
					z = gz > 0 ? zz + 1 : zz - 1;
					if (gx < 0) gx *= -1;
					if (gy < 0) gy *= -1;
					if (gz < 0) gz *= -1;
					val = gx * BUF_IDX3D(trf_buffer, xs, ys, zs, x, yy, zz) + gy * BUF_IDX3D(trf_buffer, xs, ys, zs, xx, y, zz) + gz * BUF_IDX3D(trf_buffer, xs, ys, zs, xx, yy, z);
					val /= gx + gy + gz + alpha;// zero: at points where calculation is needed, 1 read only boundary points
					if (abs(val - BUF_IDX3D(trft_buffer, xs, ys, zs, xx, yy, zz)) > abs(BUF_IDX3D(trft_buffer, xs, ys, zs, xx, yy, zz)*0.1))
						updated = true;
					BUF_IDX3D(trft_buffer, xs, ys, zs, xx, yy, zz) = val;  // update wherever calculation is needed
				}
			}
		}
	}
	return updated;
}

void TransportFunction::Calculate(sitk::Image distanceMap, const sitk::Image& initialSlice, int initialSliceIdx, realnum distanceMax)
{
	_PROFILING;
	Initialize(distanceMap, initialSlice, initialSliceIdx, distanceMax);
	if (!m_active) return;

	for (int ii = 0; ii < mIterationsMax; ++ii) {
		bool updated = Iterate();
		if (!updated)
			break;
	}
}

void TransportFunction::RotateTransportFunction(const std::vector<double>& rotation_matrix, vector<unsigned> sample_size) {
	resample_img(mTransportFunctionMap[0], mTransportFunctionMap[0], rotation_matrix, sample_size, true);
	resample_img(mTransportFunctionMap[1], mTransportFunctionMap[1], rotation_matrix, sample_size, true);
}

const sitk::Image& TransportFunction::GetTransportFunction() const {
	return mTransportFunctionMap[mIterationCount % 2];
}

sitk::Image& TransportFunction::GetTransportFunction() {
	return mTransportFunctionMap[mIterationCount % 2];
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
