#include "AreaEikonal.h"
#include "Utils.h"
#include "math.h"
#include <omp.h>
#include "SimpleITK.h"
#include <vector>
#include <sitkImageOperators.h>
#include "FMSigned.h"

#ifdef DEBUG_CURVATURE
#include "Debug.h"
#include <string>
#endif
#define MAX_THREADS 32
//							x  210210210210210210210210210
//							y  222111000222111000222111000
//							z  222222222111111111000000000
#define CONNECTIVITY_MASK_19 0b010111010111111111010111010
#define CONNECTIVITY_MASK_6  0b000010000010101010000010000
#define REFRESH_ITERS 50
#define MAX_DIST_FROM_PLANE 10
#define MAX_FIELD_VAL 0.95
#define M -1


inline double calculate_speed(double curv, double phi, double S = 1) {
	if (curv < -phi / (4 * S))
		//return 1 / phi;
		return M * (curv + phi / (4 * S)) + 2*S / phi;
	if (abs(curv) < 1e-7)
		return S / phi;
	return (1 - sqrt(1 + (4 * curv * S) / phi)) / (-2 * curv);
}
inline bool is_outlier(double curv, double s = 2) {
	return curv < -1 / (2 * s);
}
inline double calculate_S(double curv, double phi, double s=2) {
	if (curv < -1 / (2 * s)) {
		double smc = s - M * curv;
		return phi / 4 * (smc + sqrt(smc * smc - 2 * M));
		//return (M * phi) / (4 * (s - M * curv - 2 / phi));
	}
	return (curv*s*s + s) * phi;
}

template<typename T>
inline T squared_sum(T x, T y, T z) {
	return x * x + y * y + z * z;
}
namespace sitk = itk::simple;
using namespace std;


void AreaEikonal::SmoothDistances() {
	mDistanceMap[0] = sitk::DiscreteGaussian(mDistanceMap[0]);
	mDistanceMap[1] = sitk::DiscreteGaussian(mDistanceMap[1]);
}

void AreaEikonal::CombineDistance() {
	vector<unsigned> size = mMeetingPointsMap.GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	double* distance0_buffer = mDistanceMap[0].GetBufferAsDouble();
	double* distance1_buffer = mDistanceMap[1].GetBufferAsDouble();
	double* combined_buffer = mCombinedDistanceMap.GetBufferAsDouble();
OMP_PARALLEL_FOR
	for (int zyx = 0; zyx < zs * ys * xs; zyx++)
	{
		int xx = zyx % xs;
		double d0 = distance0_buffer[zyx];
		double d1 = distance1_buffer[zyx];
		double newval;
		if (d1 < 0 || d0 < d1 && d0 >= -1e-5)
			newval = d0;
		else
			newval = d1;
		combined_buffer[zyx] = newval;
	}
}

void AreaEikonal::FitMeetingPlane() {
	_PROFILING;
	vector<unsigned> size= mMeetingPointsMap.GetSize();
	double* distance0_buffer = mDistanceMap[0].GetBufferAsDouble();
	double* distance1_buffer = mDistanceMap[1].GetBufferAsDouble();
	int* meeting_plane_buffer = mMeetingPointsMap.GetBufferAsInt32();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	// The points each thread finds. They go into mMeetingPoints in thread order below, not in
	// the order the threads finish: the plane is fitted to the points in the set's iteration
	// order, which depends on the order of insertion, so a varying order changed the result in
	// the last bits from run to run. Each thread's points are a contiguous range of the volume,
	// so this is the order a single thread would insert them in.
	vector<vector<POINT3D>> found(omp_get_max_threads());
OMP_PARALLEL
	{
		vector<POINT3D>& temp = found[omp_get_thread_num()];
OMP_FOR
		for (int zyx = 0; zyx < zs * ys * xs; zyx++)
		{
			int zz = zyx / (ys * xs);
			{
				int yy = (zyx / xs) % ys;
				{
					int xx = zyx % xs;
					if (xx == 0 || xx == xs - 1 || yy == 0 || yy == ys - 1 || zz == 0 || zz == zs - 1) {
						int sample_x(xx), sample_y(yy), sample_z(zz);
						if (xx == 0)
							sample_x = 1;
						else if (xx == xs - 1)
							sample_x = xs - 2;
						if (yy == 0)
							sample_y = 1;
						else if (yy == ys - 1)
							sample_y = ys - 2;
						if (zz == 0)
							sample_z = 1;
						else if (zz == zs - 1)
							sample_z = zs - 2;
						distance0_buffer[zyx] = distance0_buffer[sample_x + xs * (sample_y + ys * sample_z)];
						distance1_buffer[zyx] = distance1_buffer[sample_x + xs * (sample_y + ys * sample_z)];
					}
					if (meeting_plane_buffer[zyx] > 0)
						temp.push_back(point_to_representation(xx, yy, zz));


				}
			}
		}
	}
	for (const vector<POINT3D>& temp : found)
		mMeetingPoints.insert(temp.begin(), temp.end());
	vector<vector<double>> meeting_plane_vec;
	meeting_plane_vec.reserve(mMeetingPoints.size());

	std::transform(mMeetingPoints.begin(), mMeetingPoints.end(), back_inserter(meeting_plane_vec),
#ifdef STORE_POINT_AS_INTEGER
		[](unsigned long long val) {
			
			auto [x, y, z] = representation_to_point<double>(val);
			vector<double> v({x, y, z});
			return v;
		}
#else
		[](const Vec3<int>& val) {

			return vector<double>({ (double)val.x, (double)val.y, (double)val.z });
		}
#endif
	);
	std::pair<Vec3<double>, Vec3<double>> plane_info = best_plane_from_points(meeting_plane_vec);
	mMeetingPlaneCenter = plane_info.first;
	mMeetingPlaneNormal = plane_info.second;
	double max_norm_val = abs(mMeetingPlaneNormal.x());
	int max_norm_sign = sgn(mMeetingPlaneNormal.x());
	if (abs(mMeetingPlaneNormal.y()) > max_norm_val) {
		max_norm_val = abs(mMeetingPlaneNormal.y());
		max_norm_sign = sgn(mMeetingPlaneNormal.y());
	}
	if (abs(mMeetingPlaneNormal.z()) > max_norm_val) {
		max_norm_val = abs(mMeetingPlaneNormal.z());
		max_norm_sign = sgn(mMeetingPlaneNormal.z());
	}
	mMeetingPlaneNormal *= max_norm_sign;
	UpdateMeetingPlaneInfo();

	vector<double> sample_plane_normal = { 1, 0, 0 };
	// calculate rotation between the two normals
	vector<double> meetingPlaneNormal = vector<double>(this->mMeetingPlaneNormal.begin(), this->mMeetingPlaneNormal.end());
}

void AreaEikonal::UpdateMeetPoints() {
	_PROFILING;
	vector<unsigned> size = mPhaseFieldMap[0].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	int* meeting_plane_buffer = mMeetingPointsMap.GetBufferAsInt32();
	for (int ii = 0; ii < 2; ii++) {
		sitk::Image& distance = mDistanceMap[ii];
		sitk::Image& counterd = mDistanceMap[(ii + 1) % 2];
		double* distance_buffer = distance.GetBufferAsDouble();
		double* counterd_buffer = counterd.GetBufferAsDouble();
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
						int zp = zz + 1, zm = zz - 1;
						int yp = yy + 1, ym = yy - 1;
						int xp = xx + 1, xm = xx - 1;


						if (counterd_buffer[xx + xs * (yy + ys * zz)] >= 0 && !BUF_IDX3D(meeting_plane_buffer, xs, ys, zs, xx, yy, zz)) {
							if (distance_buffer[xx + xs * (yy + ys * zz)] - counterd_buffer[xx + xs * (yy + ys * zz)] >= 0 &&
								((distance_buffer[xx + xs * (yp + ys * zz)] > 0 && distance_buffer[xx + xs * (yp + ys * zz)] - counterd_buffer[xx + xs * (yp + ys * zz)] < 0) ||
									(distance_buffer[xx + xs * (ym + ys * zz)] > 0 && distance_buffer[xx + xs * (ym + ys * zz)] - counterd_buffer[xx + xs * (ym + ys * zz)] < 0) ||
									(distance_buffer[xp + xs * (yy + ys * zz)] > 0 && distance_buffer[xp + xs * (yy + ys * zz)] - counterd_buffer[xp + xs * (yy + ys * zz)] < 0) ||
									(distance_buffer[xm + xs * (yy + ys * zz)] > 0 && distance_buffer[xm + xs * (yy + ys * zz)] - counterd_buffer[xm + xs * (yy + ys * zz)] < 0) ||
									(distance_buffer[xx + xs * (yy + ys * zp)] > 0 && distance_buffer[xx + xs * (yy + ys * zp)] - counterd_buffer[xx + xs * (yy + ys * zp)] < 0) ||
									(distance_buffer[xx + xs * (yy + ys * zm)] > 0 && distance_buffer[xx + xs * (yy + ys * zm)] - counterd_buffer[xx + xs * (yy + ys * zm)] < 0))) {
								meeting_plane_buffer[xx + xs * (yy + ys * zz)] = 1;
OMP_ATOMIC
								mMeetingPointsCount++;
							}
						}
					}
				}
			}
		}
	}
}

void AreaEikonal::UpdateMeetingPlane() {
	int prev_point_count = mMeetingPointsCount >= mMeetingPointsMin ? mMeetingPointsCount : mMeetingPointsMin;
	UpdateMeetPoints();
	if (mMeetingPointsCount > prev_point_count)
		FitMeetingPlane();
	else if (mMeetingPointsCount == prev_point_count && mMeetingPointsCount > 0) {
		mPlaneFinalized = true;
	}
}

bool AreaEikonal::IsInitialized(){
    return mInitialized;
}

void AreaEikonal::InitializeContainers(const sitk::Image& image) {
	vector<uint32_t> size = image.GetSize();
	mMeetingPointsMap = sitk::Image(size, sitk::sitkInt32);


	mPhaseFieldMap[0] = sitk::Image(size, sitk::sitkFloat64) - 1;
	mPhaseFieldMap[1] = sitk::Image(size, sitk::sitkFloat64) - 1;
	mDistanceMap[0] = sitk::Image(size, sitk::sitkFloat64) - 1;
	mDistanceMap[1] = sitk::Image(size, sitk::sitkFloat64) - 1;
	mCombinedDistanceMap = sitk::Image(size, sitk::sitkFloat64) - 1;

	mSampleImage = sitk::Image(size, sitk::sitkFloat64);

	mCurvatureMap[0] = sitk::Image(size, sitk::sitkFloat64);
	mCurvatureMap[1] = sitk::Image(size, sitk::sitkFloat64);
	mCurvatureBand[0].clear();
	mCurvatureBand[1].clear();
#ifdef DEBUG_CURVATURE

	m_new_update[0] = sitk::Image(size, sitk::sitkUInt8);
	m_new_update[1] = sitk::Image(size, sitk::sitkUInt8);
#endif
    mInitialized = true;
}

void AreaEikonal::InitializeDistanceMap(int idx, Vec3<double> center, int initialRadius) {
	vector<uint32_t> size = mPhaseFieldMap[0].GetSize();
	int spacex(size[0]), spacey(size[1]), spacez(size[2]);
	int hs = initialRadius + 2;
	mCurrentDistance = 0;
	double* field_buffer = mPhaseFieldMap[idx].GetBufferAsDouble();
	double* distance_buffer = mDistanceMap[idx].GetBufferAsDouble();
	for (int zz = (int)center.z() - hs; zz < (int)center.z() + hs; ++zz) {
		if (zz < 0 || zz >= spacez) continue;
		for (int yy = (int)center.y() - hs; yy < (int)center.y() + hs; ++yy) {
			if (yy < 0 || yy >= spacey) continue;
			for (int xx = (int)center.x() - hs; xx < (int)center.x() + hs; ++xx) {
				if (xx < 0 || xx >= spacex) continue;
				int dx = xx - (int)center.x(), dy = yy - (int)center.y();
				int dz = zz - (int)center.z();
				int dd = dx * dx + dy * dy + dz * dz;
				if (dd <= initialRadius * initialRadius) {
					double dd_sqrt = sqrt(dd);
					if (dd > (initialRadius - 1) * (initialRadius - 1)) {
						field_buffer[xx + spacex * (yy + spacey * zz)] = initialRadius - dd_sqrt;
#ifdef USE_VECTOR_AS_SET
						mActivePoints[idx].push_back(point_to_representation(xx, yy, zz));
#else
						mActivePoints[idx].insert(point_to_representation(xx, yy, zz));
#endif
					}
					else {
						field_buffer[xx + spacex * (yy + spacey * zz)] = 1.0;
					}
					distance_buffer[xx + spacex * (yy + spacey * zz)] = dd_sqrt;
					if (mCurrentDistance < dd_sqrt) mCurrentDistance = dd_sqrt;
				}
				else if (dd < (initialRadius + 1) * (initialRadius + 1)) {
#ifdef USE_VECTOR_AS_SET
					mActivePoints[idx].push_back(point_to_representation(xx, yy, zz));
#else
					mActivePoints[idx].insert(point_to_representation(xx, yy, zz));
#endif
					BUF_IDX3D(field_buffer, spacex, spacey, spacez, xx, yy, zz) = initialRadius - sqrt(dd);
				}
			}
		}
	}
}

int AreaEikonal::GetSmallestPlaneSize() {
	vector<uint32_t> size = mPhaseFieldMap[0].GetSize();
	unsigned int xs(size[0]), ys(size[1]), zs(size[2]);
	return xs * ys < xs* zs ? xs * ys :
		xs * zs < ys* zs ? xs * zs : ys * zs;
}

void AreaEikonal::InitializeMeetingPlaneFromInitPoints(){
    mMeetingPlaneCenter = (mStartPoint+mEndPoint)/2;
    mPlaneFinalized = true;
    mMeetingPlaneNormal = mEndPoint - mStartPoint;
    mMeetingPlaneNormal /= mMeetingPlaneNormal.Norm();
    mMeetingPlaneNormal = StandardizeVector(mMeetingPlaneNormal);
	mMeetingPlaneOffset = -(mMeetingPlaneCenter*mMeetingPlaneNormal).Sum();
}

void AreaEikonal::Initialize(const sitk::Image& phi, Vec3<double>& startPoint, Vec3<double>& endPoint) {
	_PROFILING;
	mMeetingPointsCount = 0;
    mStartPoint = startPoint;
	mEndPoint = endPoint;
    if(mUseMeetingPoints){
        mMeetingPlaneCenter = Vec3<double>({ -1, -1, -1 });
        mMeetingPlaneNormal = Vec3<double>({ -2, -2, -2 });
        mPlaneFinalized = false;
	    mMeetingPlaneOffset = 1e11;
    } else {
        InitializeMeetingPlaneFromInitPoints();
    }
	mIterationCount = 0;
	mPhiMap = phi;
	InitializeContainers(phi);
	mActivePoints[0].clear();
	mActivePoints[1].clear();
	mInactivePoints[0].clear();
	mInactivePoints[1].clear();
	mVelocities[0].clear();
	mVelocities[1].clear();
	mMeetingPoints.clear();
	//g_cyc = 0;
	
	Vec3<double> point_diff(startPoint - endPoint);
	double point_dist = sqrt((point_diff*point_diff).Sum());
	double init_dist = sqrt(point_dist);
	InitializeDistanceMap(0, startPoint, init_dist);
	InitializeDistanceMap(1, endPoint, init_dist);
	mMeetingPointsMin = GetSmallestPlaneSize() / 4;
}

AreaEikonal AreaEikonal::Rotate(vector<double>& rotation_matrix, bool inverse) const {
	_PROFILING;
	AreaEikonal rotated;
	rotated.mMeetingPointsCount = 0;
	rotated.mMeetingPlaneCenter = Vec3<double>({ -1, -1, -1 });
	rotated.mIterationCount = 0;
	rotated.mPlaneFinalized = false;
	rotated.mMeetingPlaneOffset = 1e11;
	rotated.mUsesCorrection = mUsesCorrection;
	resample_img<sitk::sitkNearestNeighbor>(mMeetingPointsMap, rotated.mMeetingPointsMap, rotation_matrix, inverse);
	resample_img(mPhaseFieldMap[0], rotated.mPhaseFieldMap[0], rotation_matrix, inverse);
	neg_to_minus1(rotated.mPhaseFieldMap[0]);
	resample_img(mPhaseFieldMap[1], rotated.mPhaseFieldMap[1], rotation_matrix, inverse);
	neg_to_minus1(rotated.mPhaseFieldMap[1]);
#ifdef ROTATED_DISTMAP_W_FM
    sitk::Image mask = sitk::Image(mMeetingPointsMap.GetSize(), mDistanceMap[0].GetPixelID());
	sitk::Image rotated_mask = sitk::Image();
	resample_img<sitk::sitkNearestNeighbor>(mask, rotated_mask, rotation_matrix, inverse);
	rotated_mask = sitk::Cast(sitk::BinaryErode(sitk::Cast(rotated_mask+1, sitk::sitkUInt8)), mDistanceMap[0].GetPixelID());
	resample_img<sitk::sitkLinear>(mDistanceMap[0], rotated.mDistanceMap[0], rotation_matrix, inverse);
	neg_to_minus1(rotated.mDistanceMap[0]);
	sitk::Image frozen = sitk::Cast(sitk::Greater(rotated_mask * rotated.mDistanceMap[0], 0), sitk::sitkUInt8);
	vector<POINT3D> boundary_points = GetBoundaryPixels<sitk::sitkUInt8>(frozen, 1);
	FMSigned::build(0, rotated.mDistanceMap[0], frozen, rotated.mDistanceMap[0], boundary_points);
	resample_img<sitk::sitkLinear>(mDistanceMap[1], rotated.mDistanceMap[1], rotation_matrix, inverse);
	neg_to_minus1(rotated.mDistanceMap[1]);
	frozen = sitk::Cast(sitk::Greater(rotated_mask * rotated.mDistanceMap[1], 0), sitk::sitkUInt8);
	boundary_points = GetBoundaryPixels<sitk::sitkUInt8>(frozen, 1);
	FMSigned::build(1, rotated.mDistanceMap[1], frozen, rotated.mDistanceMap[1], boundary_points);
#else
    resample_img<sitk::sitkLinear>(mDistanceMap[0], rotated.mDistanceMap[0], rotation_matrix, inverse, true);
    neg_to_minus1(rotated.mDistanceMap[0]);
	resample_img<sitk::sitkLinear>(mDistanceMap[1], rotated.mDistanceMap[1], rotation_matrix, inverse, true);
	neg_to_minus1(rotated.mDistanceMap[1]);
#endif
	rotated.mSampleImage = resample_img(mCombinedDistanceMap, rotated.mCombinedDistanceMap, rotation_matrix, inverse, true);
	neg_to_minus1(rotated.mCombinedDistanceMap);

	resample_img(mPhiMap, rotated.mPhiMap, rotation_matrix, inverse, true);

	rotated.mCurrentDistance = mCurrentDistance;

	rotated.mCurvatureMap[0] = sitk::Image(rotated.mPhiMap.GetSize(), sitk::sitkFloat64);
	rotated.mCurvatureMap[1] = sitk::Image(rotated.mPhiMap.GetSize(), sitk::sitkFloat64);
#ifdef DEBUG_CURVATURE
	
	rotated.m_new_update[0] = sitk::Image(rotated.mPhiMap.GetSize(), sitk::sitkUInt8);
	rotated.m_new_update[1] = sitk::Image(rotated.mPhiMap.GetSize(), sitk::sitkUInt8);
#endif
	vector<unsigned> new_size = rotated.mSampleImage.GetSize();
	rotated.mMeetingPointsMin = rotated.GetSmallestPlaneSize()/4;

	for (int ii = 0; ii < 2; ii++) {
		for (auto it = mActivePoints[ii].begin(); it != mActivePoints[ii].end(); it++) {
			auto [xx, yy, zz] = representation_to_point<double>(*it);
			vector<double> pos({ xx, yy, zz });
			pos = mPhiMap.TransformContinuousIndexToPhysicalPoint(pos);
			pos = rotated.mSampleImage.TransformPhysicalPointToContinuousIndex(pos);
#ifdef USE_VECTOR_AS_SET
			rotated.mActivePoints[ii].push_back(point_to_representation((int)pos[0], (int)pos[1], (int)pos[2]));
#else
			mActivePoints[ii].insert(point_to_representation(pos[0], pos[1], pos[2]));
#endif
		}
		for (auto it = mInactivePoints[ii].begin(); it != mInactivePoints[ii].end(); it++) {
			auto [xx, yy, zz] = representation_to_point<double>(*it);
			vector<double> pos({ xx, yy, zz });
			pos = mPhiMap.TransformContinuousIndexToPhysicalPoint(pos);
			pos = rotated.mSampleImage.TransformPhysicalPointToContinuousIndex(pos);
#ifdef USE_VECTOR_AS_SET
			rotated.mActivePoints[ii].push_back(point_to_representation((int)pos[0], (int)pos[1], (int)pos[2]));
#else
			mActivePoints[ii].insert(point_to_representation(pos[0], pos[1], pos[2]));
#endif
		}
	}
	return rotated;
}

std::vector<double> AreaEikonal::TransformPhysicalPointToContinuousIndex(const std::vector<double>& point) const {
	return mSampleImage.TransformPhysicalPointToContinuousIndex(point);
}

std::vector<double> AreaEikonal::TransformContinuousIndexToPhysicalPoint(const std::vector<double>& point) const {
	return mSampleImage.TransformContinuousIndexToPhysicalPoint(point);
}
// The number of threads for a parallel loop over `work` points: the default team, but not more
// than one thread for every kPointsPerThread points, so that the loops over the few points of the
// first iterations do not wake the whole team.
static const size_t kPointsPerThread = 256;
static int TeamSize(size_t work) {
	size_t by_work = work / kPointsPerThread;
	return (int)std::max<size_t>(1, std::min<size_t>(std::max(1, omp_get_max_threads()), by_work));
}

int current_iteration = 0;

// The steps of an iteration below go over the points of both sides in one loop: the points of side
// 0 first, then those of side 1 (index k; side s = k >= n0, point k - n0 of side 1). All of the
// threads of the team work on both sides at the same time, without nested parallel regions.

void AreaEikonal::UpdateCurvatures() {
	const vector<unsigned>& size = mPhaseFieldMap[0].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	const int n0 = (int)mActivePoints[0].size();
	const int total = n0 + (int)mActivePoints[1].size();
	const int team = TeamSize(total);
	_SECTION_BEGIN(sdist_timer, "evolve.1a signed distance");
	const sitk::Image* inputs[2] = { &mPhaseFieldMap[0], &mPhaseFieldMap[1] };
	const vector<POINT3D>* bands[2] = { &mActivePoints[0], &mActivePoints[1] };
	const sitk::Image* sdist_maps[2];
	FMSigned::build_for_neighbors_both(inputs, bands, 2, 1, team, sdist_maps);
	_SECTION_END(sdist_timer);
	const double* sdist_buffer[2];
	double* curvature_buffer[2];
	_SECTION("evolve.1b curvature");
	for (int s = 0; s < 2; s++) {
		sdist_buffer[s] = sdist_maps[s]->GetBufferAsDouble();
		curvature_buffer[s] = mCurvatureMap[s].GetBufferAsDouble();
	}
	// Only the points of the band are written below, and the map is zero everywhere else. So
	// instead of clearing the whole volume, the points the last call wrote are cleared: the map
	// ends up the same. The old and the new list of points live side by side until both loops
	// are done.
	const int c0 = (int)mCurvatureBand[0].size();
	const int clear_total = c0 + (int)mCurvatureBand[1].size();
OMP_PARALLEL_NUM_THREADS(team)
	{
OMP_FOR
		for (int k = 0; k < clear_total; k++) {
			const int s = k >= c0 ? 1 : 0;
			POINT3D pn = mCurvatureBand[s][s ? k - c0 : k];
			auto [xx, yy, zz] = representation_to_point<int>(pn);
			if (zz < 0 || zz >= zs || yy < 0 || yy >= ys || xx < 0 || xx >= xs) continue;
			BUF_IDX3D(curvature_buffer[s], xs, ys, zs, xx, yy, zz) = 0;
		}
OMP_FOR_NOWAIT
		for (int k = 0; k < total; k++) {
			const int s = k >= n0 ? 1 : 0;
			POINT3D pn = mActivePoints[s][s ? k - n0 : k];
			auto [xx, yy, zz] = representation_to_point<int>(pn);
			if (zz < 0 || zz >= zs) continue;
			if (yy < 0 || yy >= ys) continue;
			if (xx < 0 || xx >= xs) continue;

			int zp = zz < zs - 1 ? zz + 1: zs - 1;
			int zm = zz > 0 ? zz - 1 : 0;
			int yp = yy < ys - 1 ? yy + 1 : ys - 1;
			int ym = yy > 0 ? yy - 1 : 0;
			int xp = xx < xs - 1 ? xx + 1 : xs - 1;
			int xm = xx > 0 ? xx - 1 : 0;
			const double* sdist = sdist_buffer[s];
			realnum sumcur =
				BUF_IDX3D(sdist, xs, ys, zs, xp, yy, zz) +
				BUF_IDX3D(sdist, xs, ys, zs, xm, yy, zz) +
				BUF_IDX3D(sdist, xs, ys, zs, xx, yp, zz) +
				BUF_IDX3D(sdist, xs, ys, zs, xx, ym, zz) +
				BUF_IDX3D(sdist, xs, ys, zs, xx, yy, zp) +
				BUF_IDX3D(sdist, xs, ys, zs, xx, yy, zm) -
				6 * BUF_IDX3D(sdist, xs, ys, zs, xx, yy, zz);
			//if (sumcur < -2) sumcur == -2;
			BUF_IDX3D(curvature_buffer[s], xs, ys, zs, xx, yy, zz) = sumcur;
		}
	}
	for (int s = 0; s < 2; s++)
		mCurvatureBand[s].assign(mActivePoints[s].begin(), mActivePoints[s].end());
	OnCurvatureUpdated(0);
	OnCurvatureUpdated(1);
}

// The idx-th smallest of the n_top smallest values of S, which is the idx-th smallest value
// overall: nth_element finds it without sorting them all. (This used to sort per-thread lists and
// merge their first n_top values, a merge that stopped at the first value not below DBL_MAX; such
// values, and NaN, are left out here too.) DBL_MAX if there is none.
static double SelectMinS(vector<double>& s_values, int n_outliers, int n_top) {
	if (s_values.empty())
		return DBL_MAX;
	auto valid_end = std::partition(s_values.begin(), s_values.end(), [](double v) { return v < DBL_MAX; });
	size_t valid_count = valid_end - s_values.begin();
	if (valid_count == 0)
		return DBL_MAX;
	size_t top_count = min((size_t)n_top, valid_count);
	size_t idx = min((size_t)max(200, n_outliers / 20), top_count - 1);
	std::nth_element(s_values.begin(), s_values.begin() + idx, valid_end);
	return s_values[idx];
}

void AreaEikonal::GetMinS(double min_s[2]) const {
	_SECTION("evolve.2 minS");
	int n_top = 5000;
	const double* phi_buffer = mPhiMap.GetBufferAsDouble();
	const double* field_buffer[2];
	const double* curvature_buffer[2];
	const vector<unsigned>& size = mPhaseFieldMap[0].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	const int n0 = (int)mActivePoints[0].size();
	const int total = n0 + (int)mActivePoints[1].size();
	// S of every point of the bands, in band order. Points outside the volume get DBL_MAX, which
	// the selection below leaves out.
	vector<double> s_values[2] = { vector<double>(mActivePoints[0].size(), DBL_MAX), vector<double>(mActivePoints[1].size(), DBL_MAX) };
	double* s_data[2] = { s_values[0].data(), s_values[1].data() };
	int n_outliers[2] = { 0, 0 };
	for (int s = 0; s < 2; s++) {
		field_buffer[s] = mPhaseFieldMap[s].GetBufferAsDouble();
		curvature_buffer[s] = mCurvatureMap[s].GetBufferAsDouble();
	}

OMP_PARALLEL_NUM_THREADS(TeamSize(total))
	{
		int outliers[2] = { 0, 0 };
OMP_FOR_NOWAIT
		for (int k = 0; k < total; k++) {
			const int s = k >= n0 ? 1 : 0;
			const int j = s ? k - n0 : k;
			POINT3D pn = mActivePoints[s][j];
			auto [xx, yy, zz] = representation_to_point<int>(pn);
			if (zz < 0 || zz >= zs) continue;
			if (yy < 0 || yy >= ys) continue;
			if (xx < 0 || xx >= xs) continue;

			int zp = zz < zs - 1 ? zz + 1 : zs - 1;
			int zm = zz > 0 ? zz - 1 : 0;
			int yp = yy < ys - 1 ? yy + 1 : ys - 1;
			int ym = yy > 0 ? yy - 1 : 0;
			int xp = xx < xs - 1 ? xx + 1 : xs - 1;
			int xm = xx > 0 ? xx - 1 : 0;

			const double* field = field_buffer[s];
			realnum phi = phi_buffer[xx + xs * (yy + ys * zz)];
			if (phi < 1e-33) phi = 1e-33;
			double sumcur = BUF_IDX3D(curvature_buffer[s], xs, ys, zs, xx, yy, zz);
			realnum&& xn = field[xp + xs * (yy + ys * zz)] - field[xm + xs * (yy + ys * zz)];
			realnum&& yn = field[xx + xs * (yp + ys * zz)] - field[xx + xs * (ym + ys * zz)];
			realnum&& zn = field[xx + xs * (yy + ys * zp)] - field[xx + xs * (yy + ys * zm)];
			realnum gradlen = sqrt(xn * xn + yn * yn + zn * zn + 1e-99);

			if (is_outlier(sumcur, 2/gradlen))
			{
				outliers[s]++;
			}
			s_data[s][j] = calculate_S(sumcur, phi, 2/gradlen);
		}
OMP_CRITICAL_NO_TAG
		{
			n_outliers[0] += outliers[0];
			n_outliers[1] += outliers[1];
		}
	}
	// The selections for the two sides are independent: one thread for each.
OMP_PARALLEL_SECTIONS_NUM_THREADS(2)
	{
OMP_SECTION
		min_s[0] = SelectMinS(s_values[0], n_outliers[0], n_top);
OMP_SECTION
		min_s[1] = SelectMinS(s_values[1], n_outliers[1], n_top);
	}
}

void AreaEikonal::UpdateVelo(double S, double max_velocity[2]) {
	_SECTION("evolve.3 velocity");
	const vector<unsigned>& size = mPhaseFieldMap[0].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	// Both sides read mPhiMap. The non-const GetBufferAsDouble() makes the image unique first
	// (copy-on-write): if its data is still shared with another copy, such as the caller's, it
	// clones the data and swaps its own internal pointer, which races with another thread doing
	// the same or reading the buffer. Going through a const reference never modifies the image.
	const double* phi_buffer = std::as_const(mPhiMap).GetBufferAsDouble();
	const double* field_buffer[2];
	const double* curvature_buffer[2];
	for (int s = 0; s < 2; s++) {
		field_buffer[s] = mPhaseFieldMap[s].GetBufferAsDouble();
		curvature_buffer[s] = mCurvatureMap[s].GetBufferAsDouble();
		mVelocities[s].resize(mActivePoints[s].size());
	}
	const int n0 = (int)mActivePoints[0].size();
	const int total = n0 + (int)mActivePoints[1].size();
	max_velocity[0] = max_velocity[1] = 0;

	// Inhomogeneous
	// Calculate the velocity at every point (velo[zz][yy][xx]), and the maximum velocity (maxv)
OMP_PARALLEL_NUM_THREADS(TeamSize(total))
	{
		double thread_max[2] = { 0, 0 };
OMP_FOR_NOWAIT
		for (int k = 0; k < total; k++) {
			const int s = k >= n0 ? 1 : 0;
			const int j = s ? k - n0 : k;
			POINT3D pn = mActivePoints[s][j];
			auto [xx, yy, zz] = representation_to_point<int>(pn);
			if (zz < 0 || zz >= zs) continue;
			if (yy < 0 || yy >= ys) continue;
			if (xx < 0 || xx >= xs) continue;

			int zp = zz < zs - 1 ? zz + 1 : zs - 1;
			int zm = zz > 0 ? zz - 1 : 0;
			int yp = yy < ys - 1 ? yy + 1 : ys - 1;
			int ym = yy > 0 ? yy - 1 : 0;
			int xp = xx < xs - 1 ? xx + 1 : xs - 1;
			int xm = xx > 0 ? xx - 1 : 0;

			const double* field = field_buffer[s];
			realnum&& xn = field[xp + xs * (yy + ys * zz)] - field[xm + xs * (yy + ys * zz)];
			realnum&& yn = field[xx + xs * (yp + ys * zz)] - field[xx + xs * (ym + ys * zz)];
			realnum&& zn = field[xx + xs * (yy + ys * zp)] - field[xx + xs * (yy + ys * zm)];
			realnum gradlen = sqrt(xn * xn + yn * yn + zn * zn + 1e-99);

			realnum eikon;

			realnum phi = phi_buffer[xx + xs * (yy + ys * zz)];
			if (phi < 1e-33) phi = 1e-33;
			if (mUsesCorrection) {
				realnum sumcur = BUF_IDX3D(curvature_buffer[s], xs, ys, zs, xx, yy, zz);
				eikon = calculate_speed(sumcur, phi, S);
			}
			else {
				eikon = 1 / phi;
			}
			eikon *= gradlen; // normalize
			if (eikon < 1e-11) eikon = 1e-11;

			mVelocities[s][j] = make_pair(pn, eikon);
			if (eikon > thread_max[s]) thread_max[s] = eikon;
		}
OMP_CRITICAL(maxS_selection)
		{
			for (int s = 0; s < 2; s++)
				if (thread_max[s] > max_velocity[s]) max_velocity[s] = thread_max[s];
		}
	}
}

void AreaEikonal::UpdateFields(double maxv) {
	_SECTION("evolve.4 field");
	const vector<unsigned>& size = mPhaseFieldMap[0].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	double* field_buffer[2] = { mPhaseFieldMap[0].GetBufferAsDouble(), mPhaseFieldMap[1].GetBufferAsDouble() };
	const int n0 = (int)mVelocities[0].size();
	const int total = n0 + (int)mVelocities[1].size();
OMP_PARALLEL_FOR_NUM_THREADS(TeamSize(total))
	for (int k = 0; k < total; k++) {
		const int s = k >= n0 ? 1 : 0;
		const int j = s ? k - n0 : k;
		POINT3D pn = mVelocities[s][j].first;
		double velo = mVelocities[s][j].second;

		auto [xx, yy, zz] = representation_to_point<int>(pn);
		BUF_IDX3D(field_buffer[s], xs, ys, zs, xx, yy, zz) += velo;
		if (BUF_IDX3D(field_buffer[s], xs, ys, zs, xx, yy, zz) >= MAX_FIELD_VAL) {
			mActivePoints[s][j] = REMOVABLE_POINT;
			if (BUF_IDX3D(field_buffer[s], xs, ys, zs, xx, yy, zz) >= 1.25)
				BUF_IDX3D(field_buffer[s], xs, ys, zs, xx, yy, zz) = 1.25;
		}
	}
	OnFieldUpdated(0);
	OnFieldUpdated(1);
}

void AreaEikonal::UpdateDistances(double current_distance) {
	_SECTION_BEGIN(loop_timer, "evolve.5a distance");
	const vector<unsigned>& size = mPhaseFieldMap[0].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	const Vec3<double> origin_point[2] = { mStartPoint, mEndPoint };
	double* field_buffer[2];
	double* distance_buffer[2];
	for (int s = 0; s < 2; s++) {
		field_buffer[s] = mPhaseFieldMap[s].GetBufferAsDouble();
		distance_buffer[s] = mDistanceMap[s].GetBufferAsDouble();
	}
	const int n0 = (int)mVelocities[0].size();
	const int total = n0 + (int)mVelocities[1].size();
	const int team = TeamSize(total);
	// The new points each thread finds, for each side, appended to the bands in thread order below.
	vector<vector<POINT3D>> found(2 * team);
OMP_PARALLEL_NUM_THREADS(team)
	{
OMP_FOR_NOWAIT
		for (int k = 0; k < total; k++) {
			const int s = k >= n0 ? 1 : 0;
			vector<POINT3D>& temp = found[2 * omp_get_thread_num() + s];
			POINT3D pn = mVelocities[s][s ? k - n0 : k].first;
			auto [xx, yy, zz] = representation_to_point<int>(pn);
			if (field_buffer[s][xx + xs * (yy + ys * zz)] > 0 && distance_buffer[s][xx + xs * (yy + ys * zz)] < -0.5f) {
				distance_buffer[s][xx + xs * (yy + ys * zz)] = current_distance;
				for (int dk = -1; dk <= 1; dk++) {
					if (zz + dk >= zs || zz + dk < 0) continue;
					for (int dj = -1; dj <= 1; dj++) {
						if (yy + dj >= ys || yy + dj < 0) continue;
						for (int di = -1; di <= 1; di++) {
							if (xx + di >= xs || xx + di < 0) continue;
							if ((di + 1) % 2 + (dj + 1) % 2 + (dk + 1) % 2 == 2) {
								if (BUF_IDX3D(field_buffer[s], xs, ys, zs, xx + di, yy + dj, zz + dk) < MAX_FIELD_VAL
									&& (mMeetingPlaneCenter.x() < -0.5 || -sgn((mMeetingPlaneNormal * origin_point[s]).Sum() + mMeetingPlaneOffset) * ((mMeetingPlaneNormal * Vec3<double>({ (double)xx + di, (double)yy + dj, (double)zz + dk })).Sum() + mMeetingPlaneOffset) < MAX_DIST_FROM_PLANE))
								{
									if (xx + di > 0 && xx + di < xs - 1 && yy + dj > 0 && yy + dj < ys - 1 && zz + dk > 0 && zz + dk < zs - 1) {
										temp.push_back(point_to_representation(xx + di, yy + dj, zz + dk));
									}
								} // end if
							} // end if
						} // end for di
					} // end for dj
				} // end for dk
			}
		}
	}
	for (int t = 0; t < team; t++)
		for (int s = 0; s < 2; s++)
			mActivePoints[s].insert(mActivePoints[s].end(), found[2 * t + s].begin(), found[2 * t + s].end());

	_SECTION_END(loop_timer);
	_SECTION("evolve.5b remove duplicates");
	mVelocities[0].clear();
	mVelocities[1].clear();
	_RECORD_VALUE("evolve.5b size before removing duplicates", mActivePoints[0].size() + mActivePoints[1].size());
	// The two sides are independent: one thread for each.
OMP_PARALLEL_SECTIONS_NUM_THREADS(2)
	{
OMP_SECTION
		RemoveDuplicatesFromBand(0);
OMP_SECTION
		RemoveDuplicatesFromBand(1);
	}
}

void AreaEikonal::RemoveDuplicatesFromBand(int idx) {
	// The band now holds the points UpdateField marked as REMOVABLE_POINT, and new neighbours,
	// many of them more than once or already in the band. Drop both. The order of the band does
	// not matter to any step of the evolution (each treats its points independently; the only
	// reductions are a maximum and an order statistic), so instead of sorting it, keep the first
	// occurrence of each point, found with a per-voxel flag that is cleared again afterwards.
	// The band ends up with the same points as a sort and unique would leave.
	auto& act_set = mActivePoints[idx];
	vector<unsigned> size = mPhaseFieldMap[idx].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	size_t voxel_count = (size_t)xs * ys * zs;
	for (POINT3D pn : act_set) {
		if (pn == REMOVABLE_POINT) continue;
		auto [xx, yy, zz] = representation_to_point<int>(pn);
		if (xx < 0 || xx >= xs || yy < 0 || yy >= ys || zz < 0 || zz >= zs) {
			// Not a voxel (a band rotated from another volume can have such points): fall back
			// to sorting, which keeps them as they were kept before.
			sort(act_set.begin(), act_set.end());
			act_set.erase(unique(act_set.begin(), act_set.end()), act_set.end());
			if (!act_set.empty() && act_set.back() == REMOVABLE_POINT)
				act_set.pop_back();
			return;
		}
	}
	vector<unsigned char>& in_band = mInBand[idx];
	if (in_band.size() != voxel_count)
		in_band.assign(voxel_count, 0);
	size_t kept = 0;
	for (POINT3D pn : act_set) {
		if (pn == REMOVABLE_POINT) continue;
		auto [xx, yy, zz] = representation_to_point<int>(pn);
		unsigned char& seen = in_band[xx + (size_t)xs * (yy + (size_t)ys * zz)];
		if (seen) continue;
		seen = 1;
		act_set[kept++] = pn;
	}
	act_set.resize(kept);
	for (POINT3D pn : act_set) {
		auto [xx, yy, zz] = representation_to_point<int>(pn);
		in_band[xx + (size_t)xs * (yy + (size_t)ys * zz)] = 0;
	}
}

bool AreaEikonal::IsDone() const {
	_PROFILING;
	_SECTION("iter.IsDone");
	if (mActivePoints[0].empty() && mActivePoints[1].empty())
		return true;
	vector<unsigned> size = mPhaseFieldMap[0].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	const double* dist0_buffer = mDistanceMap[0].GetBufferAsDouble();
	const double* dist1_buffer = mDistanceMap[1].GetBufferAsDouble();
	// Done when no voxel is left that neither side has reached. Look for one, starting at the
	// one found last time and wrapping around: that one is usually still unreached, so most
	// calls look at a single voxel instead of the whole volume. The answer is the same as a
	// scan from the start, since every voxel is looked at before deciding that it is done.
	int voxel_count = zs * ys * xs;
	int start = mUnreachedVoxel < voxel_count ? mUnreachedVoxel : 0;
	for (int offset = 0; offset < voxel_count; offset++) {
		int zyx = start + offset < voxel_count ? start + offset : start + offset - voxel_count;
		if (dist0_buffer[zyx] < 0 && dist1_buffer[zyx] < 0) {
			mUnreachedVoxel = zyx;
			return false;
		}
	}
	return true;
}
void AreaEikonal::Iterate()
{
	_PROFILING;
	_SECTION("iter.Iterate");
	_RECORD_VALUE("band.side0", mActivePoints[0].size());
	_RECORD_VALUE("band.side1", mActivePoints[1].size());
	// mAreaEikonal.m_field[i]
	if (mUseMeetingPoints && !mPlaneFinalized && ++mIterationCount > REFRESH_ITERS) {
		_SECTION("iter.UpdateMeetingPlane");
		mIterationCount = 0;
		UpdateMeetingPlane();
	}
	//Calculate the value of capital S
	UpdateCurvatures();
	double min_capital_s_vals[2];
	GetMinS(min_capital_s_vals);
	realnum min_capital_s = min_capital_s_vals[0] < min_capital_s_vals[1] ? min_capital_s_vals[0] : min_capital_s_vals[1];
	if (min_capital_s < 1e-11) min_capital_s = 1e-11;

	//Calculate the velocity along the front
	double max_velocities[2];
	UpdateVelo(min_capital_s, max_velocities);
	double max_velocity = max_velocities[0] > max_velocities[1] ? max_velocities[0] : max_velocities[1];
	realnum mdatspedmax = 2;
	max_velocity = mdatspedmax/max_velocity;
	mCurrentDistance += max_velocity;

	// update phase field and distance values
	UpdateFields(max_velocity);
	UpdateDistances(mCurrentDistance);
}
void AreaEikonal::Calculate(const sitk::Image& image, Vec3<double>& point1, Vec3<double>& point2) {
	Initialize(image, point1, point2);
	Calculate();
}
void AreaEikonal::Calculate() {
	_SECTION("eik.Calculate");
	int iterations = 0;
	while (!IsDone()) {
		Iterate();
		iterations++;
	}
	_RECORD_VALUE("eik.iterations per Calculate", iterations);
}
void AreaEikonal::CombineDistance(int slice, double p0_x, double p1_x) {
	_PROFILING;
	vector<unsigned> size = mPhiMap.GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	mCombinedDistanceMap = sitk::Image({ (unsigned) xs, (unsigned)ys, (unsigned)zs }, sitk::sitkFloat64) - 1;
	double* dist0_buffer = mDistanceMap[0].GetBufferAsDouble();
	double* dist1_buffer = mDistanceMap[1].GetBufferAsDouble();
	double* combined_dist_buffer = mCombinedDistanceMap.GetBufferAsDouble();
OMP_PARALLEL_FOR
	for (int zyx = 0; zyx < zs*ys*xs; zyx++) {
		int xx = zyx % xs;
		int xx_sgn = sgn(slice - xx);
		int p0_sgn = sgn(slice - p0_x);
		int p1_sgn = sgn(slice - p1_x);
		if (abs(slice - xx) < 0.5) {
			combined_dist_buffer[zyx] = (dist0_buffer[zyx + p0_sgn] + dist1_buffer[zyx + p1_sgn]);
		}
		else if (p0_sgn == xx_sgn) {
			combined_dist_buffer[zyx] = dist0_buffer[zyx];
		}
		else if (p1_sgn == xx_sgn) {
			combined_dist_buffer[zyx] = dist1_buffer[zyx];
		}
	}
}

void AreaEikonal::SetUsesCorrection(bool useCorrection) {
	mUsesCorrection = useCorrection;
}

void AreaEikonal::SetUsingMeetingPoints(bool useMeetingPoints){
    mUseMeetingPoints = useMeetingPoints;
}

bool AreaEikonal::IsUsingMeetingPoints(){
    return mUseMeetingPoints;
}

Vec3<double> AreaEikonal::GetMeetingPlaneCenter() const {
	return mMeetingPlaneCenter;
}

void AreaEikonal::SetMeetingPlaneCenter(Vec3<double> center){
    mMeetingPlaneCenter = center;
    if(mMeetingPlaneNormal.x() < -1){
        return;
    }
    UpdateMeetingPlaneInfo();
}

Vec3<double> AreaEikonal::GetMeetingPlaneNormal() const {
	return mMeetingPlaneNormal;
}

void AreaEikonal::SetMeetingPlaneNormal(Vec3<double> normal){
    mMeetingPlaneNormal = normal;
    if(mMeetingPlaneCenter.x() < 0){
        return;
    }
    UpdateMeetingPlaneInfo();
}

double AreaEikonal::GetMeetingPlaneOffset() const {
	return mMeetingPlaneOffset;
}

void AreaEikonal::UpdateMeetingPlaneOffset(){
    mMeetingPlaneOffset = -(mMeetingPlaneCenter*mMeetingPlaneNormal).Sum();
}

void AreaEikonal::UpdateRotationMatrix(){
    vector<double> meeting_plane_normal(mMeetingPlaneNormal);
    mRotationMatrix = rotation_matrix_from_vectors<double>(vector<double>({ 1, 0, 0 }), meeting_plane_normal);
}

void AreaEikonal::UpdateMeetingPlaneInfo(){
    UpdateMeetingPlaneOffset();
    UpdateRotationMatrix();
}

const sitk::Image& AreaEikonal::GetPhiMap() const {
	return mPhiMap;
}

const sitk::Image& AreaEikonal::GetMeetingPointsMap() const {
	return mMeetingPointsMap;
}

const sitk::Image& AreaEikonal::GetSampleImage() const
{
	return mSampleImage;
}

double AreaEikonal::GetCurrentDistance() const {
	return mCurrentDistance;
}

vector<double> AreaEikonal::GetRotationMatrix() {
    if(mRotationMatrix.size() == 0){
        cout << "Rotation matrix missing!!!" << endl;
    }
    return mRotationMatrix;
}

const sitk::Image& AreaEikonal::GetDistanceMap(int i) const {
	return mDistanceMap[i];
}

const sitk::Image& AreaEikonal::GetCombinedDistanceMap() const {
	return mCombinedDistanceMap;
}

const sitk::Image& AreaEikonal::GetPhaseFieldMap(int i) const {
	return mPhaseFieldMap[i];
}

vector<Vec3<int>> AreaEikonal::ResolvePath(Vec3<int> point, int idx) const {
	return ::ResolvePath(point, mDistanceMap[idx]);
}