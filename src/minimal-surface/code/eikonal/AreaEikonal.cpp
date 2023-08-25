#include "stdafx.h"
#include "AreaEikonal.h"
#include "Utils.h"
#include "math.h"
#include <omp.h>
#include "SimpleITK.h"
#include <vector>
#include <sitkImageOperators.h>
#include "FMSigned.h"
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
OMP_PARALLEL
	{
		vector<POINT3D> temp;
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
OMP_CRITICAL(fill_meeting_plane)
		mMeetingPoints.insert(temp.begin(), temp.end());
	}
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
	mMeetingPlaneOffset = -(mMeetingPlaneCenter*mMeetingPlaneNormal).Sum();

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
#ifdef DEBUG_CURVATURE

	m_new_update[0] = sitk::Image(size, sitk::sitkUInt8);
	m_new_update[1] = sitk::Image(size, sitk::sitkUInt8);
#endif
}

void AreaEikonal::InitializeDistanceMap(int idx, Vec3<double> center, int initialRadius) {
	vector<uint32_t> size = mPhaseFieldMap[0].GetSize();
	unsigned int spacex(size[0]), spacey(size[1]), spacez(size[2]); 
	int hs = initialRadius + 2;
	mCurrentDistance = 0;
	double* field_buffer = mPhaseFieldMap[idx].GetBufferAsDouble();
	double* distance_buffer = mDistanceMap[idx].GetBufferAsDouble();
	for (int zz = center.z() - hs; zz < center.z() + hs; ++zz) {
		if (zz < 0 || zz >= spacez) continue;
		for (int yy = center.y() - hs; yy < center.y() + hs; ++yy) {
			if (yy < 0 || yy >= spacey) continue;
			for (int xx = center.x() - hs; xx < center.x() + hs; ++xx) {
				if (xx < 0 || xx >= spacex) continue;
				int dx = xx - center.x(), dy = yy - center.y();
				int dz = zz - center.z();
				int dd = dx * dx + dy * dy + dz * dz;
				if ((int)dd <= initialRadius * initialRadius) {
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
	mMeetingPlaneOffset = -(mMeetingPlaneCenter*mMeetingPlaneNormal).Sum();
}

void AreaEikonal::Initialize(const sitk::Image& image, Vec3<double>& startPoint, Vec3<double>& endPoint, double beta, double alpha) {
	_PROFILING;
	mMeetingPointsCount = 0;
    mStartPoint = startPoint;
	mEndPoint = endPoint;
#ifdef MEETING_PLANE_FROM_INIT_POINTS
    InitializeMeetingPlaneFromInitPoints();
#else
	mMeetingPlaneCenter = Vec3<double>({ -1, -1, -1 });
    mPlaneFinalized = false;
	mMeetingPlaneOffset = 1e11;
#endif
	mIterationCount = 0;
	mPhiMap = CalculatePhi(image, beta, alpha);
	InitializeContainers(image);
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
	sitk::Image mask = sitk::Image(mMeetingPointsMap.GetSize(), mDistanceMap[0].GetPixelID());
	sitk::Image rotated_mask = sitk::Image();
	resample_img<sitk::sitkNearestNeighbor>(mask, rotated_mask, rotation_matrix, inverse);
	rotated_mask = sitk::Cast(sitk::BinaryErode(sitk::Cast(rotated_mask+1, sitk::sitkUInt8)), mDistanceMap[0].GetPixelID());
	resample_img(mPhaseFieldMap[0], rotated.mPhaseFieldMap[0], rotation_matrix, inverse);
	neg_to_minus1(rotated.mPhaseFieldMap[0]);
	resample_img(mPhaseFieldMap[1], rotated.mPhaseFieldMap[1], rotation_matrix, inverse);
	neg_to_minus1(rotated.mPhaseFieldMap[1]);
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
			rotated.mActivePoints[ii].push_back(point_to_representation(pos[0], pos[1], pos[2]));
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
			rotated.mActivePoints[ii].push_back(point_to_representation(pos[0], pos[1], pos[2]));
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
int current_iteration = 0;

void AreaEikonal::UpdateCurvature(int i) {
	static int n_threads = omp_get_max_threads();
	sitk::Image& field = mPhaseFieldMap[i];
	auto& act_set = mActivePoints[i];
	const vector<unsigned>& size = mPhaseFieldMap[i].GetSize();
	int xs = size[0], ys = size[1], zs = size[2]; 
	sitk::Image& temp_sdist = FMSigned::build_for_neighbors(i, field, act_set, 2, 1);
	double* sdist_buffer = temp_sdist.GetBufferAsDouble();
	double* curvature_buffer = mCurvatureMap[i].GetBufferAsDouble();
	memset(curvature_buffer, 0, xs * ys * zs * sizeof(double));
OMP_PARALLEL_FOR_NUM_THREADS(n_threads/2)
#ifdef USE_VECTOR_AS_SET
	for (int i = 0; i < act_set.size(); i++)
	{
		{
			POINT3D pn = act_set[i];
#else
	for (int b = 0; b < act_set.bucket_count(); b++)
	{
		for (auto bi = act_set.begin(b); bi != act_set.end(b); bi++) {
			POINT3D pn = *bi;
#endif
			{

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
				realnum sumcur =
					BUF_IDX3D(sdist_buffer, xs, ys, zs, xp, yy, zz) +
					BUF_IDX3D(sdist_buffer, xs, ys, zs, xm, yy, zz) +
					BUF_IDX3D(sdist_buffer, xs, ys, zs, xx, yp, zz) +
					BUF_IDX3D(sdist_buffer, xs, ys, zs, xx, ym, zz) +
					BUF_IDX3D(sdist_buffer, xs, ys, zs, xx, yy, zp) +
					BUF_IDX3D(sdist_buffer, xs, ys, zs, xx, yy, zm) -
					6 * BUF_IDX3D(sdist_buffer, xs, ys, zs, xx, yy, zz);
				//if (sumcur < -2) sumcur == -2;
				BUF_IDX3D(curvature_buffer, xs, ys, zs, xx, yy, zz) = sumcur;
			}
		} 
	}
}

realnum AreaEikonal::GetMinS(int i) const {

	static int n_threads = omp_get_max_threads();
	const sitk::Image& field = mPhaseFieldMap[i];
	auto& act_set = mActivePoints[i];
	int n_top = 5000;
	//int n_top = act_set.size();
	vector<double> topSs;
	topSs.reserve(n_top + 1);
	const double* field_buffer = field.GetBufferAsDouble();
	const double* phi_buffer = mPhiMap.GetBufferAsDouble();
	const double* curvature_buffer = mCurvatureMap[i].GetBufferAsDouble();
	const vector<unsigned>& size = mPhaseFieldMap[i].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	int n_outliers(0);
	realnum minS(DBL_MAX);
	if (act_set.size() == 0)
		return minS;
	vector<vector<double>> heaps(n_threads / 2);

OMP_PARALLEL_NUM_THREADS(n_threads/2)
	{
		//double temp_minS(DBL_MAX);
		vector<double>& temp = heaps[omp_get_thread_num()];
		temp.reserve(act_set.size()/(n_threads/2));
OMP_FOR_REDUCTION(+, n_outliers)
#ifdef USE_VECTOR_AS_SET
		for (int i = 0; i < act_set.size(); i++)
		{
			{
				POINT3D pn = act_set[i];
#else
		for (int b = 0; b < act_set.bucket_count(); b++)
		{
			for (auto bi = act_set.begin(b); bi != act_set.end(b); bi++) {
				POINT3D pn = *bi;
#endif
				{
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

					realnum phi = phi_buffer[xx + xs * (yy + ys * zz)];
					if (phi < 1e-33) phi = 1e-33;
					double sumcur = BUF_IDX3D(curvature_buffer, xs, ys, zs, xx, yy, zz);
					realnum&& xn = field_buffer[xp + xs * (yy + ys * zz)] - field_buffer[xm + xs * (yy + ys * zz)];
					realnum&& yn = field_buffer[xx + xs * (yp + ys * zz)] - field_buffer[xx + xs * (ym + ys * zz)];
					realnum&& zn = field_buffer[xx + xs * (yy + ys * zp)] - field_buffer[xx + xs * (yy + ys * zm)];
					realnum gradlen = sqrt(xn * xn + yn * yn + zn * zn + 1e-99);

					if (is_outlier(sumcur, 2/gradlen))
					{
						n_outliers++;
					}
					realnum S = calculate_S(sumcur, phi, 2/gradlen);
					temp.push_back(S);
				}
			} // for zz
		}
		std::sort(temp.begin(), temp.end());
	}
	vector<int> indices(heaps.size());
	for (int i = 0; i < n_top; i++) {
		double next_elem(DBL_MAX);
		int heap_idx;
		for (int h = 0; h < heaps.size(); h++) {
			int idx = indices[h];
			if (idx >= heaps[h].size()) continue;
			double top = heaps[h][idx];
			if (top < next_elem) {
				heap_idx = h;
				next_elem = top;
			}
		}
		if (next_elem == DBL_MAX) break;
		topSs.push_back(next_elem);
		indices[heap_idx]++;
	}
	if (topSs.empty()) {
		return DBL_MAX;
	}
	/*minS = 0;
	std::for_each(topSs.begin(), topSs.end(), [&minS](double val) {minS += val; });
	minS /= topSs.size();*/
	/*if (n_outliers == 0)
		minS = topSs[0];
	else
		minS = topSs[n_outliers*2];*/
	size_t idx = min((size_t)max(200, n_outliers/20), topSs.size()-1);
	minS = topSs[idx];
	return minS;
}
realnum AreaEikonal::UpdateVelo(int i, double S) {
	static int n_threads = omp_get_max_threads();
	//_PROFILING;
	sitk::Image& field = mPhaseFieldMap[i];
	double* field_buffer = field.GetBufferAsDouble();

	double* phi_buffer = mPhiMap.GetBufferAsDouble();
	

	const vector<unsigned>& size = mPhaseFieldMap[i].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];

	////////////////////////////////////////////////////////////////////
	// evolution logic
	////////////////////////////////////////////////////////////////////
	auto& act_set = mActivePoints[i];
	
	//realnum max_velocity(0);
	double* curvature_buffer = mCurvatureMap[i].GetBufferAsDouble();

	// Inhomogeneous
	// Calculate the velocity at every point (velo[zz][yy][xx]), and the maximum velocity (maxv)
	auto& changed_velo = mVelocities[i];
#ifdef USE_VECTOR_AS_SET
	changed_velo.resize(act_set.size());
#else
	changed_velo.reserve(act_set.size());
#endif
	int c = 0;

	double maxv(0);

OMP_PARALLEL_NUM_THREADS(n_threads/2)
	{
		double temp_maxv(0);
OMP_FOR
#ifdef USE_VECTOR_AS_SET
		for (int i = 0; i < act_set.size(); i++)
		{
			{
				POINT3D pn = act_set[i];
#else
		for (int b = 0; b < act_set.bucket_count(); b++)
		{
			for (auto bi = act_set.begin(b); bi != act_set.end(b); bi++) {
				POINT3D pn = *bi;
#endif
				{

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

					realnum&& xn = field_buffer[xp + xs * (yy + ys * zz)] - field_buffer[xm + xs * (yy + ys * zz)];
					realnum&& yn = field_buffer[xx + xs * (yp + ys * zz)] - field_buffer[xx + xs * (ym + ys * zz)];
					realnum&& zn = field_buffer[xx + xs * (yy + ys * zp)] - field_buffer[xx + xs * (yy + ys * zm)];
					realnum gradlen = sqrt(xn * xn + yn * yn + zn * zn + 1e-99);


					realnum eikon;

					realnum phi = phi_buffer[xx + xs * (yy + ys * zz)];
					if (phi < 1e-33) phi = 1e-33;
					if (mUsesCorrection) {
						realnum sumcur = BUF_IDX3D(curvature_buffer, xs, ys, zs, xx, yy, zz);
						eikon = calculate_speed(sumcur, phi, S);
					}
					else {
						eikon = 1 / phi;
					}
					eikon *= gradlen; // normalize
					if (eikon < 1e-11) eikon = 1e-11;

#ifdef USE_VECTOR_AS_SET
					changed_velo[i] = make_pair(pn, eikon);
#else
OMP_CRITICAL_NO_TAG
					changed_velo[pn] = eikon;
#endif
					if (eikon > temp_maxv) temp_maxv = eikon;
				}
			} // for zz
		}
OMP_CRITICAL(maxS_selection)
		{
			if (temp_maxv > maxv) maxv = temp_maxv;
		}
	}
#ifdef DEBUG_CURVATURE
	{
		if (current_iteration % 50 == 0) {
			save_image("Y:/BIOMAG/shortest path/curv_debug/dbg_" + std::to_string(i) + "_" + std::to_string(current_iteration) + "_curv.tif", mCurvatureMap[i]);
			save_image("Y:/BIOMAG/shortest path/curv_debug/dbg_" + std::to_string(i) + "_" + std::to_string(current_iteration) + "_dist.tif", m_distance[i]);
			save_image("Y:/BIOMAG/shortest path/curv_debug/dbg_" + std::to_string(i) + "_" + std::to_string(current_iteration) + "_field.tif", field);
			//save_image("Y:/BIOMAG/shortest path/curv_debug/dbg_" + std::to_string(i) + "_" + std::to_string(current_iteration) + "_sdist.tif", temp_sdist);
		}
	}
#endif
	return maxv;
}
void AreaEikonal::UpdateField(int idx, double maxv) {
	static int n_threads = omp_get_max_threads();
	_PROFILING;
	double* field_buffer = mPhaseFieldMap[idx].GetBufferAsDouble();
	vector<unsigned> size = mPhaseFieldMap[idx].GetSize();
	int xs = size[0], ys = size[1], zs = size[2]; 
	auto& changed_velo = mVelocities[idx];
	auto& act_set = mActivePoints[idx];
#ifdef DEBUG_CURVATURE
	uint8_t* new_update_buffer = m_new_update[idx].GetBufferAsUInt8();
	memset(new_update_buffer, 0, xs * ys * zs);
	double* distance_buffer = m_distance[idx].GetBufferAsDouble();
#endif
OMP_PARALLEL_NUM_THREADS(n_threads/2)
	{
OMP_FOR
#ifdef USE_VECTOR_AS_SET
		for (int i = 0; i < changed_velo.size(); i++) {

			POINT3D pn = changed_velo[i].first;
			double velo = changed_velo[i].second;
#else
		for (int b = 0; b < changed_velo.bucket_count(); b++)
			for (auto bi = changed_velo.begin(b); bi != changed_velo.end(b); bi++) {
				POINT3D pn = bi->first;
				double velo = bi->second;
#endif

				auto [xx, yy, zz] = representation_to_point<int>(pn);
#ifdef DEBUG_CURVATURE
				if (BUF_IDX(distance_buffer, xs, ys, zs, xx, yy, zz) < -0.5)
					BUF_IDX(new_update_buffer, xs, ys, zs, xx, yy, zz) = 1;
#endif
				BUF_IDX3D(field_buffer, xs, ys, zs, xx, yy, zz) += velo;
				if (BUF_IDX3D(field_buffer, xs, ys, zs, xx, yy, zz) >= MAX_FIELD_VAL) {

#ifdef USE_VECTOR_AS_SET
					act_set[i] = REMOVABLE_POINT;
#else
OMP_CRITICAL_NO_TAG
					act_set.erase(pn);
#endif
					if (BUF_IDX3D(field_buffer, xs, ys, zs, xx, yy, zz) >= 1.25)
						BUF_IDX3D(field_buffer, xs, ys, zs, xx, yy, zz) = 1.25;
				}
		}
	}
#ifdef DEBUG_CURVATURE
OMP_CRITICAL_NO_TAG
	{
		if (current_iteration % 50 == 0)
			save_image("Y:/BIOMAG/shortest path/curv_debug/dbg_" + std::to_string(idx) + "_" + std::to_string(current_iteration) + "_updated_pos.tif", m_new_update[idx]);
		if (idx)
			current_iteration++;
	}
#endif
}
void AreaEikonal::UpdateDistance(int idx, realnum current_distance) {

	static int n_threads = omp_get_max_threads();
	_PROFILING;
	Vec3<double> origin_point = !idx ? mStartPoint : mEndPoint;
	auto& changed_velo = mVelocities[idx];
	auto& act_set = mActivePoints[idx];
	double* field_buffer = mPhaseFieldMap[idx].GetBufferAsDouble();
	// mAreaEikonal.m_distance[i]
	double* distance_buffer = mDistanceMap[idx].GetBufferAsDouble();

	// mAreaEikonal.m_distance[(i+1)&1]
	double* counterd_buffer = mDistanceMap[(idx + 1) & 1].GetBufferAsDouble();
	//unsigned int* connectivity_buffer = m_distance_connectivity.GetBufferAsUInt32();
	vector<unsigned> size = mPhaseFieldMap[idx].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
OMP_PARALLEL_NUM_THREADS(n_threads/2)
	{
#ifdef USE_VECTOR_AS_SET
		vector<POINT3D> temp;
		vector<POINT3D> temp_inact;
#else
		unordered_set<POINT3D_SET> temp;
#endif
		temp.reserve(1000);

OMP_FOR_NOWAIT
#ifdef USE_VECTOR_AS_SET
		for (int i = 0; i < changed_velo.size(); i++) {
			POINT3D pn = changed_velo[i].first;
#else
		for (int b = 0; b < changed_velo.bucket_count(); b++)
			for (auto bi = changed_velo.begin(b); bi != changed_velo.end(b); bi++) {
				POINT3D pn = bi->first;

#endif
				auto [xx, yy, zz] = representation_to_point<int>(pn);
				if (field_buffer[xx + xs * (yy + ys * zz)] > 0 && distance_buffer[xx + xs * (yy + ys * zz)] < -0.5f) {
					distance_buffer[xx + xs * (yy + ys * zz)] = current_distance;
					for (int k = -1; k <= 1; k++) {
						if (zz + k >= zs || zz + k < 0) continue;
						for (int j = -1; j <= 1; j++) {
							if (yy + j >= ys || yy + j < 0) continue;
							for (int i = -1; i <= 1; i++) {
								if (xx + i >= xs || xx + i < 0) continue;
								if ((i + 1) % 2 + (j + 1) % 2 + (k + 1) % 2 == 2) {
									if (BUF_IDX3D(field_buffer, xs, ys, zs, xx + i, yy + j, zz + k) < MAX_FIELD_VAL
										&& (mMeetingPlaneCenter.x() < -0.5 || -sgn((mMeetingPlaneNormal * origin_point).Sum() + mMeetingPlaneOffset) * ((mMeetingPlaneNormal * Vec3<double>({ (double)xx + i, (double)yy + j, (double)zz + k })).Sum() + mMeetingPlaneOffset) < MAX_DIST_FROM_PLANE))
									{
										if (xx + i > 0 && xx + i < xs - 1 && yy + j > 0 && yy + j < ys - 1 && zz + k > 0 && zz + k < zs - 1) {

#ifdef USE_VECTOR_AS_SET
											temp.push_back(point_to_representation(xx + i, yy + j, zz + k));
#else
											temp.insert(point_to_representation(xx + i, yy + j, zz + k));
#endif
										}
										else {
											temp_inact.push_back(point_to_representation(xx + i, yy + j, zz + k));
										}
									} // end if
								} // end if
							} // end for i
						} // end for j
					} // end for k
				}
		} // end for i
OMP_CRITICAL_NO_TAG
		{
#ifdef USE_VECTOR_AS_SET
			act_set.insert(act_set.end(), temp.begin(), temp.end());
#else
			act_set.insert(temp.begin(), temp.end());
#endif
		}
	}


	changed_velo.clear();
#ifdef USE_VECTOR_AS_SET
	sort(act_set.begin(), act_set.end());
	act_set.erase(unique(act_set.begin(), act_set.end()), act_set.end());
	if (!act_set.empty() && act_set.back() == REMOVABLE_POINT)
		act_set.pop_back();
#endif
}

bool AreaEikonal::IsDone() const {
	_PROFILING;
	if (mActivePoints[0].empty() && mActivePoints[1].empty())
		return true;
	vector<unsigned> size = mPhaseFieldMap[0].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	const double* dist0_buffer = mDistanceMap[0].GetBufferAsDouble();
	const double* dist1_buffer = mDistanceMap[1].GetBufferAsDouble();
OMP_PARALLEL_FOR
	for (int zyx = 0; zyx < zs*ys*xs; zyx++) {
		if (dist0_buffer[zyx] < 0 && dist1_buffer[zyx] < 0)
			return false;
	}
	return true;
}
void AreaEikonal::Iterate()
{
	_PROFILING;
	// mAreaEikonal.m_field[i]
#ifndef MEETING_PLANE_FROM_INIT_POINTS
	if (++mIterationCount % REFRESH_ITERS == 0 && !mPlaneFinalized) {
		UpdateMeetingPlane();
	}
#endif
	//Calculate the value of capital S
	realnum min_capital_s = 1e-11;
	vector<double> min_capital_s_vals(2);
	_BLOCK_PROFILING;
	OMP_PARALLEL_FOR_NUM_THREADS(2)
		for (int i = 0; i < 2; i++) {
			UpdateCurvature(i);
			min_capital_s_vals[i] = GetMinS(i);
		}
	_UNBLOCK_PROFILING;
	min_capital_s = min_capital_s_vals[0] < min_capital_s_vals[1] ? min_capital_s_vals[0] : min_capital_s_vals[1];
	if (min_capital_s < 1e-11) min_capital_s = 1e-11;
	
	//Calculate the velocity along the front
	vector<double> max_velocities(2);
	_BLOCK_PROFILING;
OMP_PARALLEL_FOR_NUM_THREADS(2)
	for (int i = 0; i < 2; i++) {
		max_velocities[i] = UpdateVelo(i, min_capital_s);
	}
	_UNBLOCK_PROFILING;
	double max_velocity = max_velocities[0] > max_velocities[1] ? max_velocities[0] : max_velocities[1];
	realnum mdatspedmax = 2; 
	max_velocity = mdatspedmax/max_velocity;
	mCurrentDistance += max_velocity;

	// update phase field and distance values
	_BLOCK_PROFILING;
OMP_PARALLEL_FOR_NUM_THREADS(2)
	for (int i = 0; i < 2; i++) {
		UpdateField(i, max_velocity);
		UpdateDistance(i, mCurrentDistance);
	}
	_UNBLOCK_PROFILING;
}
void AreaEikonal::Calculate(const sitk::Image& image, Vec3<double>& point1, Vec3<double>& point2, double beta, double alpha) {
	Initialize(image, point1, point2, beta, alpha);
	Calculate();
}
void AreaEikonal::Calculate() {
	while (!IsDone()) {
		Iterate();
	}
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
Vec3<double> AreaEikonal::GetMeetingPlaneCenter() const {
	return mMeetingPlaneCenter;
}

Vec3<double> AreaEikonal::GetMeetingPlaneNormal() const {
	return mMeetingPlaneNormal;
}

double AreaEikonal::GetMeetingPlaneOffset() const {
	return mMeetingPlaneOffset;
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