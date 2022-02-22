#include "stdafx.h"
#include "AreaEikonal.h"
#include "Utils.h"
#include "math.h"
#include <omp.h>
#include "SimpleITK.h"
#include <vector>
#include <sitkImageOperators.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkFlatStructuringElement.h>
#include "FMSigned.h"
#define BIG_NUMBER 1e11
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

inline unsigned int connectivity_bit(bool val, int x, int y, int z) {

	return (unsigned int)val << (x + 1) + 3 * ((y + 1) + 3 * (z + 1));
}

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

//const int max_threads = omp_get_max_threads();
int a = 0;
CCurvEikonal::CCurvEikonal(void)
{
	

}

CCurvEikonal::~CCurvEikonal(void)
{
}

void CPhaseContainer::CombineDistance() {
	vector<unsigned> size = meeting_plane_positions.GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	double* distance0_buffer = m_distance[0].GetBufferAsDouble();
	double* distance1_buffer = m_distance[1].GetBufferAsDouble();
	double* combined_buffer = m_combined_distance.GetBufferAsDouble();
#pragma omp parallel for FIND_MEET_POINTS_SCHEDULE
	for (int zyx = 0; zyx < zs * ys * xs; zyx++)
	{
		int zz = zyx / (ys * xs);
		{
			int yy = (zyx / xs) % ys;
			{
				int xx = zyx % xs;
				double d0 = distance0_buffer[xx + xs * (yy + ys * zz)];
				double d1 = distance1_buffer[xx + xs * (yy + ys * zz)];
				double newval;
				if (d1 < 0 || d0 < d1 && d0 >= -1e-5)
					newval = d0;
				else
					newval = d1;
				combined_buffer[xx + xs * (yy + ys * zz)] = newval;
			}
		}
	}
}

void CPhaseContainer::ExtractMeetingPlane() {
	_PROFILING;
	vector<unsigned> size= meeting_plane_positions.GetSize();
	double* distance0_buffer = m_distance[0].GetBufferAsDouble();
	double* distance1_buffer = m_distance[1].GetBufferAsDouble();
	int* meeting_plane_buffer = meeting_plane_positions.GetBufferAsInt32();
	int xs(size[0]), ys(size[1]), zs(size[2]);
#pragma omp parallel
	{
		vector<POINT3D> temp;
#pragma omp for
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
						distance0_buffer[xx + xs * (yy + ys * zz)] = distance0_buffer[sample_x + xs * (sample_y + ys * sample_z)];
						distance1_buffer[xx + xs * (yy + ys * zz)] = distance1_buffer[sample_x + xs * (sample_y + ys * sample_z)];
					}
					if (meeting_plane_buffer[xx + xs * (yy + ys * zz)] > 0)
						temp.push_back(point_to_representation(xx, yy, zz));


				}
			}
		}
#pragma omp critical (fill_meeting_plane)
		meeting_plane.insert(temp.begin(), temp.end());
	}
	vector<vector<double>> meeting_plane_vec;
	meeting_plane_vec.reserve(meeting_plane.size());

	std::transform(meeting_plane.begin(), meeting_plane.end(), back_inserter(meeting_plane_vec),
#ifdef STORE_POINT_AS_INTEGER
		[](unsigned long long val) {
			
			auto [x, y, z] = representation_to_point<double>(val);
			vector<double> v({x, y, z});
			return v;
		}
#else
		[](const IPoi3<int>& val) {

			return vector<double>({ (double)val.x, (double)val.y, (double)val.z });
		}
#endif
	);
	std::pair<IPoi3<double>, IPoi3<double>> plane_info = best_plane_from_points(meeting_plane_vec);
	plane_center = plane_info.first;
	plane_normal = plane_info.second;
	double max_norm_val = abs(plane_normal.x);
	int max_norm_sign = sgn(plane_normal.x);
	if (abs(plane_normal.y) > max_norm_val) {
		max_norm_val = abs(plane_normal.y);
		max_norm_sign = sgn(plane_normal.y);
	}
	if (abs(plane_normal.z) > max_norm_val) {
		max_norm_val = abs(plane_normal.z);
		max_norm_sign = sgn(plane_normal.z);
	}
	plane_normal.x *= max_norm_sign;
	plane_normal.y *= max_norm_sign;
	plane_normal.z *= max_norm_sign;
	plane_offset = -(plane_center.x * plane_normal.x + plane_center.y * plane_normal.y + plane_center.z * plane_normal.z);

	vector<double> sample_plane_normal = { 1, 0, 0 };
	// calculate rotation between the two normals
	vector<double> plane_normal = vector<double>({ this->plane_normal.x, this->plane_normal.y, this->plane_normal.z });
	rotation_matrix = rotation_matrix_from_vectors<double>(sample_plane_normal, plane_normal);
}

// Image prep&check

#define GAUTOO 1

// Phasefield stuff

realnum g_w = 2.75;


void CPhaseContainer::FindMeetPoints() {
	_PROFILING;
	vector<unsigned> size = m_field[0].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	int* meeting_plane_buffer = meeting_plane_positions.GetBufferAsInt32();
	for (int ii = 0; ii < 2; ii++) {
		sitk::Image& distance = m_distance[ii];
		sitk::Image& counterd = m_distance[(ii + 1) % 2];
		double* distance_buffer = distance.GetBufferAsDouble();
		double* counterd_buffer = counterd.GetBufferAsDouble();
#pragma omp parallel for FIND_MEET_POINTS_SCHEDULE
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


						if (counterd_buffer[xx + xs * (yy + ys * zz)] >= 0 && !BUF_IDX(meeting_plane_buffer, xs, ys, zs, xx, yy, zz)) {
							if (distance_buffer[xx + xs * (yy + ys * zz)] - counterd_buffer[xx + xs * (yy + ys * zz)] >= 0 &&
								((distance_buffer[xx + xs * (yp + ys * zz)] > 0 && distance_buffer[xx + xs * (yp + ys * zz)] - counterd_buffer[xx + xs * (yp + ys * zz)] < 0) ||
									(distance_buffer[xx + xs * (ym + ys * zz)] > 0 && distance_buffer[xx + xs * (ym + ys * zz)] - counterd_buffer[xx + xs * (ym + ys * zz)] < 0) ||
									(distance_buffer[xp + xs * (yy + ys * zz)] > 0 && distance_buffer[xp + xs * (yy + ys * zz)] - counterd_buffer[xp + xs * (yy + ys * zz)] < 0) ||
									(distance_buffer[xm + xs * (yy + ys * zz)] > 0 && distance_buffer[xm + xs * (yy + ys * zz)] - counterd_buffer[xm + xs * (yy + ys * zz)] < 0) ||
									(distance_buffer[xx + xs * (yy + ys * zp)] > 0 && distance_buffer[xx + xs * (yy + ys * zp)] - counterd_buffer[xx + xs * (yy + ys * zp)] < 0) ||
									(distance_buffer[xx + xs * (yy + ys * zm)] > 0 && distance_buffer[xx + xs * (yy + ys * zm)] - counterd_buffer[xx + xs * (yy + ys * zm)] < 0))) {
								meeting_plane_buffer[xx + xs * (yy + ys * zz)] = 1;
#pragma omp atomic
								n_meet_points++;
							}
						}
					}
				}
			}
		}
	}
}


//-------------------------------------------------------------------------------------------------


void CPhaseContainer::Initialize(sitk::Image data, CVec3& start_point, CVec3& end_point) {
	_PROFILING;
	vector<uint32_t> size = data.GetSize();
	unsigned int spacex(size[0]), spacey(size[1]), spacez(size[2]);
	//g_cyc = 0;
	m_thickstate = sitk::Image(size, sitk::sitkInt32) - 1;
	m_Sumcurvature = sitk::Image(size, sitk::sitkFloat64);
	meeting_plane_positions = sitk::Image(size, sitk::sitkInt32);
	// expansion

	m_temp_sdist[0] = sitk::Image(size, sitk::sitkFloat64);
	m_temp_sdist[1] = sitk::Image(size, sitk::sitkFloat64);
	m_smoothstate = sitk::Image(size, sitk::sitkInt32) - 1;
	// expansion


	m_field[0] = sitk::Image(size, sitk::sitkFloat64) - 1;
	m_field[1] = sitk::Image(size, sitk::sitkFloat64) - 1;
	m_distance[0] = sitk::Image(size, sitk::sitkFloat64) - 1;
	m_distance[1] = sitk::Image(size, sitk::sitkFloat64) - 1;
	m_combined_distance = sitk::Image(size, sitk::sitkFloat64) - 1;
	m_flow_idx = sitk::Image(size, sitk::sitkInt32) - 1;


	m_velo[0] = sitk::Image(size, sitk::sitkFloat64);
	m_velo[1] = sitk::Image(size, sitk::sitkFloat64);
	m_aux = sitk::Image(size, sitk::sitkFloat64);
	m_smoothdist = sitk::Image(size, sitk::sitkFloat64) - 1;
	m_smoothaux = sitk::Image(size, sitk::sitkFloat64);
	m_smoothaux2 = sitk::Image(size, sitk::sitkFloat64);

	m_data = data;
	m_sample_image = sitk::Image(size, sitk::sitkFloat64);
	m_start_point = start_point;
	m_end_point = end_point;

	m_curvature[0] = sitk::Image(size, sitk::sitkFloat64);
	m_curvature[1] = sitk::Image(size, sitk::sitkFloat64);
#ifdef DEBUG_CURVATURE
	
	m_new_update[0] = sitk::Image(size, sitk::sitkUInt8);
	m_new_update[1] = sitk::Image(size, sitk::sitkUInt8);
#endif
	double xdist = start_point.x - end_point.x;
	double ydist = start_point.y - end_point.y;
	double zdist = start_point.z - end_point.z;
	double point_dist = sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
	double init_dist = sqrt(point_dist);
	int hs = INIT_DIST+2;//5 (int)(1.5f*g_w/2);
	m_currentdistance = 0;
	for (int ii = 0; ii < 2; ++ii) {
		CVec3 point = !ii ? start_point : end_point;
		double* field_buffer = m_field[ii].GetBufferAsDouble();
		double* distance_buffer = m_distance[ii].GetBufferAsDouble();
		for (int zz = point.z - hs; zz < point.z + hs; ++zz) {
			if (zz < 0 || zz >= spacez) continue;
			for (int yy = point.y - hs; yy < point.y + hs; ++yy) {
				if (yy < 0 || yy >=spacey) continue;
				for (int xx = point.x - hs; xx < point.x + hs; ++xx) {
					if (xx < 0 || xx >= spacex) continue;
					int dx = xx - point.x, dy = yy - point.y;
					int dz = zz - point.z;
					int dd = dx * dx + dy * dy + dz * dz;
					if ((int)dd <= INIT_DIST * INIT_DIST) {
						double dd_sqrt = sqrt(dd);
						if (dd > (INIT_DIST-1)*(INIT_DIST-1)) {
							field_buffer[xx + spacex * (yy + spacey * zz)] = INIT_DIST-dd_sqrt;
#ifdef USE_VECTOR_AS_SET
							active_set[ii].push_back(point_to_representation(xx, yy, zz));
#else
							active_set[ii].insert(point_to_representation(xx, yy, zz));
#endif
						}
						else {
							field_buffer[xx + spacex * (yy + spacey * zz)] = 1.0;
						}
						distance_buffer[xx + spacex * (yy + spacey * zz)] = dd_sqrt;
						if (m_currentdistance < dd_sqrt) m_currentdistance = dd_sqrt;

						// curvatures
						//if (dd < 1e-11) dd = 1e-11; dd = 1.0/dd;
						//m_phasefield.m_thickstate[ii][zz][yy][xx] = 2;
						//m_phasefield.m_Sumcurvature[ii][zz][yy][xx] = dd+dd;
					}
					else if (dd < (INIT_DIST+1)*(INIT_DIST+1)) {
#ifdef USE_VECTOR_AS_SET
						active_set[ii].push_back(point_to_representation(xx, yy, zz));
#else
						active_set[ii].insert(point_to_representation(xx, yy, zz));
#endif
						BUF_IDX(field_buffer, spacex, spacey, spacez, xx, yy, zz) = INIT_DIST-sqrt(dd);
					}
				}
			}
		}

	}
	n_min_meet_points = spacex * spacey < spacex* spacez ? spacex * spacey :
		spacex * spacez < spacey* spacez ? spacex * spacez : spacey * spacez;
	n_min_meet_points /= 4;

	m_bdone = false;
}

void CPhaseContainer::Initialize(CPhaseContainer& phasefield, vector<double>& rotation_matrix, bool inverse) {
	_PROFILING;
	vector<unsigned int> sample_size;
	resample_img(phasefield.m_thickstate, m_thickstate, rotation_matrix, sample_size, inverse);
	resample_img(phasefield.m_Sumcurvature, m_Sumcurvature, rotation_matrix, sample_size, inverse);
	neg_to_minus1(m_Sumcurvature);

	resample_img<sitk::sitkNearestNeighbor>(phasefield.meeting_plane_positions, meeting_plane_positions, rotation_matrix, sample_size, inverse);
	// expansion

	m_temp_sdist[0] = sitk::Image(m_thickstate.GetSize(), sitk::sitkFloat64);
	m_temp_sdist[1] = sitk::Image(m_thickstate.GetSize(), sitk::sitkFloat64);
	resample_img<sitk::sitkNearestNeighbor>(phasefield.m_smoothstate, m_smoothstate, rotation_matrix, sample_size, inverse);
	// expansion


	resample_img(phasefield.m_field[0], m_field[0], rotation_matrix, sample_size, inverse);
	neg_to_minus1(m_field[0]);
	resample_img(phasefield.m_field[1], m_field[1], rotation_matrix, sample_size, inverse);
	neg_to_minus1(m_field[1]);
	resample_img(phasefield.m_distance[0], m_distance[0], rotation_matrix, sample_size, inverse, true);
	neg_to_minus1(m_distance[0]);
	resample_img(phasefield.m_distance[1], m_distance[1], rotation_matrix, sample_size, inverse, true);
	neg_to_minus1(m_distance[1]);
	m_sample_image = resample_img(phasefield.m_combined_distance, m_combined_distance, rotation_matrix, sample_size, inverse, true);
	neg_to_minus1(m_combined_distance);
	resample_img<sitk::sitkNearestNeighbor>(phasefield.m_flow_idx, m_flow_idx, rotation_matrix, sample_size, inverse, true);

	resample_img(phasefield.m_velo[0], m_velo[0], rotation_matrix, sample_size, inverse);
	resample_img(phasefield.m_velo[1], m_velo[1], rotation_matrix, sample_size, inverse);
	resample_img(phasefield.m_aux, m_aux, rotation_matrix, sample_size, inverse);
	resample_img(phasefield.m_smoothdist, m_smoothdist, rotation_matrix, sample_size, inverse);
	neg_to_minus1(m_smoothdist);
	resample_img(phasefield.m_smoothaux, m_smoothaux, rotation_matrix, sample_size, inverse);
	resample_img(phasefield.m_smoothaux2, m_smoothaux2, rotation_matrix, sample_size, inverse);

	resample_img(phasefield.m_data, m_data, rotation_matrix, sample_size, inverse, true);

	m_currentdistance = phasefield.m_currentdistance;

	m_curvature[0] = sitk::Image(m_thickstate.GetSize(), sitk::sitkFloat64);
	m_curvature[1] = sitk::Image(m_thickstate.GetSize(), sitk::sitkFloat64);
#ifdef DEBUG_CURVATURE
	
	m_new_update[0] = sitk::Image(m_thickstate.GetSize(), sitk::sitkUInt8);
	m_new_update[1] = sitk::Image(m_thickstate.GetSize(), sitk::sitkUInt8);
#endif
	m_bdone = phasefield.m_bdone;
	vector<unsigned> new_size = m_sample_image.GetSize();
	n_min_meet_points = new_size[0] * new_size[1] < new_size[0] * new_size[2] ? new_size[0] * new_size[1] :
		new_size[0] * new_size[2] < new_size[1] * new_size[2] ? new_size[0] * new_size[2] : new_size[1] * new_size[2];
	n_min_meet_points /= 4;

	for (int ii = 0; ii < 2; ii++) {
		for (auto it = phasefield.active_set[ii].begin(); it != phasefield.active_set[ii].end(); it++) {
			auto [xx, yy, zz] = representation_to_point<double>(*it);
			vector<double> pos({ xx, yy, zz });
			pos = phasefield.m_data.TransformContinuousIndexToPhysicalPoint(pos);
			pos = m_sample_image.TransformPhysicalPointToContinuousIndex(pos);
#ifdef USE_VECTOR_AS_SET
			active_set[ii].push_back(point_to_representation(pos[0], pos[1], pos[2]));
#else
			active_set[ii].insert(point_to_representation(pos[0], pos[1], pos[2]));
#endif
		}
	}
}


bool g_modeswitch = false;
int current_iteration = 0;

void CPhaseContainer::UpdateCurvature(int i) {
	static int n_threads = omp_get_max_threads();
	sitk::Image& field = m_field[i];
	auto& act_set = active_set[i];
	const vector<unsigned>& size = m_field[i].GetSize();
	int xs = size[0], ys = size[1], zs = size[2]; 
	sitk::Image& temp_sdist = FMSigned::build_for_neighbors(i, field, 2, 1, act_set);
	double* sdist_buffer = temp_sdist.GetBufferAsDouble();
	double* curvature_buffer = m_curvature[i].GetBufferAsDouble();
	memset(curvature_buffer, 0, xs * ys * zs * sizeof(double));
#pragma omp parallel for num_threads(n_threads/2)
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
				if (zz > 0 && zz < zs - 1)
				{
					if (yy > 0 && yy < ys - 1)
					{
						if (xx > 0 && xx < xs - 1)
							//for (int zz = 1; zz < zs - 1; ++zz) {
							//	for (int yy = 1; yy < ys - 1; ++yy) {
							//		for (int xx = 1; xx < xs - 1; ++xx) {
						{
							int zp = zz + 1, zm = zz - 1;
							int yp = yy + 1, ym = yy - 1;
							int xp = xx + 1, xm = xx - 1;


							realnum sumcur =
								BUF_IDX(sdist_buffer, xs, ys, zs, xp, yy, zz) +
								BUF_IDX(sdist_buffer, xs, ys, zs, xm, yy, zz) +
								BUF_IDX(sdist_buffer, xs, ys, zs, xx, yp, zz) +
								BUF_IDX(sdist_buffer, xs, ys, zs, xx, ym, zz) +
								BUF_IDX(sdist_buffer, xs, ys, zs, xx, yy, zp) +
								BUF_IDX(sdist_buffer, xs, ys, zs, xx, yy, zm) -
								6 * BUF_IDX(sdist_buffer, xs, ys, zs, xx, yy, zz);
							//if (sumcur < -2) sumcur == -2;
							BUF_IDX(curvature_buffer, xs, ys, zs, xx, yy, zz) = sumcur;
						}
					} 
				} 
			}
		} 
	}
}

realnum CPhaseContainer::GetMinS(int i) {

	static int n_threads = omp_get_max_threads();
	sitk::Image& field = m_field[i];
	auto& act_set = active_set[i];
	//int n_top = 100;
	int n_top = act_set.size();
	vector<double> topSs;
	topSs.reserve(n_top + 1);
	double* field_buffer = field.GetBufferAsDouble();
	double* data_buffer = m_data.GetBufferAsDouble();
	double* curvature_buffer = m_curvature[i].GetBufferAsDouble();
	const vector<unsigned>& size = m_field[i].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	int n_outliers(0);
	realnum minS(DBL_MAX);
	if (act_set.size() == 0)
		return minS;
	vector<vector<double>> heaps(n_threads / 2);

#pragma omp parallel num_threads(n_threads/2)
	{
		//double temp_minS(DBL_MAX);
		vector<double>& temp = heaps[omp_get_thread_num()];
		temp.reserve(act_set.size()/(n_threads/2));
#pragma omp for reduction(+: n_outliers)
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
					if (zz > 0 && zz < zs - 1)
					{
						if (yy > 0 && yy < ys - 1)
						{
							if (xx > 0 && xx < xs - 1)
								//for (int zz = 1; zz < zs - 1; ++zz) {
								//	for (int yy = 1; yy < ys - 1; ++yy) {
								//		for (int xx = 1; xx < xs - 1; ++xx) {
							{
								int zp = zz + 1, zm = zz - 1;
								int yp = yy + 1, ym = yy - 1;
								int xp = xx + 1, xm = xx - 1;

								realnum d0 = data_buffer[xx + xs * (yy + ys * zz)];
								if (d0 < 1e-33) d0 = 1e-33;
								double sumcur = BUF_IDX(curvature_buffer, xs, ys, zs, xx, yy, zz);
								realnum&& xn = field_buffer[xp + xs * (yy + ys * zz)] - field_buffer[xm + xs * (yy + ys * zz)];
								realnum&& yn = field_buffer[xx + xs * (yp + ys * zz)] - field_buffer[xx + xs * (ym + ys * zz)];
								realnum&& zn = field_buffer[xx + xs * (yy + ys * zp)] - field_buffer[xx + xs * (yy + ys * zm)];
								realnum gradlen = sqrt(xn * xn + yn * yn + zn * zn + 1e-99);

								if (is_outlier(sumcur, 2/gradlen))
								{
									n_outliers++;
								}
								realnum S = calculate_S(sumcur, 1 / d0, 2/gradlen);
								temp.push_back(S);
								//loop body
							}
						} //for xx
					} // for yy
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
	/*minS = 0;
	std::for_each(topSs.begin(), topSs.end(), [&minS](double val) {minS += val; });
	minS /= topSs.size();*/
	/*if (n_outliers == 0)
		minS = topSs[0];
	else
		minS = topSs[n_outliers*2];*/
	size_t idx = min((size_t)100, topSs.size()-1);
	minS = topSs[idx];
	return minS;
}
realnum CPhaseContainer::UpdateVelo(int i, bool use_correction, double S) {
	static int n_threads = omp_get_max_threads();
	//_PROFILING;
	sitk::Image& field = m_field[i];
	double* field_buffer = field.GetBufferAsDouble();

	double* data_buffer = m_data.GetBufferAsDouble();
	

	const vector<unsigned>& size = m_field[i].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];

	////////////////////////////////////////////////////////////////////
	// evolution logic
	////////////////////////////////////////////////////////////////////
	auto& act_set = active_set[i];
	
	//realnum maxv(0);
	double* curvature_buffer = m_curvature[i].GetBufferAsDouble();

	// Inhomogeneous
	// Calculate the velocity at every point (velo[zz][yy][xx]), and the maximum velocity (maxv)
	auto& changed_velo = m_changed_velo[i];
#ifdef USE_VECTOR_AS_SET
	changed_velo.resize(act_set.size());
#else
	changed_velo.reserve(act_set.size());
#endif
	int c = 0;

	double maxv(0);

#pragma omp parallel num_threads(n_threads/2)
	{
		double temp_maxv(0);
#pragma omp for
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
					if (zz > 0 && zz < zs - 1)
					{
						if (yy > 0 && yy < ys - 1)
						{
							if (xx > 0 && xx < xs - 1)
								//for (int zz = 1; zz < zs - 1; ++zz) {
								//	for (int yy = 1; yy < ys - 1; ++yy) {
								//		for (int xx = 1; xx < xs - 1; ++xx) {
							{
								int zp = zz + 1, zm = zz - 1;
								int yp = yy + 1, ym = yy - 1;
								int xp = xx + 1, xm = xx - 1;

								realnum&& xn = field_buffer[xp + xs * (yy + ys * zz)] - field_buffer[xm + xs * (yy + ys * zz)];
								realnum&& yn = field_buffer[xx + xs * (yp + ys * zz)] - field_buffer[xx + xs * (ym + ys * zz)];
								realnum&& zn = field_buffer[xx + xs * (yy + ys * zp)] - field_buffer[xx + xs * (yy + ys * zm)];
								realnum gradlen = sqrt(xn * xn + yn * yn + zn * zn + 1e-99);

								realnum sumcur = BUF_IDX(curvature_buffer, xs, ys, zs, xx, yy, zz);

								realnum eikon;

								realnum d0 = data_buffer[xx + xs * (yy + ys * zz)];
								if (d0 < 1e-33) d0 = 1e-33;

								eikon = calculate_speed(sumcur, 1 / d0, S);
								eikon *= gradlen; // normalize
								if (eikon < 1e-11) eikon = 1e-11;

#ifdef USE_VECTOR_AS_SET
								changed_velo[i] = make_pair(pn, eikon);
#else
#pragma omp critical
								changed_velo[pn] = eikon;
#endif
								if (eikon > temp_maxv) temp_maxv = eikon;
								//loop body
							}
						} //for xx
					} // for yy
				}
			} // for zz
		}
#pragma omp critical (maxS_selection)
		{
			if (temp_maxv > maxv) maxv = temp_maxv;
		}
	}
#ifdef DEBUG_CURVATURE
	{
		if (current_iteration % 50 == 0) {
			save_image("Y:/BIOMAG/shortest path/curv_debug/dbg_" + std::to_string(i) + "_" + std::to_string(current_iteration) + "_curv.tif", m_curvature[i]);
			save_image("Y:/BIOMAG/shortest path/curv_debug/dbg_" + std::to_string(i) + "_" + std::to_string(current_iteration) + "_dist.tif", m_distance[i]);
			save_image("Y:/BIOMAG/shortest path/curv_debug/dbg_" + std::to_string(i) + "_" + std::to_string(current_iteration) + "_field.tif", field);
			//save_image("Y:/BIOMAG/shortest path/curv_debug/dbg_" + std::to_string(i) + "_" + std::to_string(current_iteration) + "_sdist.tif", temp_sdist);
		}
	}
#endif
	return maxv;
}
void CPhaseContainer::UpdateField(int idx, double maxv) {
	static int n_threads = omp_get_max_threads();
	_PROFILING;
	double* field_buffer = m_field[idx].GetBufferAsDouble();
	vector<unsigned> size = m_field[idx].GetSize();
	int xs = size[0], ys = size[1], zs = size[2]; 
	auto& changed_velo = m_changed_velo[idx];
	auto& act_set = active_set[idx];
#ifdef DEBUG_CURVATURE
	uint8_t* new_update_buffer = m_new_update[idx].GetBufferAsUInt8();
	memset(new_update_buffer, 0, xs * ys * zs);
	double* distance_buffer = m_distance[idx].GetBufferAsDouble();
#endif
#pragma omp parallel num_threads(n_threads/2)
	{
#pragma omp for
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
				BUF_IDX(field_buffer, xs, ys, zs, xx, yy, zz) += velo;
				if (BUF_IDX(field_buffer, xs, ys, zs, xx, yy, zz) >= MAX_FIELD_VAL) {

#ifdef USE_VECTOR_AS_SET
					act_set[i] = REMOVABLE_POINT;
#else
#pragma omp critical
					act_set.erase(pn);
#endif
					if (BUF_IDX(field_buffer, xs, ys, zs, xx, yy, zz) >= 1.25)
						BUF_IDX(field_buffer, xs, ys, zs, xx, yy, zz) = 1.25;
				}
		}
	}
#ifdef DEBUG_CURVATURE
#pragma omp critical
	{
		if (current_iteration % 50 == 0)
			save_image("Y:/BIOMAG/shortest path/curv_debug/dbg_" + std::to_string(idx) + "_" + std::to_string(current_iteration) + "_updated_pos.tif", m_new_update[idx]);
		if (idx)
			current_iteration++;
	}
#endif
}
void CPhaseContainer::UpdateDistance(int idx, realnum current_distance) {

	static int n_threads = omp_get_max_threads();
	_PROFILING;
	CVec3 origin_point = !idx ? m_start_point : m_end_point;
	auto& changed_velo = m_changed_velo[idx];
	auto& act_set = active_set[idx];
	double* field_buffer = m_field[idx].GetBufferAsDouble();
	// m_phasefield.m_distance[i]
	double* distance_buffer = m_distance[idx].GetBufferAsDouble();

	// m_phasefield.m_distance[(i+1)&1]
	double* counterd_buffer = m_distance[(idx + 1) & 1].GetBufferAsDouble();
	//unsigned int* connectivity_buffer = m_distance_connectivity.GetBufferAsUInt32();
	vector<unsigned> size = m_field[idx].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
#pragma omp parallel num_threads(n_threads/2)
	{
#ifdef USE_VECTOR_AS_SET
		vector<POINT3D> temp;
#else
		unordered_set<POINT3D_SET> temp;
#endif
		temp.reserve(1000);

#pragma omp for nowait
#ifdef USE_VECTOR_AS_SET
		for (int i = 0; i < changed_velo.size(); i++) {
			POINT3D pn = changed_velo[i].first;
#else
		for (int b = 0; b < changed_velo.bucket_count(); b++)
			for (auto bi = changed_velo.begin(b); bi != changed_velo.end(b); bi++) {
				POINT3D pn = bi->first;

#endif
				auto [xx, yy, zz] = representation_to_point<int>(pn);
				if (zz == 0 || zz >= zs - 1) continue;
				{
					if (yy == 0 || yy >= ys - 1) continue;
					{
						if (xx == 0 || xx >= xs - 1) continue;
						/*for (int zz = 1; zz < zs - 1; ++zz) {
							for (int yy = 0 + 1; yy < ys - 1; ++yy) {
								for (int xx = 0 + 1; xx < xs - 1; ++xx) {*/
						if (field_buffer[xx + xs * (yy + ys * zz)] > 0 && distance_buffer[xx + xs * (yy + ys * zz)] < -0.5f) {
							distance_buffer[xx + xs * (yy + ys * zz)] = current_distance;
							for (int k = -1; k <= 1; k++) {
								for (int j = -1; j <= 1; j++) {
									for (int i = -1; i <= 1; i++) {
										if ((i + 1) % 2 + (j + 1) % 2 + (k + 1) % 2 == 2) {
											if (BUF_IDX(field_buffer, xs, ys, zs, xx + i, yy + j, zz + k) < MAX_FIELD_VAL
												&& (plane_center.x < -0.5 || sgn(plane_normal.x * origin_point.x + plane_normal.y * origin_point.y + plane_normal.z * origin_point.z + plane_offset) * (plane_normal.x * xx + i + plane_normal.y * yy + j + plane_normal.z * zz + k + plane_offset) + MAX_DIST_FROM_PLANE > 0)
												&& xx+i > 0 && xx+i < xs-1 && yy + j > 0 && yy + j < ys - 1 && zz + k > 0 && zz + k < zs - 1
												) {

#ifdef USE_VECTOR_AS_SET
												temp.push_back(point_to_representation(xx + i, yy + j, zz + k));
#else
												temp.insert(point_to_representation(xx + i, yy + j, zz + k));
#endif
											} // end if
										} // end if
									} // end for i
								} // end for j
							} // end for k
						} // end if xx
					} // end if yy
				} // end if zz
		} // end for i
#pragma omp critical
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

bool CPhaseContainer::IsDone() {
	_PROFILING;
	if (active_set[0].empty() && active_set[1].empty())
		return true;
	vector<unsigned> size = m_field[0].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	double* dist0_buffer = m_distance[0].GetBufferAsDouble();
	double* dist1_buffer = m_distance[1].GetBufferAsDouble();
	for (int zz = 1; zz < zs - 1; zz++) {
		for (int yy = 1; yy < ys - 1; yy++) {
			for (int xx = 1; xx < xs - 1; xx++) {
				if (BUF_IDX(dist0_buffer, xs, ys, zs, xx, yy, zz) < 0 && BUF_IDX(dist1_buffer, xs, ys, zs, xx, yy, zz) < 0)
					return false;
			}
		}
	}
	return true;
}
void CPhaseContainer::Iterate(bool use_correction)
{
	_PROFILING;
	// m_phasefield.m_field[i]
	if(!n_plane_finalized)
		if (++m_counter % REFRESH_ITERS == 0) {
			int n_prev_meet_points = n_meet_points >= 100? n_meet_points : 100;
			FindMeetPoints();
			if(n_meet_points > n_min_meet_points)
				ExtractMeetingPlane();
			if (n_meet_points == n_prev_meet_points && n_meet_points > 0) {
				n_plane_finalized = true;
			}
			
		}
	vector<double> maxvs(2);
	vector<double> minSs(2);
	_BLOCK_PROFILING;
#pragma omp parallel for num_threads(2)
	for (int i = 0; i < 2; i++) {
		UpdateCurvature(i);
		minSs[i] = GetMinS(i);
	}
	_UNBLOCK_PROFILING;
	realnum minS = minSs[0] < minSs[1] ? minSs[0] : minSs[1];
	if (minS < 1e-11) minS = 1e-11;
	_BLOCK_PROFILING;
#pragma omp parallel for num_threads(2)
	for (int i = 0; i < 2; i++) {
		maxvs[i] = UpdateVelo(i, use_correction, minS);
	}
	_UNBLOCK_PROFILING;
	double maxv = maxvs[0] > maxvs[1] ? maxvs[0] : maxvs[1];
	realnum mdatspedmax = 2; //*2 set max speed
	maxv = mdatspedmax/maxv;
	m_currentdistance += maxv;
	//m_currentdistance += 1;

	// update phase field values
	_BLOCK_PROFILING;
#pragma omp parallel for num_threads(2)
	for (int i = 0; i < 2; i++) {
		UpdateField(i, maxv);
		UpdateDistance(i, m_currentdistance);
	}
	_UNBLOCK_PROFILING;
	m_bdone = IsDone();
	// TODO: check for collision
	// updating the distance map to the current value, and checking if it is complete
	
	////////////////////////////////////////////////////////////////////
	// phase field regularization w/o velocity zeroed out
	////////////////////////////////////////////////////////////////////


}
void CPhaseContainer::CalculateAlignedCombinedDistance(double p0_x, double p1_x) {
	_PROFILING;
	vector<unsigned> size = m_data.GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	m_combined_distance = sitk::Image({ (unsigned) xs, (unsigned)ys, (unsigned)zs }, sitk::sitkFloat64) - 1;
	m_flow_idx = sitk::Image({ (unsigned)xs, (unsigned)ys, (unsigned)zs }, sitk::sitkInt32) - 1;
	double* dist0_buffer = m_distance[0].GetBufferAsDouble();
	double* dist1_buffer = m_distance[1].GetBufferAsDouble();
	double* combined_dist_buffer = m_combined_distance.GetBufferAsDouble();
	int* flow_idx_buffer = m_flow_idx.GetBufferAsInt32();

	for (int zz = 0; zz < zs; zz++) {
		for (int yy = 0; yy < ys; yy++) {
			for (int xx = 0; xx < xs; xx++) {
				int xx_sgn = sgn(m_plane_slice - xx);
				int p0_sgn = sgn(m_plane_slice - p0_x);
				int p1_sgn = sgn(m_plane_slice - p1_x);
				if (abs(m_plane_slice - xx) < 0.5) {
					combined_dist_buffer[xx + xs * (yy + ys * zz)] = (dist0_buffer[xx + p0_sgn + xs * (yy + ys * zz)] + dist1_buffer[xx + p1_sgn + xs * (yy + ys * zz)]);
				}
				else if (p0_sgn == xx_sgn) {
					combined_dist_buffer[xx + xs * (yy + ys * zz)] = dist0_buffer[xx + xs * (yy + ys * zz)];
					flow_idx_buffer[xx + xs * (yy + ys * zz)] = 0;
				}
				else if (p1_sgn == xx_sgn) {
					combined_dist_buffer[xx + xs * (yy + ys * zz)] = dist1_buffer[xx + xs * (yy + ys * zz)];
					flow_idx_buffer[xx + xs * (yy + ys * zz)] = 1;
				}
			}
		}
	}
}

// only for xslice
void CPlanePhaseField::GetDistancemean(SVoxImg<SWorkImg<realnum>>& distance1, SVoxImg<SWorkImg<realnum>>& distance2,  int xslice)
{
	_PROFILING;
	//TODO: instead of searching for the selected x slice, locate the meeting surface of the two processes, and use that
	int xs = distance1.xs, ys = distance1.ys, zs = distance1.zs;
	m_nall = m_npos = m_since_last_update = 0;
	m_bInited = false;
	if (!xslice) return;

	m_gx2D.Set(ys, zs, 0.0f);
	m_gy2D.Set(ys, zs, 0.0f);
	m_field2D.Set(ys, zs, 1.0f);
	m_velo2D.Set(ys, zs, 0.0f);

	for (int zz = 1; zz < zs-1; ++zz) {
		for (int yy = 1; yy < ys-1; ++yy) {
			{
				int xx = xslice;
				int yp = yy + 1; if (yp > ys - 2) yp = ys - 2;
				int ym = yy - 1; if (ym < 1) ym = 1;
				int zp = zz + 1; if (zp > zs - 2) zp = zs - 2;
				int zm = zz - 1; if (zm < 1) zm = 1;
				m_gx2D[zz][yy] = (distance1[zz][yp][xx] - distance1[zz][ym][xx] + distance2[zz][yp][xx] - distance2[zz][ym][xx]);
				m_gy2D[zz][yy] = (distance1[zp][yy][xx] - distance1[zm][yy][xx] + distance2[zp][yy][xx] - distance2[zm][yy][xx]);
			}
			//if current point is not on the edge of the data
			if (!(zz < 5 || zz >= zs - 5 || yy < 5 || yy >= ys - 5))
				m_field2D[zz][yy] = -1.0f;
			else ++m_npos;
			++m_nall;
		}
	}

	m_nect = 100;
	m_nfct = 20000;
	m_bInited = true;

}

void CPlanePhaseField::GetDistancemean(sitk::Image& distance, int xslice)
{
	vector<unsigned> size = distance.GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	double* distance_buffer = distance.GetBufferAsDouble();
	m_nall = m_npos = m_since_last_update = 0;
	m_bInited = false;
	if (!xslice) return;

	m_gx2D.Set(ys, zs, 0.0f);
	m_gy2D.Set(ys, zs, 0.0f);
	m_field2D.Set(ys, zs, 1.0f);
	m_velo2D.Set(ys, zs, 0.0f);
	m_distance2D[0].Set(ys, zs, -1);
	m_distance2D[1].Set(ys, zs, -1);


	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			{
				int xx = xslice;
				m_distance2D[0][zz][yy] = distance_buffer[xx - 1 + xs * (yy + ys * zz)];
				m_distance2D[1][zz][yy] = distance_buffer[xx + 1 + xs * (yy + ys * zz)];
				if (zz == 0 || yy == 0 || zz == zs - 1 || yy == ys - 1) continue;
				int yp = yy + 1; if (yp > ys - 2) yp = ys - 2;
				int ym = yy - 1; if (ym < 1) ym = 1;
				int zp = zz + 1; if (zp > zs - 2) zp = zs - 2;
				int zm = zz - 1; if (zm < 1) zm = 1;
				m_gx2D[zz][yy] = distance_buffer[xx + xs * (yp + ys * zz)] - distance_buffer[xx + xs * (ym + ys * zz)];
				m_gy2D[zz][yy] = distance_buffer[xx + xs * (yy + ys * zp)] - distance_buffer[xx + xs * (yy + ys * zm)];
			}
			//if current point is not on the edge of the data
			if (!(zz < 5 || zz >= zs - 5 || yy < 5 || yy >= ys - 5))
				m_field2D[zz][yy] = -1.0f;
			else ++m_npos;
			++m_nall;
		}
	}

	m_nect = 100;
	m_nfct = 20000;
	m_bInited = true;

}

void CPlanePhaseField::Iterate() 
{
	_PROFILING;
	if (!m_bInited) return;

	int xs = m_field2D.xs, ys = m_field2D.ys;
	m_velo2D.GetLaplace(m_field2D);
	m_velo2D *= 4.0f;

	realnum fac = 1.0f;
	realnum maxv(0.0f);
	#pragma omp parallel
	{
		double temp_maxv(0);
#pragma omp for
		for (int yx = 0; yx < xs * ys; yx++) {
			int yy = yx / xs;
			if (yy > 0 && yy < ys - 1) {
				int xx = yx % xs;
				if (xx > 0 && xx < xs - 1) {
					realnum fmm = m_field2D[yy][xx];
					m_velo2D[yy][xx] -= fac * fmm * fmm * fmm;
					m_velo2D[yy][xx] += fac * fmm;
					realnum gfx = 0.5f * (m_field2D[yy][xx + 1] - m_field2D[yy][xx - 1]);
					realnum gfy = 0.5f * (m_field2D[yy + 1][xx] - m_field2D[yy - 1][xx]);

					m_velo2D[yy][xx] += 10.0f * (gfx * m_gx2D[yy][xx] + gfy * m_gy2D[yy][xx]);
					realnum cv = m_velo2D[yy][xx];
					if (cv < 0)
						cv *= -1;
					if (cv > temp_maxv)
						temp_maxv = cv;
				}
			}
		}
#pragma omp critical
		if (temp_maxv > maxv)
			maxv = temp_maxv;
	}

	if (maxv > 1e-5f)
		maxv = 0.025f / maxv;
	int npos(0);

	#pragma omp parallel for reduction(+: npos)
	for (int yy = 1; yy < ys - 1; ++yy) {
		for (int xx = 1; xx < xs - 1; ++xx) {
			m_field2D[yy][xx] += m_velo2D[yy][xx]*maxv;
			if (m_field2D[yy][xx] > 0) ++npos;
		}
	}
	if (npos > m_npos) 
		m_since_last_update = 0;
	else if (++m_since_last_update >= 200)
		m_nect = 0;
	
	m_npos = npos;
	/*--m_nect;
	if (!m_nect) {
		realnum relpos = ((realnum)abs(npos - m_npos)) / ((realnum)m_nall);
		m_nect = 100;
		m_npos = npos;
		if (relpos < 0.0003f) {
			m_nect = 0;
		}
	}*/

	--m_nfct;
	if (m_nfct <= 0)
		m_nect = 0;

}

std::unordered_set<unsigned long>& CPlanePhaseField::RetrieveBound()
{
	_PROFILING;
	int xs = m_field2D.xs, ys = m_field2D.ys;
	m_bound.clear();
	for (int yy = 1; yy < ys - 1; ++yy) {
		for (int xx = 1; xx < xs - 1; ++xx) {

			if (m_field2D[yy][xx] > 0) {
				if (m_field2D[yy + 1][xx] <= 0 || m_field2D[yy - 1][xx] <= 0
					|| m_field2D[yy][xx + 1] <= 0 || m_field2D[yy][xx - 1] <= 0)
					m_bound.emplace((yy << 16) + xx);
			}

		}
	}
	return m_bound;
}


CVec3 operator *(double f, CVec3 &v)
{
	return CVec3(f*v.x,f*v.y,f*v.z);
}
CVec3 operator +(CVec3 &v, CVec3 &w)
{
	return CVec3(v.x+w.x,v.y+w.y,v.z+w.z);
}
