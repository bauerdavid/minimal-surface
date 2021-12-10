#include "stdafx.h"
#include "AreaEikonal.h"
#include "Utils.h"
#include "math.h"
#include <omp.h>
#include "SimpleITK.h"
#include <vector>
#include <sitkImageOperators.h>
#define BIG_NUMBER 1e11
#define MAX_THREADS 32
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

void CCurvEikonal::ExtractMeetingPlane() {
	m_phasefield.FindMeetPoints();
	vector<unsigned> size= m_phasefield.meeting_plane_positions.GetSize();
	double* distance0_buffer = m_phasefield.m_distance[0].GetBufferAsDouble();
	double* distance1_buffer = m_phasefield.m_distance[1].GetBufferAsDouble();
	double* combined_buffer = m_phasefield.m_combined_distance.GetBufferAsDouble();
	int* meeting_plane_buffer = m_phasefield.meeting_plane_positions.GetBufferAsInt32();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
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
					meeting_plane.insert(IPoi3<double>(xx, yy, zz));
				double newval;
				double d0 = distance0_buffer[xx + xs * (yy + ys * zz)];
				double d1 = distance1_buffer[xx + xs * (yy + ys * zz)];
				if (d1 < 0 || d0 < d1 && d0 >= -1e-5)
					newval = d0;
				else
					newval = d1;
				combined_buffer[xx + xs * (yy + ys * zz)] = newval;
			}
		}
	}

	std::pair<IPoi3<double>, IPoi3<double>> plane_info = best_plane_from_points(meeting_plane);
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
	vector<unsigned> size = m_field[0].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	int* meeting_plane_buffer = meeting_plane_positions.GetBufferAsInt32();
	for (int ii = 0; ii < 2; ii++) {
		/*for (int zyx = 0; zyx < zs * ys * xs; zyx++) {
			int zz = zyx / (ys * xs);
			if (zz == 0 || zz >= zs - 1) continue;
			{
				int yy = (zyx / xs) % ys;
				if (yy == 0 || yy >= ys - 1) continue;
				{
					int xx = zyx % xs;
					if (xx == 0 || xx >= xs - 1) continue;*/
		sitk::Image& distance = m_distance[ii];
		sitk::Image& counterd = m_distance[(ii + 1) % 2];
		double* distance_buffer = distance.GetBufferAsDouble();
		double* counterd_buffer = counterd.GetBufferAsDouble();
#pragma omp parallel for schedule(static)

		for (int zz = 1; zz < zs - 1; ++zz) {
			for (int yy = 1; yy < ys - 1; ++yy) {
				for (int xx = 1; xx < xs - 1; ++xx) {
					int zp = zz + 1, zm = zz - 1;
					int yp = yy + 1, ym = yy - 1;
					int xp = xx + 1, xm = xx - 1;

					
					if (counterd_buffer[xx + xs*(yy+ys*zz)] >= 0) {
						// meeting_plane.insert(IPoi3<int>(xx, yy, zz));
						/*if (!meeting_plane_positions[zz][yy][xx] &&
							((distance[zz][yp][xx] >= 0 && !meeting_plane_positions[zz][yp][xx])
							|| (distance[zz][ym][xx] >= 0 && !meeting_plane_positions[zz][ym][xx])
							|| (distance[zz][yy][xp] >= 0 && !meeting_plane_positions[zz][yy][xp])
							|| (distance[zz][yy][xm] >= 0 && !meeting_plane_positions[zz][yy][xm])
							|| (distance[zp][yy][xx] >= 0 && !meeting_plane_positions[zp][yy][xx])
							|| (distance[zm][yy][xx] >= 0 && !meeting_plane_positions[zm][yy][xx]))
							&& xm > 2 && xp < xs - 2 && ym > 2 && yp < ys - 2 && zm > 2 && zp < zs - 2) {
							meeting_plane_positions[zz][yy][xx] = 1;
						}*/
						/*if (distance[zz][yy][xx] == counterd[zz][yy][xx]) {
							meeting_plane_positions[zz][yy][xx] = 1;
						}*/
						if (distance_buffer[xx + xs * (yy + ys * zz)] - counterd_buffer[xx + xs * (yy + ys * zz)] >= 0 &&
							((distance_buffer[xx + xs * (yp + ys * zz)] > 0 && distance_buffer[xx + xs * (yp + ys * zz)] - counterd_buffer[xx + xs * (yp + ys * zz)] < 0) ||
								(distance_buffer[xx + xs * (ym + ys * zz)] > 0 && distance_buffer[xx + xs * (ym + ys * zz)] - counterd_buffer[xx + xs * (ym + ys * zz)] < 0) ||
								(distance_buffer[xp + xs * (yy + ys * zz)] > 0 && distance_buffer[xp + xs * (yy + ys * zz)] - counterd_buffer[xp + xs * (yy + ys * zz)] < 0) ||
								(distance_buffer[xm + xs * (yy + ys * zz)] > 0 && distance_buffer[xm + xs * (yy + ys * zz)] - counterd_buffer[xm + xs * (yy + ys * zz)] < 0) ||
								(distance_buffer[xx + xs * (yy + ys * zp)] > 0 && distance_buffer[xx + xs * (yy + ys * zp)] - counterd_buffer[xx + xs * (yy + ys * zp)] < 0) ||
								(distance_buffer[xx + xs * (yy + ys * zm)] > 0 && distance_buffer[xx + xs * (yy + ys * zm)] - counterd_buffer[xx + xs * (yy + ys * zm)] < 0))) {
							meeting_plane_buffer[xx + xs * (yy + ys * zz)] = 1;
						}

						//continue;
					}
				}
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------


void CPhaseContainer::Initialize(SVoxImg<SWorkImg<realnum>>& data, CVec3& start_point, CVec3& end_point) {
	unsigned int spacex(data.xs), spacey(data.ys), spacez(data.zs);
	//g_cyc = 0;
	m_thickstate = sitk::Image({ spacex, spacey, spacez }, sitk::sitkInt32) - 1;
	m_Sumcurvature = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64);
	meeting_plane_positions = sitk::Image({ spacex, spacey, spacez }, sitk::sitkInt32);
	// expansion
	unx = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64);
	uny = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64);
	unz = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64);
	m_smoothstate = sitk::Image({ spacex, spacey, spacez }, sitk::sitkInt32) - 1;
	// expansion


	m_field[0] = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64) - 1;
	m_field[1] = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64) - 1;
	m_distance[0] = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64) - 1;
	m_distance[1] = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64) - 1;
	m_combined_distance = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64) - 1;
	m_flow_idx = sitk::Image({ spacex, spacey, spacez }, sitk::sitkInt32) - 1;


	m_velo[0] = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64);
	m_velo[1] = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64);
	m_aux = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64);
	m_smoothdist = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64) - 1;
	m_smoothaux = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64);
	m_smoothaux2 = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64);

	m_data = vox_img_2_sitk(data);
	m_sample_image = sitk::Image({ spacex, spacey, spacez }, sitk::sitkFloat64);

	int hs = 11 - 0;//5 (int)(1.5f*g_w/2);
	m_currentdistance = 0;
	for (int ii = 0; ii < 2; ++ii) {
		CVec3 point = !ii ? start_point : end_point;
		double* field_buffer = m_field[ii].GetBufferAsDouble();
		double* distance_buffer = m_distance[ii].GetBufferAsDouble();
		for (int zz = point.z - hs; zz < point.z + hs; ++zz) {
			if (zz < 0) continue;
			for (int yy = point.y - hs; yy < point.y + hs; ++yy) {
				if (yy < 0) continue;
				for (int xx = point.x - hs; xx < point.x + hs; ++xx) {
					if (xx < 0) continue;
					int dx = xx - point.x, dy = yy - point.y;
					int dz = zz - point.z;
					realnum dd = (realnum)(dx * dx + dy * dy + dz * dz);
					if ((int)dd < 10 * 10 / 1) {
						//TODO: initialize the two starting points on different distance maps (instead of [0]: [ii])
						field_buffer[xx + spacex * (yy + spacey * zz)] = 1.0; // ii->0
						dd = sqrt(dd);
						distance_buffer[xx + spacex * (yy + spacey * zz)] = dd; // ii->0
						if (m_currentdistance < dd) m_currentdistance = dd;

						// curvatures
						//if (dd < 1e-11) dd = 1e-11; dd = 1.0/dd;
						//m_phasefield.m_thickstate[ii][zz][yy][xx] = 2;
						//m_phasefield.m_Sumcurvature[ii][zz][yy][xx] = dd+dd;
					}

				}
			}
		}

	}


	m_bdone = false;
}

void CPhaseContainer::Initialize(CPhaseContainer& phasefield, vector<double>& rotation_matrix, bool inverse) {
	vector<unsigned int> sample_size;
	resample_img(phasefield.m_thickstate, m_thickstate, rotation_matrix, sample_size, inverse);
	resample_img(phasefield.m_Sumcurvature, m_Sumcurvature, rotation_matrix, sample_size, inverse);
	neg_to_minus1(m_Sumcurvature);

	resample_img<sitk::sitkNearestNeighbor>(phasefield.meeting_plane_positions, meeting_plane_positions, rotation_matrix, sample_size, inverse);
	// expansion

	resample_img(phasefield.unx, unx, rotation_matrix, sample_size, inverse);
	neg_to_minus1(unx);
	resample_img(phasefield.uny, uny, rotation_matrix, sample_size, inverse);
	neg_to_minus1(uny);
	resample_img(phasefield.unz, unz, rotation_matrix, sample_size, inverse);
	neg_to_minus1(unz);
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

	
	m_bdone = phasefield.m_bdone;
}

void CPhaseContainer::SmoothMap(sitk::Image &src1, sitk::Image& src2, sitk::Image &out)
{
	int j = 0;
	// int j = i;
	//SVoxImg<SWorkImg<realnum>> &src1 = m_phasefield.m_distance[i]; // source
	//SVoxImg<SWorkImg<realnum>>& src2 = m_phasefield.m_distance[(i+1)%2]; // source
	double* aux2_buffer = m_smoothaux2.GetBufferAsDouble();
	double* aux_buffer = m_smoothaux.GetBufferAsDouble();
	//SVoxImg<SWorkImg<realnum>> &out = m_phasefield.m_smoothdist[j]; // output
	
	int* smstat_buffer = m_smoothstate.GetBufferAsInt32();

	vector<unsigned> size = src1.GetSize();
	int xs = size[0], ys = size[1], zs = size[2];

	double* src1_buffer = src1.GetBufferAsDouble();
	double* src2_buffer = src2.GetBufferAsDouble();
	double* out_buffer = out.GetBufferAsDouble();

	#pragma omp parallel for schedule(static)
	/*for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {*/
	for (int zyx = 0; zyx < zs * ys * xs; zyx++) {
		int zz = zyx / (ys * xs);
		{
			int yy = (zyx / xs) % ys;
			{
				int xx = zyx % xs;
				int smst = smstat_buffer[xx+xs*(yy+ys*zz)];
				if (smst > 0) continue;
				smstat_buffer[xx + xs * (yy + ys * zz)] = -2;
				if (src1_buffer[xx + xs * (yy + ys * zz)] < 0 && src2_buffer[xx + xs * (yy + ys * zz)] < 0) {
					aux_buffer[xx + xs * (yy + ys * zz)] = src1_buffer[xx + xs * (yy + ys * zz)];
					continue;
				}
				int xm(xx-1); if (xm == -1) xm = 1;
				int xp(xx+1); if (xp == xs) xp = xs-2;
				realnum s  = src1_buffer[xx + xs * (yy + ys * zz)] >= 0 ? src1_buffer[xx + xs * (yy + ys * zz)] : src2_buffer[xx + xs * (yy + ys * zz)];
				int w = 0;
				realnum comp = src1_buffer[xp + xs * (yy + ys * zz)] >= 0 ? src1_buffer[xp + xs * (yy + ys * zz)] : src2_buffer[xp + xs * (yy + ys * zz)];
				if (comp >= 0) { s += comp; ++w; }
				comp = src1_buffer[xm + xs * (yy + ys * zz)] >= 0 ? src1_buffer[xm + xs * (yy + ys * zz)] : src2_buffer[xm + xs * (yy + ys * zz)];
				if (comp >= 0) { s += comp; ++w; }
				if (!w) aux_buffer[xx + xs * (yy + ys * zz)] = s;
				else {
					s += src1_buffer[xx + xs * (yy + ys * zz)] >= 0 ? src1_buffer[xx + xs * (yy + ys * zz)] : src2_buffer[xx + xs * (yy + ys * zz)];
					if (w == 1) aux_buffer[xx + xs * (yy + ys * zz)] = (1.0/3.0)*s;
					else {
						aux_buffer[xx + xs * (yy + ys * zz)] = 0.25*s;
						++smstat_buffer[xx + xs * (yy + ys * zz)]; //
					}
				}
			}
		}
	}
	#pragma omp parallel for schedule(static)
	/*for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {*/
	for (int zyx = 0; zyx < zs * ys * xs; zyx++) {
		int zz = zyx / (ys * xs);
		{
			int yy = (zyx / xs) % ys;
			{
				int xx = zyx % xs;
				int ym(yy - 1); if (ym == -1) ym = 1;
				int yp(yy + 1); if (yp == ys) yp = ys - 2;
				int smst = smstat_buffer[xx + xs * (yy + ys * zz)];
				if (smst > 0) continue; 
				if (aux_buffer[xx + xs * (yy + ys * zz)] < 0) {
					aux2_buffer[xx + xs * (yy + ys * zz)] = aux_buffer[xx + xs * (yy + ys * zz)];
					continue;
				}
				realnum s  = aux_buffer[xx + xs * (yy + ys * zz)]; int w = 0;
				realnum comp = aux_buffer[xx + xs * (yp + ys * zz)];
				if (comp >= 0) { s += comp; ++w; }
				comp = aux_buffer[xx + xs * (ym + ys * zz)];
				if (comp >= 0) { s += comp; ++w; }
				if (!w) aux2_buffer[xx + xs * (yy + ys * zz)] = s;
				else {
					s += aux_buffer[xx + xs * (yy + ys * zz)];
					if (w == 1) aux2_buffer[xx + xs * (yy + ys * zz)] = (1.0/3.0)*s;
					else {
						aux2_buffer[xx + xs * (yy + ys * zz)] = 0.25*s;
						++smstat_buffer[xx + xs * (yy + ys * zz)]; //
					}
				}
			}
		}
	}
	#pragma omp parallel for schedule(static)
	/*for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {*/
	for (int zyx = 0; zyx < zs * ys * xs; zyx++) {
		int zz = zyx / (ys * xs);
		{
			int yy = (zyx / xs) % ys;
			{
				int xx = zyx % xs;
				int zm(zz - 1); if (zm == -1) zm = 1;
				int zp(zz + 1); if (zp == zs) zp = zs - 2;
				int smst = smstat_buffer[xx + xs * (yy + ys * zz)];
				if (smst > 0) continue; 
				if (aux2_buffer[xx + xs * (yy + ys * zz)] < 0) {
					out_buffer[xx + xs * (yy + ys * zz)] = aux2_buffer[xx + xs * (yy + ys * zz)];
					continue;
				}
				realnum s  = aux2_buffer[xx + xs * (yy + ys * zz)]; int w = 0;
				realnum comp = aux2_buffer[xx + xs * (yy + ys * zp)];
				if (comp >= 0) { s += comp; ++w; }
				comp = aux2_buffer[xx + xs * (yy + ys * zm)];
				if (comp >= 0) { s += comp; ++w; }
				if (!w) out_buffer[xx + xs * (yy + ys * zz)] = s;
				else {
					s += aux2_buffer[xx + xs * (yy + ys * zz)];
					if (w == 1) out_buffer[xx + xs * (yy + ys * zz)] = (1.0/3.0)*s;
					else {
						out_buffer[xx + xs * (yy + ys * zz)] = 0.25*s;
						++smstat_buffer[xx + xs * (yy + ys * zz)]; //
					}
				}
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------

void CPhaseContainer::CalculateFundQuant(int i, int test)
{
	int j = 0;
	// int j = i;
	int* thstat_buffer = m_thickstate.GetBufferAsInt32();
	double* Sum_buffer = m_Sumcurvature.GetBufferAsDouble();
	double* distance_buffer = m_smoothdist.GetBufferAsDouble(); // STAYS
	double* unx_buffer = unx.GetBufferAsDouble();
	double* uny_buffer = uny.GetBufferAsDouble();
	double* unz_buffer = unz.GetBufferAsDouble();
	SmoothMap(m_distance[i], m_distance[(i + 1) % 2], m_smoothdist); // m_distance[i];
	int* smstat_buffer = m_smoothstate.GetBufferAsInt32();
	vector<unsigned> size = m_smoothdist.GetSize();
	int xs = size[0], ys = size[1], zs = size[2];

	////////////////////////////////////////////////////////////////////////////////
	// curvature update
	////////////////////////////////////////////////////////////////////////////////
	static int nsh; nsh = 0;
	//int ign(0), ian(0);

	#pragma omp parallel for schedule(static)
	/*for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {*/
	for (int zyx = 0; zyx < zs * ys * xs; zyx++) {
		int zz = zyx / (ys * xs);
		{
			int yy = (zyx / xs) % ys;
			{
				int xx = zyx % xs;
				int zm(zz - 1); if (zm == -1) zm = 1;
				int zp(zz + 1); if (zp == zs) zp = zs - 2;
				int yp(yy + 1); if (yp == ys) yp = ys - 2;
				int ym(yy - 1); if (ym == -1) ym = 1; 
				int xp(xx+1); if (xp == xs) xp = xs-2;
				int xm(xx-1); if (xm == -1) xm = 1;

				if (thstat_buffer[xx + xs * (yy + ys * zz)] == 2) continue;

				//Gau[zz][yy][xx] = 0; Sum[zz][yy][xx] = 0;

				if (distance_buffer[xx + xs * (yy + ys * zz)] >= 0
				 &&	distance_buffer[xx + xs * (yp + ys * zz)] >= 0 && distance_buffer[xx + xs * (ym + ys * zz)] >= 0
				 && distance_buffer[xp + xs * (yy + ys * zz)] >= 0 && distance_buffer[xm + xs * (yy + ys * zz)] >= 0
				 && distance_buffer[xx + xs * (yy + ys * zp)] >= 0 && distance_buffer[xx + xs * (yy + ys * zm)] >= 0 // sum curvature can be calculated

				 && distance_buffer[xp + xs * (yp + ys * zz)] >= 0 && distance_buffer[xm + xs * (ym + ys * zz)] >= 0
				 && distance_buffer[xm + xs * (yp + ys * zz)] >= 0 && distance_buffer[xp + xs * (ym + ys * zz)] >= 0
				 &&	distance_buffer[xp + xs * (yy + ys * zp)] >= 0 && distance_buffer[xm + xs * (yy + ys * zm)] >= 0
				 && distance_buffer[xm + xs * (yy + ys * zp)] >= 0 && distance_buffer[xp + xs * (yy + ys * zm)] >= 0
				 &&	distance_buffer[xx + xs * (yp + ys * zp)] >= 0 && distance_buffer[xx + xs * (ym + ys * zm)] >= 0
				 && distance_buffer[xx + xs * (yp + ys * zm)] >= 0 && distance_buffer[xx + xs * (ym + ys * zp)] >= 0 // Gauss curvature can be calculated as well
				 ) {
					// calculating the gradients
					realnum gxp = (distance_buffer[xp + xs * (yy + ys * zz)] -distance_buffer[xx + xs * (yy + ys * zz)]);
					realnum gyp = (distance_buffer[xx + xs * (yp + ys * zz)] -distance_buffer[xx + xs * (yy + ys * zz)]);
					realnum gzp = (distance_buffer[xx + xs * (yy + ys * zp)] -distance_buffer[xx + xs * (yy + ys * zz)]);
					realnum gxm = (distance_buffer[xx + xs * (yy + ys * zz)] -distance_buffer[xm + xs * (yy + ys * zz)]);
					realnum gym = (distance_buffer[xx + xs * (yy + ys * zz)] -distance_buffer[xx + xs * (ym + ys * zz)]);
					realnum gzm = (distance_buffer[xx + xs * (yy + ys * zz)] -distance_buffer[xx + xs * (yy + ys * zm)]);

					realnum gx = 0.5*(gxp+gxm);
					realnum gy = 0.5*(gyp+gym);
					realnum gz = 0.5*(gzp+gzm);
					realnum hxx = (distance_buffer[xp + xs * (yy + ys * zz)] +distance_buffer[xm + xs * (yy + ys * zz)] -2*distance_buffer[xx + xs * (yy + ys * zz)]);
					realnum hyy = (distance_buffer[xx + xs * (yp + ys * zz)] +distance_buffer[xx + xs * (ym + ys * zz)] -2*distance_buffer[xx + xs * (yy + ys * zz)]);
					realnum hzz = (distance_buffer[xx + xs * (yy + ys * zp)] +distance_buffer[xx + xs * (yy + ys * zm)] -2*distance_buffer[xx + xs * (yy + ys * zz)]);
					realnum hxy = 0.25*(distance_buffer[xp + xs * (yp + ys * zz)] +distance_buffer[xm + xs * (ym + ys * zz)] -distance_buffer[xm + xs * (yp + ys * zz)] -distance_buffer[xp + xs * (ym + ys * zz)]);
					realnum hxz = 0.25*(distance_buffer[xp + xs * (yy + ys * zp)] +distance_buffer[xm + xs * (yy + ys * zm)] -distance_buffer[xm + xs * (yy + ys * zp)] -distance_buffer[xp + xs * (yy + ys * zm)]);
					realnum hyz = 0.25*(distance_buffer[xx + xs * (yp + ys * zp)] +distance_buffer[xx + xs * (ym + ys * zm)] -distance_buffer[xx + xs * (ym + ys * zp)] -distance_buffer[xx + xs * (yp + ys * zm)]);
					
					realnum ig2 = 1.0/(gx*gx+gy*gy+gz*gz+1e-99);
					realnum ig = sqrt(ig2);
					gx *= ig; gy *= ig; gz *= ig; // unit here
					unx_buffer[xx + xs * (yy + ys * zz)] = gx;
					uny_buffer[xx + xs * (yy + ys * zz)] = gy;
					unz_buffer[xx + xs * (yy + ys * zz)] = gz;

					/*realnum cr1 = gx * gy * (hxz * hyz - hxy * hzz) + gx * gz * (hxy * hyz - hxz * hyy) + gy * gz * (hxy * hxz - hyz * hxx); cr1 *= ig2;
					realnum cr2 = gx*gx*(hyy*hzz-hyz*hyz)+gy*gy*(hxx*hzz-hxz*hxz)+gz*gz*(hxx*hyy-hxy*hxy); cr2 *= ig2;
					Gau[zz][yy][xx] = 2*cr1+cr2;*/

					realnum cr1 = gx*gx*(hyy+hzz)+gy*gy*(hxx+hzz)+gz*gz*(hxx+hyy); cr1 *= ig;
					realnum cr2 = gx*gy*hxy+gx*gz*hxz+gy*gz*hyz; cr2 *= ig;
					Sum_buffer[xx + xs * (yy + ys * zz)] = cr1-2*cr2;

					thstat_buffer[xx + xs * (yy + ys * zz)] = 1; // calculated

					if (smstat_buffer[xx + xs * (yp + ys * zz)] == 1 && smstat_buffer[xx + xs * (ym + ys * zz)] == 1
					 && smstat_buffer[xp + xs * (yy + ys * zz)] == 1 && smstat_buffer[xm + xs * (yy + ys * zz)] == 1
					 && smstat_buffer[xx + xs * (yy + ys * zp)] == 1 && smstat_buffer[xx + xs * (yy + ys * zm)] == 1
						
						&& smstat_buffer[xp + xs * (yp + ys * zz)] == 1 && smstat_buffer[xm + xs * (ym + ys * zz)] == 1
						&& smstat_buffer[xm + xs * (yp + ys * zz)] == 1 && smstat_buffer[xp + xs * (ym + ys * zz)] == 1
						&& smstat_buffer[xp + xs * (yy + ys * zp)] == 1 && smstat_buffer[xm + xs * (yy + ys * zm)] == 1
						&& smstat_buffer[xm + xs * (yy + ys * zp)] == 1 && smstat_buffer[xp + xs * (yy + ys * zm)] == 1
						&& smstat_buffer[xx + xs * (yp + ys * zp)] == 1 && smstat_buffer[xx + xs * (ym + ys * zm)] == 1
						&& smstat_buffer[xx + xs * (yp + ys * zm)] == 1 && smstat_buffer[xx + xs * (ym + ys * zp)] == 1

					 && smstat_buffer[xx + xs * (yy + ys * zz)] == 1)
						thstat_buffer[xx + xs * (yy + ys * zz)] = 2;/**/

				}

			}
		}
	}



	//SmoothMap(Gau,i); // Gauss
	//SmoothMap(Sum,i); // Sum

	//SmoothMap(thick,i); // thickness


}

bool g_modeswitch = false;

realnum CPhaseContainer::UpdateVelo(int i, bool use_correction) {
	double* field_buffer = m_field[i].GetBufferAsDouble();

	double* counterfield_buffer = m_field[(i + 1) % 2].GetBufferAsDouble();

	// m_velo
	sitk::Image& velo = m_velo[i];

	// m_distance[i]
	double* distance_buffer = m_distance[i].GetBufferAsDouble();

	// m_distance[(i+1)&1]
	double* counterd_buffer = m_distance[(i + 1) & 1].GetBufferAsDouble();


	// m_Sumcurvature[i]
	double* Sum_buffer = m_Sumcurvature.GetBufferAsDouble();
	int* thickstate_buffer = m_thickstate.GetBufferAsInt32();
	double* data_buffer = m_data.GetBufferAsDouble();
	double* unx_buffer = unx.GetBufferAsDouble();
	double* uny_buffer = uny.GetBufferAsDouble();
	double* unz_buffer = unz.GetBufferAsDouble();

	vector<unsigned> size = m_field[i].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];

	////////////////////////////////////////////////////////////////////
	// fundamental quantities
	////////////////////////////////////////////////////////////////////

	CalculateFundQuant(i); // for object-level evolution

	////////////////////////////////////////////////////////////////////
	// evolution logic
	////////////////////////////////////////////////////////////////////
	typedef struct {
		double val;
		char pad[128];
	} tvals;
	tvals maxinfo[MAX_THREADS];
	realnum maxv(0);
	unsigned long counts[MAX_THREADS];
	for (int i = 0; i < MAX_THREADS; i++) {
		maxinfo[i].val = 0;
		counts[i] = 0;

	}
	// Inhomogeneous
	// Calculate the velocity at every point (velo[zz][yy][xx]), and the maximum velocity (maxv)
	velo = velo - velo;
	double* velo_buffer = velo.GetBufferAsDouble();
	{
#pragma omp parallel for shared(maxinfo) schedule(static)
		/*for (int zyx = 0; zyx < zs * ys * xs; zyx++) {
			int zz = zyx / (ys * xs);
			if (zz == 0 || zz >= zs - 1) continue;
			{
				int yy = (zyx / xs) % ys;
				if (yy == 0 || yy >= ys - 1) continue;
				{
					int xx = zyx % xs;
					if (xx == 0 || xx >= xs - 1) continue;*/
		for (int zz = 1; zz < zs - 1; ++zz) {
			for (int yy = 1; yy < ys - 1; ++yy) {
				for (int xx = 1; xx < xs - 1; ++xx) {
					int id = omp_get_thread_num();
					counts[id]++;
					// velo[zz][yy][xx] = 0; // zero out
					/*IPoi3<int>* pos = new IPoi3(xx, yy, zz);
					bool should_run = true;
					#pragma omp critical(meeting_plane)
					if (meeting_plane.count(*pos)) {
						delete pos;
						should_run = false;
					}
					if (!should_run) continue;
					delete pos;*/
					/*
					if (meeting_plane_positions[zz][yy][xx]) {
						continue;
					}*/
					int zp = zz + 1, zm = zz - 1;
					int yp = yy + 1, ym = yy - 1;
					int xp = xx + 1, xm = xx - 1;

					if (distance_buffer[xx + xs * (yp + ys * zz)] < 0 && distance_buffer[xx + xs * (ym + ys * zz)] < 0
						&& distance_buffer[xp + xs * (yy + ys * zz)] < 0 && distance_buffer[xm + xs * (yy + ys * zz)] < 0
						&& distance_buffer[xx + xs * (yy + ys * zp)] < 0 && distance_buffer[xx + xs * (yy + ys * zm)] < 0) continue;



					
					if (field_buffer[xx + xs * (yy + ys * zz)] >= 0.95) continue;

					realnum xn = field_buffer[xp + xs * (yy + ys * zz)] - field_buffer[xm + xs * (yy + ys * zz)];
					realnum yn = field_buffer[xx + xs * (yp + ys * zz)] - field_buffer[xx + xs * (ym + ys * zz)];
					realnum zn = field_buffer[xx + xs * (yy + ys * zp)] - field_buffer[xx + xs * (yy + ys * zm)];
					realnum gradlen = sqrt(xn * xn + yn * yn + zn * zn + 1e-99);


					realnum sumcur(0.0);

					//Sum[zz][yy][xx] = 0;
					realnum wAct(0.0);
					int iok = -1;
					for (int rr = 1; rr <= 25; ++rr) {
						for (int iz = -rr; iz <= rr; ++iz) {
							for (int iy = -rr; iy <= rr; ++iy) {
								for (int ix = -rr; ix <= rr; ++ix) {
									if (!ix && !iy && !iz) continue;
									realnum r2 = (realnum)ix * ix + (realnum)iy * iy + (realnum)iz * iz;

									if (r2 <= rr * rr) {
										int pz = zz + iz; if (pz < 1) pz = 1; if (pz > zs - 2) pz = zs - 2;
										int py = yy + iy; if (py < 1) py = 1; if (py > ys - 2) py = ys - 2;
										int px = xx + ix; if (px < 1) px = 1; if (px > xs - 2) px = xs - 2;
										if (thickstate_buffer[px + xs * (py + ys * pz)] < 0) continue;
										realnum dot = unx_buffer[px + xs * (py + ys * pz)] * ix + uny_buffer[px + xs * (py + ys * pz)] * iy + unz_buffer[px + xs * (py + ys * pz)] * iz;
										dot *= -1;
										if (dot < 0.82) continue;
										++iok; if (dot > 0.94) ++iok;
										wAct += dot / r2;
										sumcur += (dot / r2) * Sum_buffer[px + xs * (py + ys * pz)];
									}


								} // for ix
							} // for iy
						}// for iz
						if (wAct > 0 && iok > 0) {
							/*Gau[zz][yy][xx] = gcw / wAct;
							Sum[zz][yy][xx] = sumcur/wAct;
							++g_rr[rr];*/
							sumcur /= wAct;
							break;
						}

					} // for rr

					realnum eikon = data_buffer[xx + xs * (yy + ys * zz)];

					realnum d0 = data_buffer[xx + xs * (yy + ys * zz)];
					if (d0 < 1e-33) d0 = 1e-33;

					realnum squar = -1.0;
					realnum discr = 1 + 4 * sumcur * d0 * 1;// *gradlen;

					if (abs(sumcur) > 1e-33) {
						squar = -1.0 + sqrt(discr);
						squar /= 2 * sumcur;
					}
					else {
						squar = d0;
					}

					if (use_correction) {
						//realnum den = 1.0 + sumcur/1.2; if (den < 1e-6) den = 1e-6;	eikon = data[zz][yy][xx]/den;
						if (squar > 0) {
							eikon = squar;
							//eikon *= exp(-sumcur); // best
						}
						else {
							eikon = data_buffer[xx + xs * (yy + ys * zz)];
						}
					}
					else
						eikon = data_buffer[xx + xs * (yy + ys * zz)];

					eikon *= gradlen; // normalize
					if (eikon < 1e-11) eikon = 1e-11;
					velo_buffer[xx + xs * (yy + ys * zz)] = eikon; // dS = 1

					/*#pragma omp critical(meeting_plane)
					{
						if (eikon > maxv) maxv = eikon;
						if (counterd[zz][yy][xx] >= 0) {
							// meeting_plane.insert(IPoi3<int>(xx, yy, zz));
							meeting_plane_positions[zz][yy][xx] = 1;
						}
					}*/
					if (eikon > maxinfo[id].val) maxinfo[id].val = eikon;
					/*if (counterd[zz][yy][xx] >= 0) {
						// meeting_plane.insert(IPoi3<int>(xx, yy, zz));
						meeting_plane_positions[zz][yy][xx] = 1;
					}*/

					//loop body
				} //for xx
			} // for yy
		} // for zz
		for (int i = 0; i < MAX_THREADS; i++) {
			if (maxinfo[i].val > maxv)
				maxv = maxinfo[i].val;
		}
		for (int i = 0; i < MAX_THREADS; i++)
			counts[i] = 0;
	}



	////////////////////////////////////////////////////////////////////
	// conditioning (speed limitation)
	////////////////////////////////////////////////////////////////////
	return maxv;
}
void CPhaseContainer::UpdateField(int i, realnum maxv) {
	double* field_buffer = m_field[i].GetBufferAsDouble();
	double* velo_buffer = m_velo[i].GetBufferAsDouble();
	int* flow_idx_buffer = m_flow_idx.GetBufferAsInt32();
	vector<unsigned> size = m_field[i].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
#pragma omp parallel for schedule(static)
	for (int zyx = 0; zyx < zs * ys * xs; zyx++) {
		int zz = zyx / (ys * xs);
		if (zz == 0 || zz >= zs - 1) continue;
		{
			int yy = (zyx / xs) % ys;
			if (yy == 0 || yy >= ys - 1) continue;
			{
				int xx = zyx % xs;
				if (xx == 0 || xx >= xs - 1) continue;
	/*for (int zz = 1; zz < zs - 1; ++zz) {
		for (int yy = 0 + 1; yy < ys - 1; ++yy) {
			for (int xx = 0 + 1; xx < xs - 1; ++xx) {*/
				field_buffer[xx + xs * (yy + ys * zz)] += velo_buffer[xx + xs * (yy + ys * zz)] * maxv;// ;
				if (field_buffer[xx + xs * (yy + ys * zz)] > 0) {
					if (field_buffer[xx + xs * (yy + ys * zz)] > 1.25) field_buffer[xx + xs * (yy + ys * zz)] = 1.25;
					if (flow_idx_buffer[xx + xs * (yy + ys * zz)] < 0) {
						flow_idx_buffer[xx + xs * (yy + ys * zz)] = i;
					}
				}

			}
		}
	}
}
void CPhaseContainer::UpdateDistance(int i, realnum current_distance) {
	double* field_buffer = m_field[i].GetBufferAsDouble();
	// m_phasefield.m_distance[i]
	double* distance_buffer = m_distance[i].GetBufferAsDouble();

	// m_phasefield.m_distance[(i+1)&1]
	double* counterd_buffer = m_distance[(i + 1) & 1].GetBufferAsDouble();
	vector<unsigned> size = m_field[i].GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
#pragma omp parallel for schedule(static)
	for (int zyx = 0; zyx < zs * ys * xs; zyx++) {
			int zz = zyx / (ys * xs);
			if (zz == 0 || zz >= zs - 1) continue;
			{
				int yy = (zyx / xs) % ys;
				if (yy == 0 || yy >= ys - 1) continue;
				{
					int xx = zyx % xs;
					if (xx == 0 || xx >= xs - 1) continue;
	/*for (int zz = 1; zz < zs - 1; ++zz) {
		for (int yy = 0 + 1; yy < ys - 1; ++yy) {
			for (int xx = 0 + 1; xx < xs - 1; ++xx) {*/
				if (field_buffer[xx + xs * (yy + ys * zz)] > 0 && distance_buffer[xx + xs * (yy + ys * zz)] < -0.5f) {
					distance_buffer[xx + xs * (yy + ys * zz)] = current_distance;
				}
				if (distance_buffer[xx + xs * (yy + ys * zz)] < -0.5f && counterd_buffer[xx + xs * (yy + ys * zz)] < -0.5f) {
					m_bdone = false;
				}
			}
		}
	}
}
void CPhaseContainer::Iterate(bool use_correction)
{
	// m_phasefield.m_field[i]
	realnum maxv0 = UpdateVelo(0, use_correction);
	realnum maxv1 = UpdateVelo(1, use_correction);
	realnum maxv = maxv0 > maxv1? maxv0 : maxv1;
	if (maxv < 1e-11) maxv = 1e-11;

	realnum mdatspedmax = 0.5 * 4; //*2 set max speed
	maxv = mdatspedmax / maxv;
	m_currentdistance += maxv;
	//m_currentdistance += 1;

	// update phase field values
	UpdateField(0, maxv);
	UpdateField(1, maxv);
	m_bdone = true;
	UpdateDistance(0, m_currentdistance);
	UpdateDistance(1, m_currentdistance);
	// TODO: check for collision
	// updating the distance map to the current value, and checking if it is complete
	
	////////////////////////////////////////////////////////////////////
	// phase field regularization w/o velocity zeroed out
	////////////////////////////////////////////////////////////////////


}
void CPhaseContainer::CalculateAlignedCombinedDistance(double p0_x, double p1_x) {
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
	//TODO: instead of searching for the selected x slice, locate the meeting surface of the two processes, and use that
	int xs = distance1.xs, ys = distance1.ys, zs = distance1.zs;
	m_nall = m_npos = 0;
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
	m_nfct = 2000;
	m_bInited = true;

}

void CPlanePhaseField::GetDistancemean(sitk::Image& distance, int xslice)
{
	//TODO: instead of searching for the selected x slice, locate the meeting surface of the two processes, and use that
	vector<unsigned> size = distance.GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	double* distance_buffer = distance.GetBufferAsDouble();
	m_nall = m_npos = 0;
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
	m_nfct = 2000;
	m_bInited = true;

}

void CPlanePhaseField::Iterate() 
{
	if (!m_bInited) return;

	int xs = m_field2D.xs, ys = m_field2D.ys;
	m_velo2D.GetLaplace(m_field2D);
	m_velo2D *= 4.0f;

	realnum fac = 1.0f;
	realnum maxv(0.0f);
	#pragma omp parallel for
	for (int yy = 1; yy < ys - 1; ++yy) {
		for (int xx = 1; xx < xs - 1; ++xx) {
				realnum fmm = m_field2D[yy][xx];
				m_velo2D[yy][xx] -= fac * fmm * fmm * fmm;
				m_velo2D[yy][xx] += fac * fmm;
				realnum gfx = 0.5f * (m_field2D[yy][xx + 1]- m_field2D[yy][xx - 1]);
				realnum gfy = 0.5f * (m_field2D[yy + 1][xx] - m_field2D[yy - 1][xx]);

				m_velo2D[yy][xx] += 10.0f * (gfx * m_gx2D[yy][xx] + gfy * m_gy2D[yy][xx]);
				realnum cv = m_velo2D[yy][xx];
				if (cv < 0) cv *= -1;
				if (cv > maxv) maxv = cv;
		}
	}

	if (maxv > 1e-11f)
		maxv = 0.025f/maxv;
	int npos(0);

	#pragma omp parallel for
	for (int yy = 1; yy < ys - 1; ++yy) {
		for (int xx = 1; xx < xs - 1; ++xx) {
			m_field2D[yy][xx] += m_velo2D[yy][xx]*maxv;
			if (m_field2D[yy][xx] > 0) ++npos;
		}
	}
	
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
