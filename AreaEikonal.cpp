#include "stdafx.h"
#include "AreaEikonal.h"
#include "math.h"
#include <omp.h>

#define BIG_NUMBER 1e11
#define MAX_THREADS 32
//const int max_threads = omp_get_max_threads();
int a = 0;
CCurvEikonal::CCurvEikonal(void)
{
	

}

CCurvEikonal::~CCurvEikonal(void)
{
}

// Image prep&check

#define GAUTOO 1

// Phasefield stuff

realnum g_w = 2.75;

void CCurvEikonal::PhaseInit(SVoxImg<SWorkImg<realnum>>& data, IPoi reginit, IPoi arrival, int zdeparture, int zarrival, int xSection)
{
	int spacex(data.xs), spacey(data.ys), spacez(data.zs);
	//g_cyc = 0;
	m_phasefield.m_thickstate.Set(spacex,spacey,spacez);
	m_phasefield.m_Sumcurvature.Set0(spacex,spacey,spacez);
	m_phasefield.meeting_plane_positions.Set0(spacex, spacey, spacez);
	// expansion
	m_phasefield.unx.Set0(spacex,spacey,spacez);
	m_phasefield.uny.Set0(spacex,spacey,spacez);
	m_phasefield.unz.Set0(spacex,spacey,spacez);
	m_phasefield.m_smoothstate.Set(spacex,spacey,spacez);
	// expansion


	m_phasefield.m_field[0].Set(spacex,spacey,spacez);
	m_phasefield.m_field[1].Set(spacex,spacey,spacez);
	m_phasefield.m_distance[0].Set(spacex,spacey,spacez);
	m_phasefield.m_distance[1].Set(spacex,spacey,spacez);
	m_phasefield.m_combined_distance.Set(spacex, spacey, spacez);

	m_phasefield.m_velo[0].Set0(spacex,spacey,spacez);
	m_phasefield.m_velo[1].Set0(spacex, spacey, spacez);
	m_phasefield.m_aux.Set0(spacex,spacey,spacez);
	m_phasefield.m_smoothdist.Set(spacex,spacey,spacez);
	m_phasefield.m_smoothaux.Set0(spacex,spacey,spacez);
	m_phasefield.m_smoothaux2.Set0(spacex,spacey,spacez);

	m_phasefield.m_data.Set0(spacex,spacey,spacez);
	m_phasefield.m_data = data;


	int hs = 11-0;//5 (int)(1.5f*g_w/2);
	m_phasefield.m_currentdistance =  0;
	for (int ii = 0; ii < 2; ++ii) {
		int &initz = !ii ? zdeparture:zarrival; // z parameter
		IPoi &rinit = !ii ? reginit:arrival; // x-y parameter

		for (int zz = initz-hs; zz < initz+hs; ++zz) {
			if (zz < 0) continue;
			for (int yy = rinit.y-hs; yy < rinit.y+hs; ++yy) {
				if (yy < 0) continue;
				for (int xx = rinit.x-hs; xx < rinit.x+hs; ++xx) {
					if (xx < 0) continue;
					int dx = xx-rinit.x, dy = yy-rinit.y;
					int dz = zz-initz;
					realnum dd = (realnum)(dx*dx+dy*dy+dz*dz);
					if ((int)dd < 10*10/1/*4 hs*hs/4*/) {
						//TODO: initialize the two starting points on different distance maps (instead of [0]: [ii])
						m_phasefield.m_field[ii][zz][yy][xx] = 1.0; // ii->0
						dd = sqrt(dd);
						m_phasefield.m_distance[ii][zz][yy][xx] = dd; // ii->0
						if (m_phasefield.m_currentdistance < dd) m_phasefield.m_currentdistance = dd;
						
						// curvatures
						//if (dd < 1e-11) dd = 1e-11; dd = 1.0/dd;
						//m_phasefield.m_thickstate[ii][zz][yy][xx] = 2;
						//m_phasefield.m_Sumcurvature[ii][zz][yy][xx] = dd+dd;
					}

				}
			}
		}

	}

	for (int ii = 0; ii < 10-10; ++ii) {
		RegularizePhaseField(m_phasefield.m_field[0],m_phasefield.m_velo[0]);
		RegularizePhaseField(m_phasefield.m_field[1],m_phasefield.m_velo[1]);
	}
	m_phasefield.m_aux[0] = 0;
	m_phasefield.m_aux[1] = 0;


	m_reference[0].x = (realnum)reginit.x;
	m_reference[0].y = (realnum)reginit.y;
	m_reference[0].z = (realnum)zdeparture;
	m_reference[1].x = (realnum)arrival.x;
	m_reference[1].y = (realnum)arrival.y;
	m_reference[1].z = (realnum)zarrival;

	m_distanceto = IPoi3<int>(arrival.x,arrival.y,zarrival);

	m_phasefield.m_bdone = false;
}

void CCurvEikonal::RegularizePhaseField(SVoxImg<SWorkImg<realnum>> &field, SVoxImg<SWorkImg<realnum>> &velo)
{
	//++g_cyc;

	int xs = field.xs, ys = field.ys, zs = field.zs;
	realnum fac = 21.0/(g_w*g_w), dfac = -(g_w*g_w)/16.0;

	SVoxImg<SWorkImg<realnum>> &aux = m_phasefield.m_aux;

	aux.GetLaplace(field);
	velo.GetLaplace(aux);
	velo *= dfac;

	#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				velo[zz][yy][xx] -= aux[zz][yy][xx];
			}
		}
	}

	#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				realnum fmm = field[zz][yy][xx];
				velo[zz][yy][xx] -= fac*fmm*fmm*fmm;
				velo[zz][yy][xx] += fac*fmm;
			}
		}
	}

	realnum totalfieldweight = 0.00175;
	#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				field[zz][yy][xx] += totalfieldweight*velo[zz][yy][xx];
				if (field[zz][yy][xx] > 2) field[zz][yy][xx] = 2; //
				else if (field[zz][yy][xx] < -2) field[zz][yy][xx] = -2; //
			}
		}
	}


}

//-------------------------------------------------------------------------------------------------


void CPhaseContainer::SmoothMap(SVoxImg<SWorkImg<realnum>> &src1, SVoxImg<SWorkImg<realnum>>& src2, SVoxImg<SWorkImg<realnum>> &out)
{
	int j = 0;
	// int j = i;
	//SVoxImg<SWorkImg<realnum>> &src1 = m_phasefield.m_distance[i]; // source
	//SVoxImg<SWorkImg<realnum>>& src2 = m_phasefield.m_distance[(i+1)%2]; // source
	SVoxImg<SWorkImg<realnum>> &aux2 = m_smoothaux2;
	SVoxImg<SWorkImg<realnum>> &aux = m_smoothaux;
	//SVoxImg<SWorkImg<realnum>> &out = m_phasefield.m_smoothdist[j]; // output
	
	SVoxImg<SWorkImg<int>> &smstat = m_smoothstate;

	int xs = src1.xs, ys = src1.ys, zs = src1.zs;

	#pragma omp parallel for collapse(3)
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				int smst = smstat[zz][yy][xx];
				if (smst > 0) continue;
				smstat[zz][yy][xx] = -2;
				if (src1[zz][yy][xx] < 0 && src2[zz][yy][xx] < 0) {
					aux[zz][yy][xx] = src1[zz][yy][xx];
					continue;
				}
				int xm(xx-1); if (xm == -1) xm = 1;
				int xp(xx+1); if (xp == xs) xp = xs-2;
				realnum s  = src1[zz][yy][xx] >= 0 ? src1[zz][yy][xx] : src2[zz][yy][xx];
				int w = 0;
				realnum comp = src1[zz][yy][xp] >= 0 ? src1[zz][yy][xp] : src2[zz][yy][xp];
				if (comp >= 0) { s += comp; ++w; }
				comp = src1[zz][yy][xm] >= 0 ? src1[zz][yy][xm] : src2[zz][yy][xm];
				if (comp >= 0) { s += comp; ++w; }
				if (!w) aux[zz][yy][xx] = s;
				else {
					s += src1[zz][yy][xx] >= 0 ? src1[zz][yy][xx] : src2[zz][yy][xx];
					if (w == 1) aux[zz][yy][xx] = (1.0/3.0)*s;
					else {
						aux[zz][yy][xx] = 0.25*s;
						++smstat[zz][yy][xx]; //
					}
				}
			}
		}
	}
	#pragma omp parallel for collapse(3)
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				int ym(yy - 1); if (ym == -1) ym = 1;
				int yp(yy + 1); if (yp == ys) yp = ys - 2;
				int smst = smstat[zz][yy][xx];
				if (smst > 0) continue; 
				if (aux[zz][yy][xx] < 0) {
					aux2[zz][yy][xx] = aux[zz][yy][xx];
					continue;
				}
				realnum s  = aux[zz][yy][xx]; int w = 0;
				realnum comp = aux[zz][yp][xx];
				if (comp >= 0) { s += comp; ++w; }
				comp = aux[zz][ym][xx];
				if (comp >= 0) { s += comp; ++w; }
				if (!w) aux2[zz][yy][xx] = s;
				else {
					s += aux[zz][yy][xx];
					if (w == 1) aux2[zz][yy][xx] = (1.0/3.0)*s;
					else {
						aux2[zz][yy][xx] = 0.25*s;
						++smstat[zz][yy][xx]; //
					}
				}
			}
		}
	}
	#pragma omp parallel for collapse(3)
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				int zm(zz - 1); if (zm == -1) zm = 1;
				int zp(zz + 1); if (zp == zs) zp = zs - 2;
				int smst = smstat[zz][yy][xx];
				if (smst > 0) continue; 
				if (aux2[zz][yy][xx] < 0) {
					out[zz][yy][xx] = aux2[zz][yy][xx];
					continue;
				}
				realnum s  = aux2[zz][yy][xx]; int w = 0;
				realnum comp = aux2[zp][yy][xx];
				if (comp >= 0) { s += comp; ++w; }
				comp = aux2[zm][yy][xx];
				if (comp >= 0) { s += comp; ++w; }
				if (!w) out[zz][yy][xx] = s;
				else {
					s += aux2[zz][yy][xx];
					if (w == 1) out[zz][yy][xx] = (1.0/3.0)*s;
					else {
						out[zz][yy][xx] = 0.25*s;
						++smstat[zz][yy][xx]; //
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
	SVoxImg<SWorkImg<int>> &thstat = m_thickstate;
	SVoxImg<SWorkImg<realnum>> &Sum = m_Sumcurvature;
	SVoxImg<SWorkImg<realnum>>& distance = m_smoothdist; // STAYS
	SmoothMap(m_distance[i], m_distance[(i + 1) % 2], distance); // m_distance[i];
	SVoxImg<SWorkImg<int>>& smstat = m_smoothstate;

	int xs = distance.xs, ys = distance.ys, zs = distance.zs;

	////////////////////////////////////////////////////////////////////////////////
	// curvature update
	////////////////////////////////////////////////////////////////////////////////
	static int nsh; nsh = 0;
	//int ign(0), ian(0);

	#pragma omp parallel for collapse(3)
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				int zm(zz - 1); if (zm == -1) zm = 1;
				int zp(zz + 1); if (zp == zs) zp = zs - 2;
				int yp(yy + 1); if (yp == ys) yp = ys - 2;
				int ym(yy - 1); if (ym == -1) ym = 1; 
				int xp(xx+1); if (xp == xs) xp = xs-2;
				int xm(xx-1); if (xm == -1) xm = 1;

				if (thstat[zz][yy][xx] == 2) continue;

				//Gau[zz][yy][xx] = 0; Sum[zz][yy][xx] = 0;

				if (distance[zz][yy][xx] >= 0
				 &&	distance[zz][yp][xx] >= 0 && distance[zz][ym][xx] >= 0 
				 && distance[zz][yy][xp] >= 0 && distance[zz][yy][xm] >= 0
				 && distance[zp][yy][xx] >= 0 && distance[zm][yy][xx] >= 0 // sum curvature can be calculated

				 && distance[zz][yp][xp] >= 0 && distance[zz][ym][xm] >= 0
				 && distance[zz][yp][xm] >= 0 && distance[zz][ym][xp] >= 0
				 &&	distance[zp][yy][xp] >= 0 && distance[zm][yy][xm] >= 0 
				 && distance[zp][yy][xm] >= 0 && distance[zm][yy][xp] >= 0
				 &&	distance[zp][yp][xx] >= 0 && distance[zm][ym][xx] >= 0 
				 && distance[zm][yp][xx] >= 0 && distance[zp][ym][xx] >= 0 // Gauss curvature can be calculated as well
				 ) {
					// calculating the gradients
					realnum gxp = (distance[zz][yy][xp]-distance[zz][yy][xx]);
					realnum gyp = (distance[zz][yp][xx]-distance[zz][yy][xx]);
					realnum gzp = (distance[zp][yy][xx]-distance[zz][yy][xx]);
					realnum gxm = (distance[zz][yy][xx]-distance[zz][yy][xm]);
					realnum gym = (distance[zz][yy][xx]-distance[zz][ym][xx]);
					realnum gzm = (distance[zz][yy][xx]-distance[zm][yy][xx]);

					realnum gx = 0.5*(gxp+gxm);
					realnum gy = 0.5*(gyp+gym);
					realnum gz = 0.5*(gzp+gzm);
					realnum hxx = (distance[zz][yy][xp]+distance[zz][yy][xm]-2*distance[zz][yy][xx]);
					realnum hyy = (distance[zz][yp][xx]+distance[zz][ym][xx]-2*distance[zz][yy][xx]);
					realnum hzz = (distance[zp][yy][xx]+distance[zm][yy][xx]-2*distance[zz][yy][xx]);
					realnum hxy = 0.25*(distance[zz][yp][xp]+distance[zz][ym][xm]-distance[zz][yp][xm]-distance[zz][ym][xp]);
					realnum hxz = 0.25*(distance[zp][yy][xp]+distance[zm][yy][xm]-distance[zp][yy][xm]-distance[zm][yy][xp]);
					realnum hyz = 0.25*(distance[zp][yp][xx]+distance[zm][ym][xx]-distance[zp][ym][xx]-distance[zm][yp][xx]);
					
					realnum ig2 = 1.0/(gx*gx+gy*gy+gz*gz+1e-99);
					realnum ig = sqrt(ig2);
					gx *= ig; gy *= ig; gz *= ig; // unit here
					unx[zz][yy][xx] = gx;
					uny[zz][yy][xx] = gy;
					unz[zz][yy][xx] = gz;

					/*realnum cr1 = gx * gy * (hxz * hyz - hxy * hzz) + gx * gz * (hxy * hyz - hxz * hyy) + gy * gz * (hxy * hxz - hyz * hxx); cr1 *= ig2;
					realnum cr2 = gx*gx*(hyy*hzz-hyz*hyz)+gy*gy*(hxx*hzz-hxz*hxz)+gz*gz*(hxx*hyy-hxy*hxy); cr2 *= ig2;
					Gau[zz][yy][xx] = 2*cr1+cr2;*/

					realnum cr1 = gx*gx*(hyy+hzz)+gy*gy*(hxx+hzz)+gz*gz*(hxx+hyy); cr1 *= ig;
					realnum cr2 = gx*gy*hxy+gx*gz*hxz+gy*gz*hyz; cr2 *= ig;
					Sum[zz][yy][xx] = cr1-2*cr2;

					thstat[zz][yy][xx] = 1; // calculated

					if (smstat[zz][yp][xx] == 1 && smstat[zz][ym][xx] == 1
					 && smstat[zz][yy][xp] == 1 && smstat[zz][yy][xm] == 1
					 && smstat[zp][yy][xx] == 1 && smstat[zm][yy][xx] == 1
						
						&& smstat[zz][yp][xp] == 1 && smstat[zz][ym][xm] == 1
						&& smstat[zz][yp][xm] == 1 && smstat[zz][ym][xp] == 1
						&& smstat[zp][yy][xp] == 1 && smstat[zm][yy][xm] == 1
						&& smstat[zp][yy][xm] == 1 && smstat[zm][yy][xp] == 1
						&& smstat[zp][yp][xx] == 1 && smstat[zm][ym][xx] == 1
						&& smstat[zm][yp][xx] == 1 && smstat[zp][ym][xx] == 1 

					 && smstat[zz][yy][xx] == 1)
						thstat[zz][yy][xx] = 2;/**/

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
	SVoxImg<SWorkImg<realnum>>& field = m_field[i];

	SVoxImg<SWorkImg<realnum>>& counterfield = m_field[(i + 1) % 2];

	// m_velo
	SVoxImg<SWorkImg<realnum>>& velo = m_velo[i];

	// m_distance[i]
	SVoxImg<SWorkImg<realnum>>& distance = m_distance[i];

	// m_distance[(i+1)&1]
	SVoxImg<SWorkImg<realnum>>& counterd = m_distance[(i + 1) & 1];


	// m_Sumcurvature[i]
	SVoxImg<SWorkImg<realnum>>& Sum = m_Sumcurvature;


	int xs = field.xs, ys = field.ys, zs = field.zs;

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
	velo.Set0(xs, ys, zs);
	{
#pragma omp parallel for collapse(3) shared(maxinfo)
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
					if (meeting_plane_positions[zz][yy][xx]) {
						continue;
					}
					int zp = zz + 1, zm = zz - 1;
					int yp = yy + 1, ym = yy - 1;
					int xp = xx + 1, xm = xx - 1;

					if (distance[zz][yp][xx] < 0 && distance[zz][ym][xx] < 0
						&& distance[zz][yy][xp] < 0 && distance[zz][yy][xm] < 0
						&& distance[zp][yy][xx] < 0 && distance[zm][yy][xx] < 0) continue;


					if (field[zz][yy][xx] >= 0.95) continue;

					if (counterd[zz][yy][xx] >= 0) {
						// meeting_plane.insert(IPoi3<int>(xx, yy, zz));
						if ((distance[zz][yp][xx] >= 0 && !meeting_plane_positions[zz][yp][xx])
							|| (distance[zz][ym][xx] >= 0 && !meeting_plane_positions[zz][ym][xx])
							|| (distance[zz][yy][xp] >= 0 && !meeting_plane_positions[zz][yy][xp])
							|| (distance[zz][yy][xm] >= 0 && !meeting_plane_positions[zz][yy][xm])
							|| (distance[zp][yy][xx] >= 0 && !meeting_plane_positions[zp][yy][xx])
							|| (distance[zm][yy][xx] >= 0 && !meeting_plane_positions[zm][yy][xx]))
							meeting_plane_positions[zz][yy][xx] = 1;
						continue;
					}
					realnum xn = field[zz][yy][xp] - field[zz][yy][xm];
					realnum yn = field[zz][yp][xx] - field[zz][ym][xx];
					realnum zn = field[zp][yy][xx] - field[zm][yy][xx];
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
										if (m_thickstate[pz][py][px] < 0) continue;
										realnum dot = unx[pz][py][px] * ix + uny[pz][py][px] * iy + unz[pz][py][px] * iz;
										dot *= -1;
										if (dot < 0.82) continue;
										++iok; if (dot > 0.94) ++iok;
										wAct += dot / r2;
										sumcur += (dot / r2) * Sum[pz][py][px];
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

					realnum eikon = m_data[zz][yy][xx];

					realnum d0 = m_data[zz][yy][xx];
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
							eikon = m_data[zz][yy][xx];
						}
					}
					else
						eikon = m_data[zz][yy][xx];

					eikon *= gradlen; // normalize
					if (eikon < 1e-11) eikon = 1e-11;
					velo[zz][yy][xx] = eikon; // dS = 1

					/*#pragma omp critical(meeting_plane)
					{
						if (eikon > maxv) maxv = eikon;
						if (counterd[zz][yy][xx] >= 0) {
							// meeting_plane.insert(IPoi3<int>(xx, yy, zz));
							meeting_plane_positions[zz][yy][xx] = 1;
						}
					}*/
					if (eikon > maxinfo[id].val) maxinfo[id].val = eikon;
					if (counterd[zz][yy][xx] >= 0) {
						// meeting_plane.insert(IPoi3<int>(xx, yy, zz));
						meeting_plane_positions[zz][yy][xx] = 1;
					}

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
	SVoxImg<SWorkImg<realnum>>& field = m_field[i];
	SVoxImg<SWorkImg<realnum>>& velo = m_velo[i];

	int xs = field.xs, ys = field.ys, zs = field.zs;
#pragma omp parallel for
	for (int zz = 1; zz < zs - 1; ++zz) {
		for (int yy = 0 + 1; yy < ys - 1; ++yy) {
			for (int xx = 0 + 1; xx < xs - 1; ++xx) {
				field[zz][yy][xx] += velo[zz][yy][xx] * maxv;// ;
				if (field[zz][yy][xx] > 1.25) field[zz][yy][xx] = 1.25;
			}
		}
	}
}
void CPhaseContainer::UpdateDistance(int i, realnum current_distance) {
	SVoxImg<SWorkImg<realnum>>& field = m_field[i];
	// m_phasefield.m_distance[i]
	SVoxImg<SWorkImg<realnum>>& distance = m_distance[i];

	// m_phasefield.m_distance[(i+1)&1]
	SVoxImg<SWorkImg<realnum>>& counterd = m_distance[(i + 1) & 1];
	int xs = field.xs, ys = field.ys, zs = field.zs;
#pragma omp parallel for
	for (int zz = 1; zz < zs - 1; ++zz) {
		for (int yy = 0 + 1; yy < ys - 1; ++yy) {
			for (int xx = 0 + 1; xx < xs - 1; ++xx) {
				if (field[zz][yy][xx] > 0 && distance[zz][yy][xx] < -0.5f) {
					distance[zz][yy][xx] = current_distance;
				}
				if (distance[zz][yy][xx] < -0.5f && counterd[zz][yy][xx] < -0.5f) {
					m_bdone = false;
				}
			}
		}
	}
}
void CPhaseContainer::Iterate(bool use_correction)
{
	SVoxImg<SWorkImg<realnum>>& field = m_field[0];
	int xs = field.xs, ys = field.ys, zs = field.zs;
	// m_phasefield.m_field[i]
	realnum maxv0 = UpdateVelo(0, use_correction);
	realnum maxv1 = UpdateVelo(1, use_correction);
	realnum maxv = max(maxv0, maxv1);
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

// only for xslice
void CPlanePhaseField::GetDistancemean(SVoxImg<SWorkImg<realnum>>& distance, int xslice)
{
	//TODO: instead of searching for the selected x slice, locate the meeting surface of the two processes, and use that
	int xs = distance.xs, ys = distance.ys, zs = distance.zs;
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
				m_gx2D[zz][yy] = distance[zz][yp][xx - 1] - distance[zz][ym][xx - 1];
				m_gx2D[zz][yy] += distance[zz][yp][xx + 1] - distance[zz][ym][xx + 1];
				m_gy2D[zz][yy] = distance[zp][yy][xx - 1] - distance[zm][yy][xx - 1];
				m_gy2D[zz][yy] += distance[zp][yy][xx + 1] - distance[zm][yy][xx + 1];
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
