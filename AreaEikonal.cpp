#include "StdAfx.h"
#include "AreaEikonal.h"
#include "math.h"

CCurvEikonal::CCurvEikonal(void)
{
	m_valid = 0;
	m_expfac = 9;//+10

	m_j[0] = m_j[1] = 0;
}

CCurvEikonal::~CCurvEikonal(void)
{
}

// Image prep&check

#define GAUTOO 1

// Phasefield stuff

realnum g_w = 2.75;

void CCurvEikonal::PhaseInit(IPoi reginit, IPoi arrival, int zdeparture, int zarrival, int xSection)
{
	if (!m_valid) return;
	if (!m_grays && !m_color) return;
	int spacex, spacey, spacez(ZS_);	
	if (m_grays) { spacex = m_work.xs; spacey = m_work.ys; }
	else { spacex = m_workr.xs; spacey = m_workr.ys; }
	//g_cyc = 0;
	m_phasefield.m_thickstate[0].Set(spacex,spacey,spacez);
	m_phasefield.m_thickstate[1].Set(spacex,spacey,spacez);
	m_phasefield.m_thicknes[0].Set(spacex,spacey,spacez);
	m_phasefield.m_thicknes[1].Set(spacex,spacey,spacez);
	m_phasefield.m_Gaucurvature[0].Set0(spacex,spacey,spacez);
	m_phasefield.m_Gaucurvature[1].Set0(spacex,spacey,spacez);
	m_phasefield.m_Sumcurvature[0].Set0(spacex,spacey,spacez);
	m_phasefield.m_Sumcurvature[1].Set0(spacex,spacey,spacez);

	// expansion
	m_phasefield.m_distgrad[0][0].Set0(spacex,spacey,spacez);
	m_phasefield.m_distgrad[0][1].Set0(spacex,spacey,spacez);
	m_phasefield.m_distgrad[0][2].Set0(spacex,spacey,spacez);
	m_phasefield.m_distgrad[1][0].Set0(spacex,spacey,spacez);
	m_phasefield.m_distgrad[1][1].Set0(spacex,spacey,spacez);
	m_phasefield.m_distgrad[1][2].Set0(spacex,spacey,spacez);
	m_phasefield.m_smoothstate[0].Set(spacex,spacey,spacez);
	m_phasefield.m_smoothstate[1].Set(spacex,spacey,spacez);
	m_phasefield.m_expdist[0].Set(spacex,spacey,spacez);
	m_phasefield.m_expdist[1].Set(spacex,spacey,spacez);
	// expansion


	m_phasefield.m_field[0].Set(spacex,spacey,spacez);
	m_phasefield.m_field[1].Set(spacex,spacey,spacez);
	m_phasefield.m_distance[0].Set(spacex,spacey,spacez);
	m_phasefield.m_distance[1].Set(spacex,spacey,spacez);

	m_phasefield.m_velo.Set0(spacex,spacey,spacez);
	m_phasefield.m_aux.Set0(spacex,spacey,spacez);
	m_phasefield.m_smoothdist[0].Set(spacex,spacey,spacez);
	m_phasefield.m_smoothdist[1].Set(spacex,spacey,spacez);
	m_phasefield.m_smoothaux[0].Set0(spacex,spacey,spacez);
	m_phasefield.m_smoothaux[1].Set0(spacex,spacey,spacez);
	m_phasefield.m_smoothaux2[0].Set0(spacex,spacey,spacez);
	m_phasefield.m_smoothaux2[1].Set0(spacex,spacey,spacez);

	m_phasefield.m_data.Set0(spacex,spacey,spacez);
	SVoxImg<SWorkImg<realnum>> &data = m_phasefield.m_data;
	data = m_imageOp.GetTestInput();

	for (int ii = 0; ii < MAXMINPATH; ++ii) {
		m_minpath[0][ii].clear();
		m_minpath[1][ii].clear();
	}

	SVoxImg<SWorkImg<realnum>> &field0 = m_phasefield.m_field[0];
	SVoxImg<SWorkImg<realnum>> &field1 = m_phasefield.m_field[1];
	SVoxImg<SWorkImg<realnum>> &distance0 = m_phasefield.m_distance[0];
	SVoxImg<SWorkImg<realnum>> &distance1 = m_phasefield.m_distance[1];
	SVoxImg<SWorkImg<realnum>> &thick0 = m_phasefield.m_thicknes[0];
	SVoxImg<SWorkImg<realnum>> &thick1 = m_phasefield.m_thicknes[1];
	SVoxImg<SWorkImg<int>> &thstat0 = m_phasefield.m_thickstate[0];
	SVoxImg<SWorkImg<int>> &thstat1 = m_phasefield.m_thickstate[1];
	#pragma omp parallel for
	for (int zz = 0; zz < field0.zs; ++zz) {
		for (int yy = 0; yy < field0.ys; ++yy) {
			for (int xx = 0; xx < field0.xs; ++xx) {
				field0[zz][yy][xx] = -1.0;
				field1[zz][yy][xx] = -1.0;
				distance0[zz][yy][xx] = -1.0;
				distance1[zz][yy][xx] = -1.0;
				thick0[zz][yy][xx] = -1.0;
				thick1[zz][yy][xx] = -1.0;
				thstat0[zz][yy][xx] = -1;
				thstat1[zz][yy][xx] = -1;
				// expansion
				m_phasefield.m_smoothstate[0][zz][yy][xx] = -4; // neither map nor gradient
				m_phasefield.m_smoothstate[1][zz][yy][xx] = -4;
				m_phasefield.m_expdist[0][zz][yy][xx] = -1.0;
				m_phasefield.m_expdist[1][zz][yy][xx] = -1.0;
				// expansion
			}
		}
	}

	int hs = 11-0;//5 (int)(1.5f*g_w/2);
	m_currentdistance[0] = m_currentdistance[1] = 0;
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
						m_phasefield.m_field[0][zz][yy][xx] = 1.0; // ii->0
						dd = sqrt(dd);
						m_phasefield.m_distance[0][zz][yy][xx] = dd; // ii->0
						m_phasefield.m_thicknes[0][zz][yy][xx] = 1.0; // ii->0
						if (m_currentdistance[ii] < dd) m_currentdistance[ii] = dd;
						
						// curvatures
						//if (dd < 1e-11) dd = 1e-11; dd = 1.0/dd;
						//m_phasefield.m_thickstate[ii][zz][yy][xx] = 2;
						//m_phasefield.m_Gaucurvature[ii][zz][yy][xx] = dd*dd;
						//m_phasefield.m_Sumcurvature[ii][zz][yy][xx] = dd+dd;
					}

				}
			}
		}

	}

	for (int ii = 0; ii < 10-10; ++ii) {
		RegularizePhaseField(m_phasefield.m_field[0],m_phasefield.m_velo);
		RegularizePhaseField(m_phasefield.m_field[1],m_phasefield.m_velo);
	}
	m_phasefield.m_aux[0] = 0;
	m_phasefield.m_aux[1] = 0;

	m_inittype = arrival.x > 0 ? 2:1;

	// new
	if (xSection > 0/*arrival.x < 0 && zarrival > 0*/) { // designated x plane

		for (int zz = 0; zz < spacez; ++zz) {
			for (int yy = 0; yy < spacey; ++yy) {
				//for (int xx = zarrival+1; xx < spacex; ++xx) {
				//	m_phasefield.m_distance[1][zz][yy][xx] = 1e11;
				//}
				m_phasefield.m_distance[1][zz][yy][xSection] = 1e11;
			}
		}
		m_inittype = 3;

	}
	// new

	m_reference[0].x = (realnum)reginit.x;
	m_reference[0].y = (realnum)reginit.y;
	m_reference[0].z = (realnum)zdeparture;
	m_reference[1].x = (realnum)arrival.x;
	m_reference[1].y = (realnum)arrival.y;
	m_reference[1].z = (realnum)zarrival;

	m_distanceto = IPoi3(arrival.x,arrival.y,zarrival);
	m_resolvready = 0; // arrival bumping

	m_bdone = false;
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

int g_ign, g_ian;

void CCurvEikonal::SmoothDistanceMap(int i)
{
	SVoxImg<SWorkImg<realnum>> &sou = m_phasefield.m_distance[i]; // source
	SVoxImg<SWorkImg<realnum>> &aux2 = m_phasefield.m_smoothaux2[i];
	SVoxImg<SWorkImg<realnum>> &aux = m_phasefield.m_smoothaux[i];
	SVoxImg<SWorkImg<realnum>> &out = m_phasefield.m_smoothdist[i]; // output
	
	SVoxImg<SWorkImg<int>> &smstat = m_phasefield.m_smoothstate[i];

	int xs = sou.xs, ys = sou.ys, zs = sou.zs;

	#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				int smst = smstat[zz][yy][xx];
				if (smst > 0) continue;
				smstat[zz][yy][xx] = -2;
				if (sou[zz][yy][xx] < 0) {
					aux[zz][yy][xx] = sou[zz][yy][xx];
					continue;
				}
				int xm(xx-1); if (xm == -1) xm = 1;
				int xp(xx+1); if (xp == xs) xp = xs-2;
				realnum s  = sou[zz][yy][xx]; int w = 0;
				realnum comp = sou[zz][yy][xp];
				if (comp >= 0) { s += comp; ++w; }
				comp = sou[zz][yy][xm];
				if (comp >= 0) { s += comp; ++w; }
				if (!w) aux[zz][yy][xx] = s;
				else {
					s += sou[zz][yy][xx];
					if (w == 1) aux[zz][yy][xx] = (1.0/3.0)*s;
					else {
						aux[zz][yy][xx] = 0.25*s;
						++smstat[zz][yy][xx]; //
					}
				}
			}
		}
	}
	#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			int ym(yy-1); if (ym == -1) ym = 1;
			int yp(yy+1); if (yp == ys) yp = ys-2;
			for (int xx = 0; xx < xs; ++xx) {
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
	#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		int zm(zz-1); if (zm == -1) zm = 1;
		int zp(zz+1); if (zp == zs) zp = zs-2;
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
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

void CCurvEikonal::CalculateFundQuant(int i, int test)
{
	SVoxImg<SWorkImg<int>> &thstat = m_phasefield.m_thickstate[i];
	SVoxImg<SWorkImg<realnum>> &Gau = m_phasefield.m_Gaucurvature[i];
	SVoxImg<SWorkImg<realnum>> &Sum = m_phasefield.m_Sumcurvature[i];
	SVoxImg<SWorkImg<realnum>> &unx = m_phasefield.m_distgrad[i][0];
	SVoxImg<SWorkImg<realnum>> &uny = m_phasefield.m_distgrad[i][1];
	SVoxImg<SWorkImg<realnum>> &unz = m_phasefield.m_distgrad[i][2];

	SmoothDistanceMap(i);	SVoxImg<SWorkImg<realnum>>& distance = m_phasefield.m_smoothdist[i]; // m_phasefield.m_distance[i];
	SVoxImg<SWorkImg<int>>& smstat = m_phasefield.m_smoothstate[i];

	int xs = distance.xs, ys = distance.ys, zs = distance.zs;

	////////////////////////////////////////////////////////////////////////////////
	// curvature update
	////////////////////////////////////////////////////////////////////////////////
	static int nsh; nsh = 0;
	//int ign(0), ian(0);

	#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		int zm(zz-1); if (zm == -1) zm = 1;
		int zp(zz+1); if (zp == zs) zp = zs-2;
		for (int yy = 0; yy < ys; ++yy) {
			int yp(yy+1); if (yp == ys) yp = ys-2;
			int ym(yy-1); if (ym == -1) ym = 1;
			for (int xx = 0; xx < xs; ++xx) {
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

int g_rr[12]={0,0,0,0,0,0,0,0,0,0,0,0};
bool g_modeswitch = false;
realnum g_t[9999];

void CCurvEikonal::Iterate(int i)
{
	static int call(0);

	SVoxImg<SWorkImg<realnum>> &field = m_phasefield.m_field[i];
	SVoxImg<SWorkImg<realnum>> &velo = m_phasefield.m_velo;
	SVoxImg<SWorkImg<realnum>> &distance = m_phasefield.m_distance[i];
	SVoxImg<SWorkImg<realnum>> &counterd = m_phasefield.m_distance[(i+1)&1];
	SVoxImg<SWorkImg<realnum>> &data = m_phasefield.m_data;
	SVoxImg<SWorkImg<realnum>> &Gau = m_phasefield.m_Gaucurvature[i];
	SVoxImg<SWorkImg<realnum>> &Sum = m_phasefield.m_Sumcurvature[i];
	SVoxImg<SWorkImg<realnum>> &thick = m_phasefield.m_thicknes[i];
	SVoxImg<SWorkImg<realnum>> &unx = m_phasefield.m_distgrad[i][0];
	SVoxImg<SWorkImg<realnum>> &uny = m_phasefield.m_distgrad[i][1];
	SVoxImg<SWorkImg<realnum>> &unz = m_phasefield.m_distgrad[i][2];
	SVoxImg<SWorkImg<int>> &thstat = m_phasefield.m_thickstate[i];

	int xs = field.xs, ys = field.ys, zs = field.zs;

	////////////////////////////////////////////////////////////////////
	// fundamental quantities
	////////////////////////////////////////////////////////////////////
	
	CalculateFundQuant(i); // for object-level evolution

	////////////////////////////////////////////////////////////////////
	// evolution logic
	////////////////////////////////////////////////////////////////////

	realnum maxv(0);
	bool modeswitch = g_modeswitch;

	int ites(0);

	{ // Inhomogeneous
	#pragma omp parallel for shared(maxv)
	for (int zz = 1; zz < zs-1; ++zz) {
		int zp = zz+1, zm = zz-1;
		for (int yy = 1; yy < ys-1; ++yy) {
			int yp = yy+1, ym = yy-1;
			for (int xx = 1; xx < xs-1; ++xx) {

				velo[zz][yy][xx] = 0; // zero out
				if (counterd[zz][yy][xx] >= 0) continue;
				int xp = xx+1, xm = xx-1;

				if (distance[zz][yp][xx] >= 0 || distance[zz][ym][xx] >= 0 
				 || distance[zz][yy][xp] >= 0 || distance[zz][yy][xm] >= 0
				 || distance[zp][yy][xx] >= 0 || distance[zm][yy][xx] >= 0) {
				 
					if (field[zz][yy][xx] < 0.95) { 

						realnum xn = field[zz][yy][xp]-field[zz][yy][xm];
						realnum yn = field[zz][yp][xx]-field[zz][ym][xx];
						realnum zn = field[zp][yy][xx]-field[zm][yy][xx];
						realnum gradlen = sqrt(xn*xn+yn*yn+zn*zn+1e-99);
						

						realnum scw(0.0), gcw(0.0);
						{
							//Gau[zz][yy][xx] = 0;
							//Sum[zz][yy][xx] = 0;
							realnum wAct(0.0);
							int iok = -1;
							for (int rr = 1; rr <= 25; ++rr) {
								for (int iz = -rr; iz <= rr; ++iz) {
									realnum z2 = (realnum)iz*iz;
									for (int iy = -rr; iy <= rr; ++iy) {
										realnum y2 = (realnum)iy*iy;
										for (int ix = -rr; ix <= rr; ++ix) {
											if (!ix && !iy && !iz) continue;
											realnum r2 = (realnum)ix*ix;
											r2 += y2+z2;

											if (r2 <= rr*rr) {
												int pz = zz+iz; if (pz < 1) pz = 1; if (pz > zs-2) pz = zs-2;
												int py = yy+iy; if (py < 1) py = 1; if (py > ys-2) py = ys-2;
												int px = xx+ix; if (px < 1) px = 1; if (px > xs-2) px = xs-2;
												if (thstat[pz][py][px] < 0) continue;
												realnum dot = unx[pz][py][px]*ix+uny[pz][py][px]*iy+unz[pz][py][px]*iz;
												dot *= -1;
												if (dot < 0.82) continue;
												++iok; if (dot > 0.94) ++iok;
												realnum w = dot/r2;
												wAct += w;
												scw += w*Sum[pz][py][px]; gcw += w*Gau[pz][py][px];
											}


										}
									}
								}
								if (wAct > 0 && iok > 0) {
									/*Gau[zz][yy][xx] = gcw / wAct;
									Sum[zz][yy][xx] = scw/wAct;
									++g_rr[rr];*/
									scw /= wAct;
									break;
								}

							}

						}
						

						realnum eikon = data[zz][yy][xx];
						realnum sumcur = scw;//Sum[zz][yy][xx];

						realnum d0 = data[zz][yy][xx];
						if (d0 < 1e-33) d0 = 1e-33;

						realnum squar = -1.0;
						realnum discr = 1 + 4 * sumcur * d0 * 1;// *gradlen;

						if (abs(sumcur) > 1e-33) {
							squar = -1.0+sqrt(discr);
							squar /= 2*sumcur;
						}
						else {
							squar = d0;
						}
						
						if (modeswitch) {
							//realnum den = 1.0 + sumcur/1.2; if (den < 1e-6) den = 1e-6;	eikon = data[zz][yy][xx]/den;
							if (squar > 0) {
								eikon = squar;
								//eikon *= exp(-sumcur); // best
							}
							else {
								eikon = data[zz][yy][xx];
							}
						}
						else 
							eikon = data[zz][yy][xx];

						eikon *= gradlen; // normalize
						if (eikon < 1e-11) eikon = 1e-11;
						velo[zz][yy][xx] = eikon; // dS = 1

						#pragma omp critical
						{
						if (eikon > maxv) maxv = eikon;
						g_t[ites] = discr;
						if (++ites >= 9999) ites = 0;
						}

					}
				} 
				//loop body
			}
		}
	}

	}
	
	////////////////////////////////////////////////////////////////////
	// conditioning (speed limitation)
	////////////////////////////////////////////////////////////////////
	if (maxv < 1e-11) maxv = 1e-11;

	realnum mdatspedmax = 0.5 *4; //*2 set max speed
	maxv = mdatspedmax/maxv;
	m_currentdistance[i] += maxv;
	//m_currentdistance[i] += 1;

	#pragma omp parallel for
	for (int zz = 1; zz < zs-1; ++zz) {
		for (int yy = 0+1; yy < ys-1; ++yy) {
			for (int xx = 0+1; xx < xs-1; ++xx) {
				field[zz][yy][xx] += velo[zz][yy][xx]*maxv;// ;
				if (field[zz][yy][xx] > 1.25) field[zz][yy][xx] = 1.25;
			}
		}
	}

	m_bdone = true;

//colltest:

	#pragma omp parallel for
	for (int zz = 1; zz < zs - 1; ++zz) {
		for (int yy = 0+1; yy < ys-1; ++yy) {
			for (int xx = 0+1; xx < xs-1; ++xx) {
				if (field[zz][yy][xx] > 0 && distance[zz][yy][xx] < -0.5f) {
					distance[zz][yy][xx] = m_currentdistance[i];
				}
				if (distance[zz][yy][xx] < -0.5f && counterd[zz][yy][xx] < -0.5f) {
					m_bdone = false;
				}
			}
		}
	}
	////////////////////////////////////////////////////////////////////
	// phase field regularization w/o velocity zeroed out
	////////////////////////////////////////////////////////////////////

	if (m_inittype == 2) { 
		if (field[m_distanceto.z][m_distanceto.y][m_distanceto.x] > 0.9f) {
			++m_resolvready; 
		}
	} 

	++call;
}

// only for xslice
void CPlanePhaseField::GetDistancemean(SVoxImg<SWorkImg<realnum>>& distance, SVoxImg<SWorkImg<realnum>>& counterdistance)
{
	int xs = distance.xs, ys = distance.ys, zs = distance.zs;
	m_nall = m_npos = 0;
	m_bInited = false;
	int xslice(0);
	for (int xx = 1; xx < xs - 1; ++xx) {
		if (counterdistance[zs / 2][ys / 2][xx] > 0) {
			xslice = xx;
			break;
		}
	}
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

			if (!(zz < 5 || zz >= zs - 5 || yy < 5 || yy >= ys - 5))
				m_field2D[zz][yy] = -1.0f;
			else ++m_npos;
			++m_nall;
		}
	}

	m_relpos = ((realnum)m_npos) / ((realnum)m_nall);
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

void CCurvEikonal::ResolvePath(realnum x, realnum y, realnum z, bool bClear, int i, int j)
{
	SVoxImg<SWorkImg<realnum>> &distance = m_phasefield.m_distance[i]; // m_smoothdist 
	int xs = distance.xs, ys = distance.ys, zs = distance.zs;

	if (bClear) m_minpath[i][j].clear();

	int ix((int)x), iy((int)y), iz((int)z);
	if (iz < 2 || iz >= ZS_-2) return;
	if (ix < 2 || ix >= xs-2) return;
	if (iy < 2 || iy >= ys-2) return;

	if (distance[iz][iy][ix] < 0) return;

	CVec3 path(ix,iy,iz);
	m_minpath[i][j].push_back(path);
	
	for (int ii = 0; ii < 11111; ++ii) {
		ix = (int)path.x; iy = (int)path.y; iz = (int)path.z;

		if (iz < 2) iz = 2; if (iz >= ZS_-2) iz = ZS_-3;
		if (ix < 2) ix = 2; if(ix >= xs-2) ix = xs-3;
		if (iy < 2) iy = 2; if(iy >= ys-2) iy = ys-3;

		{	
			realnum mmin(0);
			int pz(0), py(0), px(0);
			for (int zo = -1; zo <= 1; ++zo) {
				int zz = iz+zo;
				for (int yo = -1; yo <= 1; ++yo) {
					int yy = iy+yo;
					for (int xo = -1; xo <= 1; ++xo) {
						int xx = ix+xo;
						if (!zo && !yo && !xo) continue;


						CVec3 dir(xo,yo,zo);

						realnum d1 = distance[zz][yy][xx], d0 = distance[iz][iy][ix];
						if (d1-d0 < mmin) {
							mmin = d1-d0;
							pz = zo; py = yo; px = xo;
						}


					}
				}
			}

			path.x = ix+px; path.y = iy+py; path.z = iz+pz;

		}

		m_minpath[i][j].push_back(path);
		


		CVec3 dir = path; dir -= m_reference[i];
		if (dir.x*dir.x+dir.y*dir.y+dir.z*dir.z < 1.75f) {
			m_minpath[i][j].push_back(m_reference[i]);
			break;
		}

	}

}

void CCurvEikonal::ResolvePath(int i)
{
	m_inittype = 0; 
	m_resolvready = 0;
	ResolvePath(m_distanceto.x,m_distanceto.y,m_distanceto.z,true,i,m_j[i]);
	if (m_minpath[i][m_j[i]].size()) ++m_j[i];
}




