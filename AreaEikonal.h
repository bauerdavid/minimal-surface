#pragma once

#include "commontype.h"
#include "ImageOp.h"

#include <vector>

#define MAXMINPATH 222 

#define XS_ 220 
#define YS_ 155 
#define ZS_ 110 
#ifndef PI
#define PI 3.1415926536
#endif PI

struct IPoi {
	int x;
	int y;
	IPoi():x(0),y(0) {}
	IPoi(int xx, int yy):x(xx),y(yy) {}
};
struct IPoi3 {
	int x;
	int y;
	int z;
	IPoi3():x(0),y(0),z(0) {}
	IPoi3(int xx, int yy, int zz):x(xx),y(yy),z(zz) {}
};

struct SCurvatureTest {
	int x;
	int y;
	int z;
	realnum n[3];
	realnum gradlen;
	realnum sum;
	realnum gau;
};

class CPlanePhaseField
{
public:
	CPlanePhaseField():m_bInited(false) {}
	bool m_bInited;
	int m_npos;
	int m_nall;
	realnum m_relpos;
	int m_nect;
	int m_nfct;
	SWorkImg<realnum> m_gx2D;
	SWorkImg<realnum> m_gy2D;
	SWorkImg<realnum> m_field2D;
	SWorkImg<realnum> m_velo2D;

	void GetDistancemean(SVoxImg<SWorkImg<realnum>> &distance, SVoxImg<SWorkImg<realnum>>& counterdistance);
	void Iterate();
	std::unordered_set<unsigned long> m_bound;
	std::unordered_set<unsigned long> &RetrieveBound();
};

class CPhaseContainer
{
public:
	// Resolve path
	SVoxImg<SWorkImg<int>> m_smoothstate[2];
	SVoxImg<SWorkImg<int>> m_thickstate[2];
	SVoxImg<SWorkImg<realnum>> m_thicknes[2];
	SVoxImg<SWorkImg<realnum>> m_Gaucurvature[2];
	SVoxImg<SWorkImg<realnum>> m_Sumcurvature[2];
	//SVoxImg<SWorkImg<realnum>> m_Hessian[6];

	SVoxImg<SWorkImg<realnum>> m_distance[2];

	// expansion
	SVoxImg<SWorkImg<realnum>> m_distgrad[2][3];
	SVoxImg<SWorkImg<realnum>> m_expdist[2];

	// Working field
	SVoxImg<SWorkImg<realnum>> m_field[2];
	//SVoxImg<SWorkImg<realnum>> m_gradlen;
	//SVoxImg<SWorkImg<realnum>> m_n[3];
	SVoxImg<SWorkImg<realnum>> m_aux;
	SVoxImg<SWorkImg<realnum>> m_smoothaux[2];
	SVoxImg<SWorkImg<realnum>> m_smoothaux2[2];
	SVoxImg<SWorkImg<realnum>> m_smoothdist[2];

	SVoxImg<SWorkImg<realnum>> m_velo;
	// Data driver
	SVoxImg<SWorkImg<realnum>> m_data;
};

class CCurvEikonal
{
public:
	CCurvEikonal(void);
	~CCurvEikonal(void);

	CImageOp m_imageOp;
	CPlanePhaseField m_inicountourCalculator;

	// phase field stuff
	CPhaseContainer m_phasefield;
	void PhaseInit(IPoi reginit, IPoi arrival, int initz = 0, int arravz = 0, int xSection = -1);
	void RegularizePhaseField(SVoxImg<SWorkImg<realnum>> &field, SVoxImg<SWorkImg<realnum>> &velo);
	void SmoothDistanceMap(int i = 0);
	void SmoothMap(SVoxImg<SWorkImg<realnum>> &map, int i = 0);
	void CalculateFundQuant(int i = 0, int test = 0);
	void Iterate(int i = 0);
	realnum m_currentdistance[2];
	CVec3 m_reference[2];
	void ResolvePath(realnum x, realnum y, realnum z, bool bClear = true, int i = 0, int j = 0);
	void ResolvePath(int i = 0);
	std::vector<CVec3> m_minpath[2][MAXMINPATH];

	/**/
	int m_resolvready;
	int m_inittype;
	bool m_bdone;
	IPoi3 m_distanceto;

	int m_j[2];
	std::vector<CVec3> m_boundcontour;
	/**/

	// image general
	int m_expfac;
	int m_valid;
	bool m_grays;
	bool m_color;

	SDisImg m_img;
	SWorkImg<realnum> m_work;
	SWorkImg<realnum> m_workr;
	SWorkImg<realnum> m_workg;
	SWorkImg<realnum> m_workb;

	SWorkImg<realnum> m_intex;
	SWorkImg<realnum> m_intey;
	SWorkImg<realnum> m_divcheck;

	// test
	std::vector<SCurvatureTest> m_ctest;
	// not used

};


//SmoothMap USES thstat
void inline CCurvEikonal::SmoothMap(SVoxImg<SWorkImg<realnum>>& map, int i)
{
	SVoxImg<SWorkImg<realnum>>& sou = map;//!t ? m_phasefield.m_Gaucurvature[i]: 
		//t == 1 ? m_phasefield.m_Sumcurvature[i]:m_phasefield.m_thicknes[i]; // source/dest
	SVoxImg<SWorkImg<realnum>>& out = sou; // output
	SVoxImg<SWorkImg<realnum>>& aux1 = m_phasefield.m_aux;
	SVoxImg<SWorkImg<realnum>>& aux2 = m_phasefield.m_velo;
	SVoxImg<SWorkImg<int>>& thstat = m_phasefield.m_thickstate[i];

	int xs = sou.xs, ys = sou.ys, zs = sou.zs;

#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				if (thstat[zz][yy][xx] != 0) {
					aux1[zz][yy][xx] = sou[zz][yy][xx];
					continue;
				}
				int xm(xx - 1); if (xm == -1) xm = 1;
				int xp(xx + 1); if (xp == xs) xp = xs - 2;
				realnum s = sou[zz][yy][xx]; int w = 0;
				if (thstat[zz][yy][xp] >= 0) {
					s += sou[zz][yy][xp]; ++w;
				}
				if (thstat[zz][yy][xm] >= 0) {
					s += sou[zz][yy][xm]; ++w;
				}
				if (!w) aux1[zz][yy][xx] = s;
				else {
					s += sou[zz][yy][xx];
					if (w == 1) aux1[zz][yy][xx] = (1.0 / 3.0) * s;
					else aux1[zz][yy][xx] = 0.25 * s;
				}
			}
		}
	}
#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			int ym(yy - 1); if (ym == -1) ym = 1;
			int yp(yy + 1); if (yp == ys) yp = ys - 2;
			for (int xx = 0; xx < xs; ++xx) {
				if (thstat[zz][yy][xx] != 0) {
					aux2[zz][yy][xx] = aux1[zz][yy][xx];
					continue;
				}
				realnum s = aux1[zz][yy][xx]; int w = 0;
				if (thstat[zz][yp][xx] >= 0) {
					s += aux1[zz][yp][xx]; ++w;
				}
				if (thstat[zz][ym][xx] >= 0) {
					s += aux1[zz][ym][xx]; ++w;
				}
				if (!w) aux2[zz][yy][xx] = s;
				else {
					s += aux1[zz][yy][xx];
					if (w == 1) aux2[zz][yy][xx] = (1.0 / 3.0) * s;
					else aux2[zz][yy][xx] = 0.25 * s;
				}
			}
		}
	}
#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		int zm(zz - 1); if (zm == -1) zm = 1;
		int zp(zz + 1); if (zp == zs) zp = zs - 2;
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				if (thstat[zz][yy][xx] != 0) continue;
				realnum s = aux2[zz][yy][xx]; int w = 0;
				if (thstat[zp][yy][xx] >= 0) {
					s += aux2[zp][yy][xx]; ++w;
				}
				if (thstat[zm][yy][xx] >= 0) {
					s += aux2[zm][yy][xx]; ++w;
				}
				if (!w) out[zz][yy][xx] = s;
				else {
					s += aux2[zz][yy][xx];
					if (w == 1) out[zz][yy][xx] = (1.0 / 3.0) * s;
					else out[zz][yy][xx] = 0.25 * s;
				}
			}
		}
	}

}
