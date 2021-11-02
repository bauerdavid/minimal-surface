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
template<typename T>
struct IPoi3 {
	T x;
	T y;
	T z;
	IPoi3():x(T()),y(T()),z(T()) {}
	IPoi3(T xx, T yy, T zz):x(xx),y(yy),z(zz) {}
	bool operator==(const IPoi3<T>& p) const {
		return this->x == p.x && this->y == p.y && this->z == p.z;
	}
};

template<typename T>
struct IPoi3Hash {
public:
	size_t operator()(const IPoi3<T>& p) const {
		return ((((long) p.x) * 73856093) ^ (((long) p.y) * 19349663) ^ (((long)p.z) * 83492791)) % 1000000;
	}
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
	int m_nect;
	int m_nfct;
	SWorkImg<realnum> m_gx2D;
	SWorkImg<realnum> m_gy2D;
	SWorkImg<realnum> m_field2D;
	SWorkImg<realnum> m_velo2D;
	// get the distance gradients along the selected x slice from the two neighboring slices
	void GetDistancemean(SVoxImg<SWorkImg<realnum>> &distance, int xslice);
	void Iterate();
	std::unordered_set<unsigned long> m_bound;
	// collects bounding points into a set
	std::unordered_set<unsigned long> &RetrieveBound();
};

class CPhaseContainer
{
public:
	// Resolve path
	/*
	Smooth distance map state
	-4: neither map nor gradient
	0<: was smoothed along every dimension
	*/
	SVoxImg<SWorkImg<int>> m_smoothstate[2];
	SVoxImg<SWorkImg<int>> m_thickstate[2];
	SVoxImg<SWorkImg<realnum>> m_thicknes[2];
	SVoxImg<SWorkImg<realnum>> m_Sumcurvature[2];
	//SVoxImg<SWorkImg<realnum>> m_Hessian[6];

	SVoxImg<SWorkImg<realnum>> m_distance[2];
	SVoxImg<SWorkImg<realnum>> m_combined_distance;

	// expansion

	// gradients of the distance maps
	SVoxImg<SWorkImg<realnum>> m_distgrad[2][3];
	SVoxImg<SWorkImg<realnum>> m_expdist[2];

	// Working field
	SVoxImg<SWorkImg<realnum>> m_field[2];
	//SVoxImg<SWorkImg<realnum>> m_gradlen;
	//SVoxImg<SWorkImg<realnum>> m_n[3];
	SVoxImg<SWorkImg<realnum>> m_aux; // temp
	SVoxImg<SWorkImg<realnum>> m_smoothaux[2]; // temp
	SVoxImg<SWorkImg<realnum>> m_smoothaux2[2]; // temp
	SVoxImg<SWorkImg<realnum>> m_smoothdist[2]; // temp

	SVoxImg<SWorkImg<realnum>> m_velo[2];
	// Image data
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
	// Initialize phasefield, including the surrounding area around the initial points
	void PhaseInit(IPoi reginit, IPoi arrival, int initz = 0, int arravz = 0, int xSection = -1);
	void RegularizePhaseField(SVoxImg<SWorkImg<realnum>> &field, SVoxImg<SWorkImg<realnum>> &velo);
	// Creates a smoother version of the distance map in smoothdist
	void SmoothDistanceMap(int i = 0);
	// Calculate fundamental quantities, like distance gradient, sum curvature and thick state(?)
	void CalculateFundQuant(int i = 0, int test = 0);
	realnum CCurvEikonal::UpdateVelo(int i);
	void CCurvEikonal::UpdateField(int i, realnum maxv);
	void CCurvEikonal::UpdateDistance(int i);
	// Update velocity, then phase field based on velocity, and distance map, where phase field passes threshold value
	void Iterate();
	realnum m_currentdistance[2];
	CVec3 m_reference[2];
	void ResolvePath(realnum x, realnum y, realnum z, bool bClear = true, int i = 0, int j = 0);
	void ResolvePath(int i = 0);
	std::vector<CVec3> m_minpath[2][MAXMINPATH];

	/**/
	bool m_bdone;
	IPoi3<int> m_distanceto;

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
	std::unordered_set<IPoi3<double>, IPoi3Hash<double>> meeting_plane;
	SVoxImg<SWorkImg<int>> meeting_plane_positions;
};
