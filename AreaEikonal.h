#pragma once

#include "commontype.h"
#include "ImageOp.h"

#include <vector>
#include <sitkImage.h>

namespace sitk = itk::simple;
using namespace std;

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
	SWorkImg<realnum> m_distance2D[2];
	// get the distance gradients along the selected x slice from the two neighboring slices
	void GetDistancemean(SVoxImg<SWorkImg<realnum>> &distance1, SVoxImg<SWorkImg<realnum>>& distance2, int xslice);
	void GetDistancemean(sitk::Image& distance, int xslice);
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
	//SVoxImg<SWorkImg<int>> m_smoothstate;
	sitk::Image m_smoothstate;
	//SVoxImg<SWorkImg<int>> m_thickstate;
	sitk::Image m_thickstate;
	//SVoxImg<SWorkImg<realnum>> m_Sumcurvature;
	sitk::Image m_Sumcurvature;
	//SVoxImg<SWorkImg<realnum>> m_Hessian[6];

	//SVoxImg<SWorkImg<realnum>> m_distance[2];
	sitk::Image m_distance[2];
	//SVoxImg<SWorkImg<realnum>> m_combined_distance;
	sitk::Image m_combined_distance;
	// for every pixel stores which flow reached it first
	//SVoxImg<SWorkImg<int>> m_flow_idx;
	sitk::Image m_flow_idx;

	double m_plane_slice;
	// expansion

	// gradients of the distance maps
	//SVoxImg<SWorkImg<realnum>> unx;
	sitk::Image unx;
	//SVoxImg<SWorkImg<realnum>> uny;
	sitk::Image uny;
	//SVoxImg<SWorkImg<realnum>> unz;
	sitk::Image unz;

	// Working field
	//SVoxImg<SWorkImg<realnum>> m_field[2];
	sitk::Image m_field[2];
	//SVoxImg<SWorkImg<realnum>> m_gradlen;
	//SVoxImg<SWorkImg<realnum>> m_n[3];
	//SVoxImg<SWorkImg<realnum>> m_aux; // temp
	sitk::Image m_aux;
	//SVoxImg<SWorkImg<realnum>> m_smoothaux; // temp
	sitk::Image m_smoothaux;
	//SVoxImg<SWorkImg<realnum>> m_smoothaux2; // temp
	sitk::Image m_smoothaux2;
	//SVoxImg<SWorkImg<realnum>> m_smoothdist; // temp
	sitk::Image m_smoothdist;
	//SVoxImg<SWorkImg<realnum>> m_velo[2];
	sitk::Image m_velo[2];

	sitk::Image m_sample_image;
	// Image data
	//SVoxImg<SWorkImg<realnum>> m_data;
	sitk::Image m_data;

	//SVoxImg<SWorkImg<int>> meeting_plane_positions;
	sitk::Image meeting_plane_positions;
	bool m_bdone;
	realnum m_currentdistance;

	void Initialize(SVoxImg<SWorkImg<realnum>>& data, CVec3& start_point, CVec3& end_point);
	void Initialize(CPhaseContainer& phasefield, vector<double>& rotation_matrix, bool inverse=false);
	void CalculateAlignedCombinedDistance(double p1_x, double p2_x);
	void FindMeetPoints();
	realnum UpdateVelo(int i, bool use_correction);
	void UpdateField(int i, realnum maxv);
	void UpdateDistance(int i, realnum current_distance);
	// Calculate fundamental quantities, like distance gradient, sum curvature and thick state(?)
	void CalculateFundQuant(int i = 0, int test = 0);
	// Creates a smoother version of the distance map in smoothdist
	void SmoothMap(sitk::Image& src1, sitk::Image& src2, sitk::Image& out);
	// Update velocity, then phase field based on velocity, and distance map, where phase field passes threshold value
	void Iterate(bool use_correction);
};

class CCurvEikonal
{
public:
	CCurvEikonal(void);
	~CCurvEikonal(void);
	void ExtractMeetingPlane();
	CPlanePhaseField m_inicountourCalculator;
	IPoi3<double> plane_center;
	IPoi3<double> plane_normal;
	double plane_offset = 1e11;
	std::vector<double> rotation_matrix;
	// phase field stuff
	CPhaseContainer m_phasefield;
	CPhaseContainer m_rotated_phasefield;


	std::vector<CVec3> m_boundcontour;

	// not used
	std::unordered_set<IPoi3<double>, IPoi3Hash<double>> meeting_plane;
};
