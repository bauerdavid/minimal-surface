#pragma once

#include "commontype.h"
#include "ImageOp.h"

#include <vector>
#include <sitkImage.h>
#include <stack>
#define USE_VECTOR_AS_SET
#define STORE_POINT_AS_INTEGER
#ifdef STORE_POINT_AS_INTEGER
#define POINT3D unsigned long long
#define POINT3D_HASH std::hash<unsigned long long>
#define REMOVABLE_POINT ULLONG_MAX
#else
#define POINT3D IPoi3<int>
#define POINT3D_HASH IPoi3Hash<int>
#define REMOVABLE_POINT POINT3D(-1, -1, -1)
#endif

#define POINT3D_MAP(val_type) POINT3D, val_type, POINT3D_HASH
#define POINT3D_SET POINT3D, POINT3D_HASH

#define CHUNK_SIZE 1//34100
//#define FIND_MEET_POINTS_SCHEDULE schedule(dynamic, CHUNK_SIZE)
//#define SMOOTH_MAP_SCHEDULE schedule(dynamic, CHUNK_SIZE)
//#define CALCULATE_FUND_QUANT_SCHEDULE schedule(dynamic, CHUNK_SIZE)
#define UPDATE_VELO_SCHEDULE schedule(dynamic, CHUNK_SIZE)
//#define UPDATE_DISTANCE_SCHEDULE schedule(dynamic, CHUNK_SIZE)
#define INFO_FILENAME "info_velo_optimized.txt"

#ifndef FIND_MEET_POINTS_SCHEDULE
#define FIND_MEET_POINTS_SCHEDULE
#endif

#ifndef SMOOTH_MAP_SCHEDULE
#define SMOOTH_MAP_SCHEDULE
#endif

#ifndef CALCULATE_FUND_QUANT_SCHEDULE
#define CALCULATE_FUND_QUANT_SCHEDULE
#endif

#ifndef UPDATE_VELO_SCHEDULE
#define UPDATE_VELO_SCHEDULE
#endif

#ifndef UPDATE_DISTANCE_SCHEDULE
#define UPDATE_DISTANCE_SCHEDULE
#endif

#define GET_MACRO(_0, _1, NAME, ...) NAME
#define STRINGIFY(x) ""#x
#define STRINGIFY_EMPTY() "-"
#define TOSTRING(x) STRINGIFY(x)

//#define DEBUG_CURVATURE
//#define DEBUG_STATES
#define ITERATE_ROTATED

namespace sitk = itk::simple;

#define MAXMINPATH 222 
#define INIT_DIST init_dist

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

template<typename T> struct IPoi3;

template<typename T>
inline size_t hash_ipoi3(const IPoi3<T>& p) {
	//return ((((size_t)p.x) * 73856093) ^ (((size_t)p.y) * 19349663) ^ (((size_t)p.z) * 83492791)) % 1000000;
	return (size_t)p.x << (2 * 21) | (size_t)p.y << 21 | (size_t)p.z;
}

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

	bool operator<(const IPoi3<T>& p) const
	{
		
		return (hash_ipoi3<T>(*this) < hash_ipoi3<T>(p));
	}
};


template<typename T>
struct IPoi3Hash {
public:
	
	size_t operator()(const IPoi3<T>& p) const {
		return hash_ipoi3(p);
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
	int m_since_last_update;
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

	CVec3 m_start_point;
	CVec3 m_end_point;

	// Resolve path
	/*
	Smooth distance map state
	-4: neither map nor gradient
	0<: was smoothed along every dimension
	*/
	sitk::Image m_smoothstate;
	sitk::Image m_thickstate;
	sitk::Image m_Sumcurvature;
	sitk::Image m_distance[2];
	sitk::Image m_combined_distance;
	// for every pixel stores which flow reached it first
	sitk::Image m_flow_idx;
	sitk::Image m_distance_connectivity;

	double m_plane_slice;
	// expansion

	// gradients of the distance maps

	sitk::Image m_temp_sdist[2];

	// Working field
	sitk::Image m_field[2];
	sitk::Image m_aux;
	sitk::Image m_smoothaux;
	sitk::Image m_smoothaux2;
	sitk::Image m_smoothdist;
	sitk::Image m_velo[2];

	sitk::Image m_sample_image;
	// Image data
	sitk::Image m_data;
	sitk::Image meeting_plane_positions;
#ifdef USE_VECTOR_AS_SET
	std::vector<POINT3D> active_set[2];
	std::vector<std::pair<POINT3D, double>> m_changed_velo[2];
#else
	unordered_set<POINT3D_SET> active_set[2];
	unordered_map<POINT3D_MAP(double)> m_changed_velo[2];
#endif

	std::unordered_set<POINT3D_SET> meeting_plane;
	IPoi3<double> plane_center = IPoi3<double>(-1, -1, -1);
	IPoi3<double> plane_normal;
	int n_meet_points = 0;
	int m_counter = 0;
	int n_min_meet_points;
	bool n_plane_finalized = false;
	double plane_offset = 1e11;
	std::vector<double> rotation_matrix;
	bool m_bdone;
	realnum m_currentdistance;
	sitk::Image m_curvature[2];
#ifdef DEBUG_CURVATURE
	sitk::Image m_new_update[2];
#endif
	void Initialize(sitk::Image data, CVec3& start_point, CVec3& end_point);
	void Initialize(CPhaseContainer& phasefield, std::vector<double>& rotation_matrix, bool inverse=false);
	void CombineDistance();
	void CalculateAlignedCombinedDistance(double p1_x, double p2_x);
	void InitializeNeighbors();
	void FindMeetPoints();
	void UpdateCurvature(int i);
	realnum GetMinS(int i);
	realnum UpdateVelo(int i, bool use_correction, double S);
	void UpdateField(int i, double maxv);
	void UpdateDistance(int i, realnum current_distance);
	bool IsDone();
	// Calculate fundamental quantities, like distance gradient, sum curvature and thick state(?)
	void CalculateFundQuant(int i = 0, int test = 0);
	// Creates a smoother version of the distance map in smoothdist
	void SmoothMap(sitk::Image& src1, sitk::Image& src2, sitk::Image& out);
	// Update velocity, then phase field based on velocity, and distance map, where phase field passes threshold value
	void Iterate(bool use_correction);
	void ExtractMeetingPlane();
};

class CCurvEikonal
{
public:
	CCurvEikonal(void);
	~CCurvEikonal(void);
	CPlanePhaseField m_inicountourCalculator;

	
	// phase field stuff
	CPhaseContainer m_phasefield;
	CPhaseContainer m_rotated_phasefield;


	std::vector<CVec3> m_boundcontour;

	// not used
};
