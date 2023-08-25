#pragma once
#include <SimpleITK.h>
#include <sitkImageOperators.h>
#include "IndexedPriorityQueue.h"
#include "AreaEikonal.h"
#include "Utils.h"
#include <vector>
#include <math.h>

//#define HIGH_ACCURACY
#define INF 2*m_max_dist//0x1.f7f7f7f7f7f7fp+1016

namespace sitk = itk::simple;
class FMSigned
{
	static FMSigned& getInstance(int id){
		if (id == 0) {
			static FMSigned instance0;
			return instance0;
		}
		static FMSigned instance1;
		return instance1;
	}
	static sitk::Image& getOnes(std::vector<uint32_t>& size) {
		static sitk::Image ones;
		if (ones.GetSize() != size)
			ones = sitk::Image(size, sitk::sitkFloat64);
		return ones;
	}
	FMSigned(FMSigned const&) = delete;
	void operator=(FMSigned const&) = delete;
	FMSigned();
	sitk::Image m_distance_map;
	sitk::Image temp;
	IndexedPriorityQueue<POINT3D_MAP(double), std::greater<double>> m_narrow_band;
	std::vector<std::pair<POINT3D, double>> m_narrow_band_v;
	sitk::Image m_frozen;
	const double* m_input_buffer;
	double* m_distance_buffer;
	uint8_t* m_frozen_buffer;
	double m_max_dist;
	double mVelocityMap;
	int xs, ys, zs;
	bool finished();
	void initialize(sitk::Image& input_map, double max_distance, double velo);
	template<class RandIt>
	void initialize(const sitk::Image& distance_map, const sitk::Image& frozen, const RandIt& bound_points, double max_distance = DBL_MAX, double velo = 1., bool fix_heap = true) {
		_PROFILING;
		m_max_dist = max_distance;
		mVelocityMap = velo;
		m_distance_map = sitk::Image(distance_map);
		m_frozen = sitk::Image(frozen);
		vector<uint32_t> input_size = m_distance_map.GetSize();
		xs = input_size[0];
		ys = input_size[1];
		zs = input_size[2];
		m_distance_buffer = m_distance_map.GetBufferAsDouble();
		m_frozen_buffer = m_frozen.GetBufferAsUInt8();
		vector<double> values;
		values.reserve(bound_points.size());
		//OMP_PARALLEL_FOR
		for (auto it = bound_points.begin(); it != bound_points.end(); it++) {
			auto [xx, yy, zz] = representation_to_point<int>(*it);
			double val = BUF_IDX3D(m_distance_buffer, xs, ys, zs, xx, yy, zz);
			BUF_IDX3D(m_frozen_buffer, xs, ys, zs, xx, yy, zz) = 1;
			values.push_back(val);
		}
		m_narrow_band = IndexedPriorityQueue<POINT3D_MAP(double), std::greater<double>>(bound_points.begin(), bound_points.end(), values.begin(), fix_heap);
	}
	template<class RandIt>
	void initialize(sitk::Image& input_map, const RandIt& bound_points, double max_distance = DBL_MAX, double velo = 1., bool fix_heap=true) {
		_PROFILING;
		vector<uint32_t> input_size = input_map.GetSize();
		xs = input_size[0];
		ys = input_size[1];
		zs = input_size[2];
		double* input_buffer = input_map.GetBufferAsDouble();
		m_max_dist = max_distance;
		mVelocityMap = velo;
		vector<uint32_t> init_size = m_distance_map.GetSize();
		if (init_size != input_size) {
			m_distance_map = sitk::Image(input_size, sitk::sitkFloat64) + INF;
			m_frozen = sitk::Image(input_size, sitk::sitkUInt8);
		}
		else {
			m_distance_map = -m_distance_map + INF;
			m_frozen -= m_frozen;
		}

		m_distance_buffer = m_distance_map.GetBufferAsDouble();
		m_frozen_buffer = m_frozen.GetBufferAsUInt8();
		vector<double> values;
		values.reserve(bound_points.size());
//OMP_PARALLEL_FOR
		for (auto it = bound_points.begin(); it != bound_points.end(); it++) {
			auto [xx, yy, zz] = representation_to_point<int>(*it);
			double val = -BUF_IDX3D(input_buffer, xs, ys, zs, xx, yy, zz);
			BUF_IDX3D(m_distance_buffer, xs, ys, zs, xx, yy, zz) = val;
			BUF_IDX3D(m_frozen_buffer, xs, ys, zs, xx, yy, zz) = 1;
			values.push_back(val);
		}
		m_narrow_band = IndexedPriorityQueue<POINT3D_MAP(double), std::greater<double>>(bound_points.begin(), bound_points.end(), values.begin(), fix_heap);
	}

	template<>
	void FMSigned::initialize< std::vector<POINT3D>>(sitk::Image& input_map, const std::vector<POINT3D>& bound_points, double max_distance, double velo, bool fix_heap) {
		_PROFILING;
		std::vector<uint32_t> input_size = input_map.GetSize();
		xs = input_size[0];
		ys = input_size[1];
		zs = input_size[2];
		double* input_buffer = input_map.GetBufferAsDouble();
		m_max_dist = max_distance;
		mVelocityMap = velo;
		std::vector<uint32_t> init_size = m_distance_map.GetSize();
		if (init_size != input_size) {
			m_distance_map = sitk::Image(input_size, sitk::sitkFloat64) + INF;
			m_frozen = sitk::Image(input_size, sitk::sitkUInt8);
			m_distance_buffer = m_distance_map.GetBufferAsDouble();
			m_frozen_buffer = m_frozen.GetBufferAsUInt8();
		}
		else {
			std::fill(m_distance_buffer, m_distance_buffer + xs * ys * zs, INF);
			memset(m_frozen_buffer, 0, xs * ys * zs);
			//m_distance_map += -m_distance_map + INF;
			//m_frozen -= m_frozen;
		}


		std::vector<double> values;
		
		_SUB_PROFILE(initmaps);
		values.reserve(bound_points.size());
OMP_PARALLEL_FOR
		for (int i = 0; i < bound_points.size(); i++) {
			auto [xx, yy, zz] = representation_to_point<int>(bound_points[i]);
			double val = -BUF_IDX3D(input_buffer, xs, ys, zs, xx, yy, zz);
			BUF_IDX3D(m_distance_buffer, xs, ys, zs, xx, yy, zz) = val;
			BUF_IDX3D(m_frozen_buffer, xs, ys, zs, xx, yy, zz) = 1;
			values[i] = val;
		}
		_END_SUBPROFILE(initmaps);
		_SUB_PROFILE(fill_narrowband);
		m_narrow_band = IndexedPriorityQueue<POINT3D_MAP(double), std::greater<double>>(bound_points.begin(), bound_points.end(), values.begin(), fix_heap);
		_END_SUBPROFILE(fill_narrowband);
	}

	void FMSigned::initialize_for_neighbors(const sitk::Image& input_map, const std::vector<POINT3D>& bound_points, double max_distance = DBL_MAX, double velo=1.) {
		_PROFILING;
		m_input_buffer = input_map.GetBufferAsDouble();
		std::vector<uint32_t> input_size = input_map.GetSize();
		xs = input_size[0];
		ys = input_size[1];
		zs = input_size[2];
		m_max_dist = max_distance;
		mVelocityMap = velo;
		std::vector<uint32_t> init_size = m_distance_map.GetSize();
		if (init_size != input_size) {
			m_distance_map = sitk::Image(input_size, sitk::sitkFloat64) + INF;
			m_frozen = sitk::Image(input_size, sitk::sitkUInt8);
			m_distance_buffer = m_distance_map.GetBufferAsDouble();
			m_frozen_buffer = m_frozen.GetBufferAsUInt8();
		}
		else {
			std::fill(m_distance_buffer, m_distance_buffer + xs * ys * zs, INF);
			memset(m_frozen_buffer, 0, xs * ys * zs);
			//m_distance_map += -m_distance_map + INF;
			//m_frozen -= m_frozen;
		}


		m_narrow_band_v.clear();

		_SUB_PROFILE(initmaps);
		m_narrow_band_v.resize(bound_points.size());
OMP_PARALLEL_FOR
		for (int i = 0; i < bound_points.size(); i++) {
			auto [xx, yy, zz] = representation_to_point<int>(bound_points[i]);
			double val = -BUF_IDX3D(m_input_buffer, xs, ys, zs, xx, yy, zz);
			BUF_IDX3D(m_distance_buffer, xs, ys, zs, xx, yy, zz) = val;
			BUF_IDX3D(m_frozen_buffer, xs, ys, zs, xx, yy, zz) = 1;
			m_narrow_band_v[i] = std::make_pair(bound_points[i],val);
		}
		_END_SUBPROFILE(initmaps);
	}
	bool compute_distance(std::vector<int> p, double& out);
	void iterate();
	void calculateForNeighbors();
	void smooth_distances();
public:
	static void build(int id, sitk::Image& input_map, sitk::Image& output, double max_distance=DBL_MAX, double velo=1.);
	template<class RandIt>
	static void build(int id, sitk::Image& input_map, sitk::Image& output, RandIt& bound_points, double max_distance = DBL_MAX, double velo = 1.) {
		_PROFILING;
			FMSigned& signedDistMap = getInstance(id);
		signedDistMap.initialize(input_map, bound_points, max_distance, velo);

		while (!signedDistMap.finished()) {
			signedDistMap.iterate();
		}

		sitk::Image&& signs = sitk::Cast(sitk::GreaterEqual(input_map, 0), signedDistMap.m_distance_map.GetPixelID()) * (-2) + 1;

		output = signedDistMap.m_distance_map * signs;
	}

	template<class RandIt>
	static void build(int id, sitk::Image& distance_map, sitk::Image& frozen, sitk::Image& output, RandIt& bound_points, double max_distance = DBL_MAX, double velo = 1.) {
		_PROFILING;
		FMSigned& signedDistMap = getInstance(id);
		signedDistMap.initialize(distance_map, frozen, bound_points, max_distance, velo);

		while (!signedDistMap.finished()) {
			signedDistMap.iterate();
		}


		output = signedDistMap.m_distance_map;
	}

	template <class RandIt>
	static sitk::Image& build_for_neighbors(int id, sitk::Image& input_map, RandIt& bound_points, double max_distance = DBL_MAX, double velo = 1., bool smooth=false) {
		_PROFILING;
		FMSigned& signedDistMap = getInstance(id);
		signedDistMap.initialize_for_neighbors(input_map, bound_points, max_distance, velo);
#ifdef DEBUG_CURVATURE
		extern int current_iteration;
		if (current_iteration % 50 == 0) {
			save_image("Y:/BIOMAG/shortest path/curv_debug/dbg_" + std::to_string(id) + "_" + std::to_string(current_iteration) + "_init_sdist.tif", signedDistMap.m_distance_map);
		}
#endif
		//save_image("Y:/BIOMAG/shortest path/init_dist.tif", signedDistMap.m_distance_map);
		signedDistMap.calculateForNeighbors();
		if (smooth) {
			signedDistMap.smooth_distances();
		}
		//save_image("Y:/BIOMAG/shortest path/out_dist.tif", signedDistMap.m_distance_map);
		//sitk::Image&& signs = sitk::Cast(sitk::GreaterEqual(input_map, 0), signedDistMap.m_distance_map.GetPixelID()) * (-2) + 1;

		return signedDistMap.m_distance_map;
	}
};

bool solve_quadratic(double coeff[3], double& output);