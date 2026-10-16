#pragma once
#include <SimpleITK.h>
#include <sitkImageOperators.h>
#include "IndexedPriorityQueue.h"
#include "AreaEikonal.h"
#include "Utils.h"
#include <vector>
#include <atomic>
#include <memory>
#include <math.h>
#include <cfloat>

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
	// True when m_distance_map holds m_band_fill_value and m_frozen holds 0 everywhere except at
	// the points in m_narrow_band_v, as build_for_neighbors_both leaves them. Every other way of
	// filling the maps clears it, and then the next one fills the whole volume.
	bool m_band_maps_clean = false;
	double m_band_fill_value = 0;
	// Per voxel: has a thread of the current call taken it to calculate its distance (see
	// relax_neighbors_of)? Zero everywhere when no call is running, like the frozen map, except
	// at the points in m_narrow_band_v.
	std::unique_ptr<std::atomic<uint8_t>[]> m_claimed;
	bool finished();
	void initialize(sitk::Image& input_map, double max_distance, double velo);
	template<class RandIt>
	void initialize(const sitk::Image& distance_map, const sitk::Image& frozen, const RandIt& bound_points, double max_distance = DBL_MAX, double velo = 1., bool fix_heap = true) {
		_PROFILING;
		m_band_maps_clean = false;
		m_max_dist = max_distance;
		mVelocityMap = velo;
		m_distance_map = sitk::Image(distance_map);
		m_frozen = sitk::Image(frozen);
		std::vector<uint32_t> input_size = m_distance_map.GetSize();
		xs = input_size[0];
		ys = input_size[1];
		zs = input_size[2];
		m_distance_buffer = m_distance_map.GetBufferAsDouble();
		m_frozen_buffer = m_frozen.GetBufferAsUInt8();
		std::vector<double> values;
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
		m_band_maps_clean = false;
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
		}
		else {
			m_distance_map = -m_distance_map + INF;
			m_frozen -= m_frozen;
		}

		m_distance_buffer = m_distance_map.GetBufferAsDouble();
		m_frozen_buffer = m_frozen.GetBufferAsUInt8();
		std::vector<double> values;
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

	void initialize(sitk::Image& input_map, const std::vector<POINT3D>& bound_points, double max_distance, double velo, bool fix_heap) {
		_PROFILING;
		m_band_maps_clean = false;
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
		for (int i = 0; i < (int)bound_points.size(); i++) {
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

	// The signed distance of the neighbours of a band, in steps, so that the loops of the two
	// instances can be shared (build_for_neighbors_both): begin_prepare_for_neighbors sizes the
	// maps and, if the last call wrote only the points in m_narrow_band_v, leaves their reset to
	// reset_written_point (it returns true then); finish_prepare_for_neighbors makes the list for
	// the new band; write_band_point fills in its points.
	bool begin_prepare_for_neighbors(const sitk::Image& input_map, double max_distance, double velo) {
		bool reset_pending = false;
		m_input_buffer = input_map.GetBufferAsDouble();
		std::vector<uint32_t> input_size = input_map.GetSize();
		xs = input_size[0];
		ys = input_size[1];
		zs = input_size[2];
		m_max_dist = max_distance;
		mVelocityMap = velo;
		std::vector<uint32_t> init_size = m_distance_map.GetSize();
		_SECTION_BEGIN(fill_timer, "side.signed_distance.fill");
		if (init_size != input_size) {
			m_distance_map = sitk::Image(input_size, sitk::sitkFloat64) + INF;
			m_frozen = sitk::Image(input_size, sitk::sitkUInt8);
			m_distance_buffer = m_distance_map.GetBufferAsDouble();
			m_frozen_buffer = m_frozen.GetBufferAsUInt8();
			clear_claims();
		}
		else if (m_band_maps_clean && m_band_fill_value == INF) {
			// The last call only wrote the points in m_narrow_band_v (its band and their
			// neighbours); put those back instead of filling the whole volume.
			reset_pending = true;
		}
		else {
			std::fill(m_distance_buffer, m_distance_buffer + xs * ys * zs, INF);
			memset(m_frozen_buffer, 0, xs * ys * zs);
			clear_claims();
			//m_distance_map += -m_distance_map + INF;
			//m_frozen -= m_frozen;
		}
		m_band_fill_value = INF;
		m_band_maps_clean = false;
		_SECTION_END(fill_timer);
		return reset_pending;
	}
	void clear_claims() {
		size_t voxel_count = (size_t)xs * ys * zs;
		m_claimed.reset(new std::atomic<uint8_t>[voxel_count]);
		for (size_t i = 0; i < voxel_count; i++)
			m_claimed[i].store(0, std::memory_order_relaxed);
	}
	void reset_written_point(int i) {
		auto [xx, yy, zz] = representation_to_point<int>(m_narrow_band_v[i].first);
		BUF_IDX3D(m_distance_buffer, xs, ys, zs, xx, yy, zz) = INF;
		BUF_IDX3D(m_frozen_buffer, xs, ys, zs, xx, yy, zz) = 0;
		m_claimed[xx + (size_t)xs * (yy + (size_t)ys * zz)].store(0, std::memory_order_relaxed);
	}
	void finish_prepare_for_neighbors(const std::vector<POINT3D>& bound_points) {
		m_narrow_band_v.clear();
		m_narrow_band_v.resize(bound_points.size());
	}
	void write_band_point(const std::vector<POINT3D>& bound_points, int i) {
		auto [xx, yy, zz] = representation_to_point<int>(bound_points[i]);
		double val = -BUF_IDX3D(m_input_buffer, xs, ys, zs, xx, yy, zz);
		BUF_IDX3D(m_distance_buffer, xs, ys, zs, xx, yy, zz) = val;
		BUF_IDX3D(m_frozen_buffer, xs, ys, zs, xx, yy, zz) = 1;
		m_narrow_band_v[i] = std::make_pair(bound_points[i],val);
	}
	bool compute_distance(int x, int y, int z, double& out);
	void iterate();
	void relax_neighbors_of(int i, std::vector<std::pair<POINT3D, double>>& found);
	static void calculate_for_neighbors_both(FMSigned* f[2], int team_size);
	void smooth_distances();
public:
	// For the instances 0 and 1 together: the signed distance of the band points of each, and of
	// the neighbours of those that have none, calculated in one team of team_size threads, with
	// loops that go over the points of both. out gets their distance maps.
	static void build_for_neighbors_both(const sitk::Image* input_maps[2], const std::vector<POINT3D>* bound_points[2], double max_distance, double velo, int team_size, const sitk::Image* out[2]) {
		_PROFILING;
		FMSigned* f[2] = { &getInstance(0), &getInstance(1) };
		int n_reset[2] = { 0, 0 };
		for (int s = 0; s < 2; s++) {
			if (f[s]->begin_prepare_for_neighbors(*input_maps[s], max_distance, velo))
				n_reset[s] = (int)f[s]->m_narrow_band_v.size();
		}
		// Reset what the last call wrote, then make the lists of the new bands: the lists are
		// needed for the reset until it is done.
		const int total_reset = n_reset[0] + n_reset[1];
OMP_PARALLEL_FOR_NUM_THREADS(team_size)
		for (int k = 0; k < total_reset; k++) {
			int s = k >= n_reset[0] ? 1 : 0;
			f[s]->reset_written_point(s ? k - n_reset[0] : k);
		}
		for (int s = 0; s < 2; s++)
			f[s]->finish_prepare_for_neighbors(*bound_points[s]);
		const int n0 = (int)bound_points[0]->size(), total = n0 + (int)bound_points[1]->size();
OMP_PARALLEL_FOR_NUM_THREADS(team_size)
		for (int k = 0; k < total; k++) {
			int s = k >= n0 ? 1 : 0;
			f[s]->write_band_point(*bound_points[s], s ? k - n0 : k);
		}
		calculate_for_neighbors_both(f, team_size);
		for (int s = 0; s < 2; s++) {
			f[s]->m_band_maps_clean = true;
			out[s] = &f[s]->m_distance_map;
		}
	}
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
};

bool solve_quadratic(double coeff[3], double& output);