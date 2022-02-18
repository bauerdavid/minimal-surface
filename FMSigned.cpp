#include "FMSigned.h"
#include <sitkImageOperators.h>
//#define USE_SITK_FOR_SDIST
using namespace std;

bool solve_quadratic(double coeff[3], double& output) {
	double a = coeff[2], b = coeff[1], c = coeff[0];
	double discriminant = b * b - 4 * a * c;
	if (discriminant >= 0) {
		if(b <= 0)
			output = (-b + sqrt(discriminant)) / (2 * a);
		else
			output = (-b - sqrt(discriminant)) / (2 * a);
		return true;
	}
	return false;
}

FMSigned::FMSigned() {

}

bool FMSigned::finished() {
	return m_narrow_band.empty();
}

void FMSigned::initialize(sitk::Image& input_map, double max_distance, double velo) {
	_PROFILING;
	double* input_buffer = input_map.GetBufferAsDouble();
	m_max_dist = max_distance;
	m_velo = velo;
	vector<uint32_t> init_size = m_distance_map.GetSize();
	vector<uint32_t> input_size = input_map.GetSize();
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
	xs = input_size[0];
	ys = input_size[1];
	zs = input_size[2];
#pragma omp parallel for
	for (int zyx = 0; zyx < zs * ys * xs; zyx++) {
		int zz = zyx / (ys * xs);
		{
			int yy = (zyx / xs) % ys;
			{
				int xx = zyx % xs;
				bool&& neg = BUF_IDX(input_buffer, xs, ys, zs, xx, yy, zz) < 0;
				bool&& xp_pos = xx + 1 < xs && BUF_IDX(input_buffer, xs, ys, zs, xx + 1, yy, zz) > 0;
				bool&& xp_neg = xx + 1 < xs && BUF_IDX(input_buffer, xs, ys, zs, xx + 1, yy, zz) < 0;
				bool&& yp_pos = yy + 1 < ys && BUF_IDX(input_buffer, xs, ys, zs, xx, yy + 1, zz) > 0;
				bool&& yp_neg = yy + 1 < ys && BUF_IDX(input_buffer, xs, ys, zs, xx, yy + 1, zz) < 0;
				bool&& zp_pos = zz + 1 < zs && BUF_IDX(input_buffer, xs, ys, zs, xx, yy, zz + 1) > 0;
				bool&& zp_neg = zz + 1 < zs && BUF_IDX(input_buffer, xs, ys, zs, xx, yy, zz + 1) < 0;

				bool&& xn_pos = xx - 1 > 0 && BUF_IDX(input_buffer, xs, ys, zs, xx - 1, yy, zz) > 0;
				bool&& xn_neg = xx - 1 > 0 && BUF_IDX(input_buffer, xs, ys, zs, xx - 1, yy, zz) < 0;
				bool&& yn_pos = yy - 1 > 0 && BUF_IDX(input_buffer, xs, ys, zs, xx, yy - 1, zz) > 0;
				bool&& yn_neg = yy - 1 > 0 && BUF_IDX(input_buffer, xs, ys, zs, xx, yy - 1, zz) < 0;
				bool&& zn_pos = zz - 1 > 0 && BUF_IDX(input_buffer, xs, ys, zs, xx, yy, zz - 1) > 0;
				bool&& zn_neg = zz - 1 > 0 && BUF_IDX(input_buffer, xs, ys, zs, xx, yy, zz - 1) < 0;
				if ((xp_pos || xn_pos || yp_pos || yn_pos || zp_pos || zn_pos) && neg
					&& (xp_neg || xn_neg || yp_neg || yn_neg || zp_neg || zn_neg)) {
					BUF_IDX(m_distance_buffer, xs, ys, zs, xx, yy, zz) = 0;
					BUF_IDX(m_frozen_buffer, xs, ys, zs, xx, yy, zz) = 1;
					POINT3D point = point_to_representation(xx, yy, zz);
#pragma omp critical (build_narrow_band)
					m_narrow_band.push(point, 0);
				}

			}
		}
	}
}



bool FMSigned::compute_distance(vector<int> p, double& out) {
	//_PROFILING;
	double coeff[] = { -1. / (m_velo * m_velo), 0, 0 };
	bool neg = false;
	for (int dim = 0; dim < 3; dim++) {
		double val1 = INF;
		double val2 = INF;
		for (int i = 0; i < 2; i++) {
			const vector<int>& offset = NEIGH6_OFFSET[2 * dim + i];
			int xn(p[0] + offset[0]), yn(p[1] + offset[1]), zn(p[2] + offset[2]);
			if (xn < 0 || xn >=xs || yn < 0 || yn >=ys || zn < 0 || zn >= zs ||
				!BUF_IDX(m_frozen_buffer, xs, ys, zs, xn, yn, zn)) continue;
			double _val1 = BUF_IDX(m_distance_buffer, xs, ys, zs, xn, yn, zn);
			if (abs(_val1) < abs(val1)) {
				val1 = _val1;
#ifdef HIGH_ACCURACY
				vector<int> pn2;
				transform(p.begin(), p.end(), offset.begin(), back_inserter(pn2), [](int p, int o) { return p + 2 * o; });

				double _val2 = BUF_IDX(m_distance_buffer, xs, ys, zs, pn2[0], pn2[1], pn2[2]);
				if (BUF_IDX(m_frozen_buffer, xs, ys, zs, pn2[0], pn2[1], pn2[2]) && _val2 <= val2)
					val2 = _val2;
				else
					val2 = INF;
#endif
			}
		}
#ifdef HIGH_ACCURACY
		if (val2 < INF) {
			double tp = (4 * val1 - val2) / 3;
			double a = 9.0 / 4;
			coeff[2] += a;
			coeff[1] -= 2 * a * tp;
			coeff[0] += a * tp * tp;
		} else
#endif
		if (abs(val1) < INF) {
			coeff[2] += 1;
			coeff[1] -= 2 * val1;
			coeff[0] += val1 * val1;
		}
	}
	out = INF;
	bool solved = solve_quadratic(coeff, out);
	return solved;
}

void FMSigned::iterate() {
	_PROFILING;
	POINT3D p;
	double dist;
	tie(p, dist) = m_narrow_band.top();
	auto [x, y, z] = representation_to_point<int>(p);
	if (dist > m_max_dist) {
		m_narrow_band.clear();
		return;
	}
	m_narrow_band.pop();
	BUF_IDX(m_frozen_buffer, xs, ys, zs, x, y, z) = 1;
	BUF_IDX(m_distance_buffer, xs, ys, zs, x, y, z) = dist;
//#pragma omp parallel for
	for (int i = 0; i < 6; i++) {
		const vector<int> &offset = NEIGH6_OFFSET[i];
		int xn = x + offset[0], yn = y + offset[1], zn = z + offset[2];
		if (xn < 0 || xn >= xs || yn < 0 || yn >= ys || zn < 0 || zn >= zs
			|| BUF_IDX(m_frozen_buffer, xs, ys, zs, xn, yn, zn)) continue;
		double dist_n;
		bool solved = compute_distance({ xn, yn, zn }, dist_n);
		if (solved)
//#pragma omp critical (push_to_narrow_band)
			m_narrow_band.push_or_promote(point_to_representation(xn, yn, zn), dist_n);
	}
}

void FMSigned::calculateForNeighbors() {
	_PROFILING;
	auto& v = m_narrow_band_v;
	int n = v.size();
#pragma omp parallel
	{
		vector< std::pair<POINT3D, double>> temp;
#pragma omp for
		for (int i = 0; i < n; i++) {
			POINT3D p;
			double dist;
			tie(p, dist) = v[i];
			auto [x, y, z] = representation_to_point<int>(p);
			//#pragma omp parallel for
			for (int i = 0; i < 6; i++) {
				const vector<int>& offset = NEIGH6_OFFSET[i];
				int xn = x + offset[0], yn = y + offset[1], zn = z + offset[2];
				if (xn < 0 || xn >= xs || yn < 0 || yn >= ys || zn < 0 || zn >= zs
					|| BUF_IDX(m_distance_buffer, xs, ys, zs, xn, yn, zn) < INF) continue;
				double dist_n;
				bool solved = compute_distance({ xn, yn, zn }, dist_n);
				if (solved) {
					BUF_IDX(m_distance_buffer, xs, ys, zs, xn, yn, zn) = dist_n;
					temp.push_back(make_pair(point_to_representation(xn, yn, zn), dist_n));
				}

			}
		}
#pragma omp critical
		m_narrow_band_v.insert(m_narrow_band_v.end(), temp.begin(), temp.end());
	}
	/*sort(m_narrow_band_v.begin() + n, m_narrow_band_v.end());
	m_narrow_band_v.erase(unique(m_narrow_band_v.begin()+n, m_narrow_band_v.end()), m_narrow_band_v.end());
	m_narrow_band = IndexedPriorityQueue<POINT3D_MAP(double), std::greater<double>>(m_narrow_band_v.begin() + n, m_narrow_band_v.end());
	while (!finished()) {
		POINT3D p;
		double dist;
		tie(p, dist) = m_narrow_band.top();
		auto [x, y, z] = representation_to_point<int>(p);
		m_narrow_band.pop();
		BUF_IDX(m_frozen_buffer, xs, ys, zs, x, y, z) = 1;
		BUF_IDX(m_distance_buffer, xs, ys, zs, x, y, z) = dist;
		for (int i = 0; i < 6; i++) {
			const vector<int>& offset = NEIGH6_OFFSET[i];
			int xn = x + offset[0], yn = y + offset[1], zn = z + offset[2];
			POINT3D neighb = point_to_representation(xn, yn, zn);
			if (xn < 0 || xn >= xs || yn < 0 || yn >= ys || zn < 0 || zn >= zs
				|| BUF_IDX(m_frozen_buffer, xs, ys, zs, xn, yn, zn) || !m_narrow_band.contains(neighb)) continue;
			double dist_n;
			bool solved = compute_distance({ xn, yn, zn }, dist_n);
			if (solved && abs(m_narrow_band[neighb]) > abs(dist_n)) {
				m_narrow_band.changeAtKey(neighb, dist_n);
			}
		}
	}*/
}


void FMSigned::build(int id, sitk::Image& input_map, sitk::Image& output, double max_distance, double velo) {
	_PROFILING;
	FMSigned& signedDistMap = getInstance(id);
	signedDistMap.initialize(input_map, max_distance, velo);
	while (!signedDistMap.finished()) {
		signedDistMap.iterate();
	}
	sitk::Image&& signs = sitk::Cast(sitk::GreaterEqual(input_map, 0), signedDistMap.m_distance_map.GetPixelID()) * (-2) + 1;
	output = signedDistMap.m_distance_map * signs;
}

void FMSigned::smooth_distances() {
	if(temp.GetSize() != m_distance_map.GetSize())
		temp = sitk::Image(m_distance_map.GetSize(), m_distance_map.GetPixelID());
	double* temp_buffer = temp.GetBufferAsDouble();
	memcpy(temp_buffer, m_distance_buffer, xs * ys * zs * sizeof(double));
	auto& v = m_narrow_band_v;
#pragma omp parallel for
	for (int i = 0; i < v.size(); i++) {
		POINT3D p;
		double dist, mean_dist, sum_dist = 0;
		tie(p, dist) = v[i];
		int n_neighbors = 0;
		auto [x, y, z] = representation_to_point<int>(p);
		for (int i = 0; i < 6; i++) {
			const vector<int>& offset = NEIGH6_OFFSET[i];
			int xn = x + offset[0], yn = y + offset[1], zn = z + offset[2];
			double dist_n;
			if (xn < 0 || xn >= xs || yn < 0 || yn >= ys || zn < 0 || zn >= zs || (dist_n = BUF_IDX(m_distance_buffer, xs, ys, zs, xn, yn, zn)) >= INF) continue;
			n_neighbors++;
			sum_dist += dist_n;
		}
		mean_dist = dist*0.5+sum_dist/(2*n_neighbors);
		BUF_IDX(temp_buffer, xs, ys, zs, x, y, z) = mean_dist;
	}
	m_distance_map = temp;
	m_distance_buffer = m_distance_map.GetBufferAsDouble();
}