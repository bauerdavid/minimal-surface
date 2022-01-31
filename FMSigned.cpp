#include "FMSigned.h"
#include <sitkImageOperators.h>
#define INF (2*m_max_dist)
using namespace std;

int solve_quadratic(double coeff[3], double solutions[2]) {
	double a = coeff[2], b = coeff[1], c = coeff[0];
	double discriminant = b * b - 4 * a * c;
	if (discriminant > 0) {
		solutions[0] = (-b + sqrt(discriminant)) / (2 * a);
		solutions[1] = (-b - sqrt(discriminant)) / (2 * a);
		return 2;
	}
	else if (discriminant == 0) {
		solutions[0] = -b / (2 * a);
		return 1;
	}
	return 0;
}
FMSigned::FMSigned() {

}

bool FMSigned::finished() {
	return m_narrow_band.empty();
}

void FMSigned::initialize(sitk::Image& distance_map, sitk::Image& output, double max_distance, double velo) {
	vector<uint32_t> size = distance_map.GetSize();
	sitk::PixelIDValueEnum type = distance_map.GetPixelID();
	double* input_buffer = distance_map.GetBufferAsDouble();
	m_max_dist = max_distance;
	output = sitk::Image(size, sitk::sitkFloat64)+ INF;
	m_distance_map = output;
	m_distance_buffer = m_distance_map.GetBufferAsDouble();
	m_frozen = sitk::Image(size, sitk::sitkUInt8);
	m_frozen_buffer = m_frozen.GetBufferAsUInt8();
	xs = size[0];
	ys = size[1];
	zs = size[2];
	for (int zz = 0; zz < zs; zz++) {
		for (int yy = 0; yy < ys; yy++) {
			for (int xx = 0; xx < xs; xx++) {
				bool&& pos = BUF_IDX(input_buffer, xs, ys, zs, xx, yy, zz) > 0;
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
				if (pos && (xp_pos || xn_pos || yp_pos || yn_pos || zp_pos || zn_pos)
					&& (xp_neg || xn_neg || yp_neg || yn_neg || zp_neg || zn_neg)) {
					BUF_IDX(m_distance_buffer, xs, ys, zs, xx, yy, zz) = 0;
					BUF_IDX(m_frozen_buffer, xs, ys, zs, xx, yy, zz) = 1;
					POINT3D point = point_to_representation(xx, yy, zz);
					m_narrow_band.push(point, 0);
				}

			}
		}
	}
}

bool FMSigned::compute_distance(vector<int> p, double& out) {
	double coeff[] = {-1, 0, 0};
	for (int dim = 0; dim < 3; dim++) {
		double val1 = INF;
		double val2 = INF;
		for (int i = 0; i < 2; i++) {
			const vector<int>& offset = NEIGH6_OFFSET[2 * dim + i];
			vector<int> pn;
			transform(p.begin(), p.end(), offset.begin(), back_inserter(pn), plus<int>());
			if (!BUF_IDX(m_frozen_buffer, xs, ys, zs, pn[0], pn[1], pn[2])) continue;
			double _val1 = BUF_IDX(m_distance_buffer, xs, ys, zs, pn[0], pn[1], pn[2]);
			if (_val1 < val1) {
				val1 = _val1;
#ifdef HIGH_ACCURACY
				vector<int> pn2;
				transform(p.begin(), p.end(), offset.begin(), back_inserter(pn2), [](int p, int o) { return p + 2 * o; });

				double _val2 = BUF_IDX(distance_buffer, xs, ys, zs, pn2[0], pn2[1], pn2[2]);
				if (BUF_IDX(frozen_buffer, xs, ys, zs, pn2[0], pn2[1], pn2[2]) && _val2 <= val2)
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
		if (val1 < INF) {
			coeff[2] += 1;
			coeff[1] -= 2 * val1;
			coeff[0] += val1 * val1;
		}
	}
	out = INF;
	double sol[2] = { INF, INF };
	if (int n_solutions = solve_quadratic(coeff, sol)) {
		if (n_solutions == 2)
			out = sol[0] > sol[1] ? sol[0] : sol[1];
		else
			out = sol[0];
		return true;
	}
	return false;
}

void FMSigned::iterate() {
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
	for (int i = 0; i < 6; i++) {
		vector<int> offset = NEIGH6_OFFSET[i];
		int xn = x + offset[0], yn = y + offset[1], zn = z + offset[2];
		if (BUF_IDX(m_frozen_buffer, xs, ys, zs, xn, yn, zn)) continue;
		double dist_n;
		bool solved = compute_distance({ xn, yn, zn }, dist_n);
		if (solved)
			m_narrow_band.push_or_promote(point_to_representation(xn, yn, zn), dist_n);
		if (m_narrow_band.contains(point_to_representation(36, 84, 45)) && m_narrow_band[point_to_representation(36, 84, 45)] < 0.7) {
			double val = m_narrow_band[point_to_representation(36, 84, 45)];
			continue;
		}
	}
}

void FMSigned::build(sitk::Image& distance_map, sitk::Image& output, double max_distance, double velo) {
	FMSigned signedDistMap;
	signedDistMap.initialize(distance_map, output, max_distance, velo);
	save_image("Y:/BIOMAG/shortest path/sdist_init.tif", signedDistMap.m_distance_map);
	double d = -1;
	while (!signedDistMap.finished()) {
		if (signedDistMap.m_narrow_band.top().second > d) {
			//save_image("Y:/BIOMAG/shortest path/sdist_unsigned_temp_" + std::to_string(d) + ".tif", signedDistMap.m_distance_map);
			d = signedDistMap.m_narrow_band.top().second;
		}
		signedDistMap.iterate();
	}
	save_image("Y:/BIOMAG/shortest path/sdist_unsigned.tif", signedDistMap.m_distance_map);
	sitk::Image&& signs = sitk::Cast(sitk::GreaterEqual(distance_map, 0), signedDistMap.m_distance_map.GetPixelID()) * 2 - 1;
	output = signedDistMap.m_distance_map * signs;
}