#pragma once
#include <utility>
#include "AreaEikonal.h"
#include <numeric>
#include <algorithm>

using namespace std;

std::pair<IPoi3<double>, IPoi3<double>> best_plane_from_points(const std::unordered_set<IPoi3<double>, IPoi3Hash<double>>& c);

void rotate(const std::vector<double> rotation_matrix, const std::vector<double> value, std::vector<double>& dest);

template <typename T> int sgn(T val) {
	return (int)(T(0) < val) - (int)(val < T(0));
}

template<typename T>
void cross_product(const vector<T> vec1, const vector<T> vec2, vector<T>& out) {
	if (vec1.size() != 3 || vec2.size() != 3) {
		return;
	}
	out.push_back(vec1[1] * vec2[2] - vec1[2] * vec2[1]);
	out.push_back(vec1[0] * vec2[2] - vec1[2] * vec2[0]);
	out.push_back(vec1[0] * vec2[1] - vec1[1] * vec2[0]);
}

template<typename T>
vector<double> rotation_matrix_from_vectors(vector<T>& vec1, vector<T>& vec2) {
	T sum1 = accumulate(vec1.begin(), vec1.end(), 0.0);
	T sum2 = accumulate(vec2.begin(), vec2.end(), 0.0);
	transform(vec1.begin(), vec1.end(), vec1.begin(), [sum1](T& c) { return c / sum1; });
	transform(vec2.begin(), vec2.end(), vec2.begin(), [sum2](T& c) { return c / sum2; });
	vector<double> v;
	cross_product(vec1, vec2, v);
	double c = inner_product(vec1.begin(), vec1.end(), vec2.begin(), 0.0);
	double norm = sqrt(inner_product(vec1.begin(), vec1.end(), vec1.begin(), 0.0));
	vector<double> kmat_plus_eye = vector<double>({
		1, -v[2], v[1],
		v[2], 1, -v[0],
		-v[1], v[0], 1 });
	vector<double> kmat_2 = vector<double>({
		(-v[2] * v[2] - v[1] * v[1]),				 (v[0] * v[1]),				   (v[2] * v[0]),
					   (v[1] * v[0]), (-v[2] * v[2] - v[1] * v[1]),				   (v[2] * v[1]),
					   (v[0] * v[2]),				 (v[1] * v[2]), (-v[0] * v[0] - v[2] * v[2]) });
	transform(kmat_2.begin(), kmat_2.end(), kmat_2.begin(), [c, norm](T& val) {return val * ((1 - c) / (norm * norm)); });
	vector<double> rotation_matrix = vector<double>(9);
	transform(kmat_plus_eye.begin(), kmat_plus_eye.end(), kmat_2.begin(), rotation_matrix.begin(), plus<double>());
	return rotation_matrix;
}
