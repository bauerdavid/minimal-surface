#pragma once
#include <SimpleITK.h>
#include "IndexedPriorityQueue.h"
#include "AreaEikonal.h"
#include "Utils.h"
#include <vector>

//#define HIGH_ACCURACY

namespace sitk = itk::simple;
const std::vector<int> NEIGH6_OFFSET[] = {
	std::vector<int>({-1, 0, 0}),
	std::vector<int>({1, 0, 0}),
	std::vector<int>({0, -1, 0}),
	std::vector<int>({0, 1, 0}),
	std::vector<int>({0, 0, -1}),
	std::vector<int>({0, 0, 1})
};
class FMSigned
{
	FMSigned();
	sitk::Image m_distance_map;
	IndexedPriorityQueue<POINT3D_MAP(double), std::greater<double>> m_narrow_band;
	sitk::Image m_frozen;
	double* m_distance_buffer;
	uint8_t* m_frozen_buffer;
	double m_max_dist;
	double velo;
	int xs, ys, zs;
	bool finished();
	void initialize(sitk::Image& distance_map, double max_distance, double velo);
	bool compute_distance(std::vector<int> p, double& out);
	void iterate();
public:
	static void build(sitk::Image& distance_map, sitk::Image& output, double max_distance=DBL_MAX, double velo=1);
};

int solve_quadratic(double coeff[3], double solutions[2]);