#include "stdafx.h"
#include "FitPlane.h"
#include <utility>
#include "AreaEikonal.h"
#include <unordered_set>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>


std::pair<IPoi3<double>, IPoi3<double>> best_plane_from_points(const std::unordered_set<IPoi3<double>, IPoi3Hash<double>>& c)
{
	// copy coordinates to  matrix in Eigen format
	size_t num_atoms = c.size();
	Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > coord(3, num_atoms);
	std::unordered_set<IPoi3<double>>::iterator it = c.begin();
	for (size_t i = 0; i < num_atoms; ++i, ++it) {
		coord.col(i)[0] = it->x;
		coord.col(i)[1] = it->y;
		coord.col(i)[2] = it->z;
	}
	// calculate centroid
	IPoi3<double> centroid(coord.row(0).mean(), coord.row(1).mean(), coord.row(2).mean());

	// subtract centroid
	coord.row(0).array() -= centroid.x;
	coord.row(1).array() -= centroid.y;
	coord.row(2).array() -= centroid.z;

	// we only need the left-singular matrix here
	//  http://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
	auto svd = coord.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
	auto U = svd.matrixU().data();
	IPoi3<double> plane_normal(U[6], U[7], U[8]);
	return std::make_pair(centroid, plane_normal);
}