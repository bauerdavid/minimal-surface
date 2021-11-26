#include "stdafx.h"
#include "Utils.h"
#include <utility>
#include "AreaEikonal.h"
#include <unordered_set>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <sitkImageFileWriter.h>
#include <sitkAdditionalProcedures.h>

using namespace std;
namespace sitk = itk::simple;

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


void rotate(const std::vector<double> rotation_matrix, const std::vector<double> value, std::vector<double>& dest) {
	dest.clear();
	for (int i = 0; i < 3; i++) {
		double val(0.0);
		for (int j = 0; j < 3; j++) {
			val += rotation_matrix[i * 3 + j] * value[j];
		}
		dest.push_back(val);
	}
}




void save_image(string filename, sitk::Image img) {
	sitk::ImageFileWriter writer;
	writer.SetFileName(filename);
	writer.Execute(img);
}

void find_rotated_size(vector<unsigned int>& original_size, vector<double>& rotation_matrix, vector<unsigned int>& rotated_size) {
	rotated_size.clear();
	int xs(original_size[0]), ys(original_size[1]), zs(original_size[2]);
	vector<double> upper_bound = { -DBL_MAX, -DBL_MAX, -DBL_MAX };
	vector<double> lower_bound = { DBL_MAX, DBL_MAX, DBL_MAX };
	for (double z = 0; z < zs + 1; z += zs) {
		for (double y = 0; y < ys + 1; y += ys) {
			for (double x = 0; x < xs + 1; x += xs) {
				vector<double> point = { x, y, z };
				vector<double> rotated;
				rotate(rotation_matrix, point, rotated);
				std::transform(rotated.begin(), rotated.end(), upper_bound.begin(), upper_bound.begin(), [](double a, double b) { return (a > b) ? a : b; });
				std::transform(rotated.begin(), rotated.end(), lower_bound.begin(), lower_bound.begin(), [](double a, double b) { return (a < b) ? a : b; });
			}
		}
	}
	std::transform(upper_bound.begin(), upper_bound.end(), lower_bound.begin(), std::back_inserter(rotated_size), [](double l, double r) {return ceil(l - r); });
}
