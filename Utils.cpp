#include "stdafx.h"
#include "Utils.h"
#include <utility>
#include "AreaEikonal.h"
#include <unordered_set>
#include <vector>
#include <sitkImageFileWriter.h>
#include <sitkAdditionalProcedures.h>
#include <sitkCastImageFilter.h>
#include <sitkImageOperators.h>
#include "SimpleITK.h"

using namespace std;
namespace sitk = itk::simple;



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




void save_image(string filename, sitk::Image& img) {
	sitk::ImageFileWriter writer;
	sitk::Image&& cast_img = sitk::Cast(img, sitk::sitkFloat32);
	writer.SetFileName(filename);
	writer.Execute(cast_img);
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

void neg_to_minus1(sitk::Image& img) {
	img = sitk::Cast(sitk::GreaterEqual(img, 0), img.GetPixelID()) * img - sitk::Cast(sitk::Less(img, 0), img.GetPixelID());
}
