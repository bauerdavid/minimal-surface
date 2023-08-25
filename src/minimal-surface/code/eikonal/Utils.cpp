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

#ifdef PROFILE_FUNCTIONS
profiler::profiler(std::string name, bool push) : name{ name }, start_time{ std::chrono::steady_clock::now() }, ended_profiling(false) {
	if (push)
		profile_manager::push_func(name);
	else
		ended_profiling = true;
}
profiler::~profiler() {
	end_profiling();
}
void profiler::end_profiling() {
	if (!ended_profiling) {
		profile_manager::record_info(name, get_execution_time());
		profile_manager::pop_func();
		std::string called_from = profile_manager::current_func();
		profile_manager::add_func_as_child(called_from, name);
		ended_profiling = true;
	}
}

double profiler::get_execution_time() {
	auto now = std::chrono::steady_clock::now();
	return std::chrono::duration_cast<std::chrono::nanoseconds>(now - start_time).count();
}
#endif

void rotate(const std::vector<double>& rotation_matrix, const std::vector<double> value, std::vector<double>& dest) {
	dest.clear();
	for (int i = 0; i < 3; i++) {
		double val(0.0);
		for (int j = 0; j < 3; j++) {
			val += rotation_matrix[i * 3 + j] * value[j];
		}
		dest.push_back(val);
	}
}




void save_image(string filename, const sitk::Image& img) {
	sitk::ImageFileWriter writer;
	sitk::Image&& cast_img = sitk::Cast(img, sitk::sitkFloat32);
	writer.SetFileName(filename);
	writer.Execute(cast_img);
}

std::vector<unsigned> find_rotated_size(const vector<unsigned int>& original_size, const vector<double>& rotation_matrix) {
	vector<unsigned int> rotated_size;
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
	return rotated_size;
}

void neg_to_minus1(sitk::Image& img) {
	img = sitk::Cast(sitk::GreaterEqual(img, 0), img.GetPixelID()) * img - sitk::Cast(sitk::Less(img, 0), img.GetPixelID());
}


sitk::Image CalculatePhi(const sitk::Image& image, double beta, double alpha){
	sitk::Image phi = sitk::Cast(sitk::GradientMagnitude(image, false), sitk::sitkFloat64);
	vector<unsigned> size = image.GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	double* phi_buffer = phi.GetBufferAsDouble();

OMP_PARALLEL_FOR
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				double grad = BUF_IDX3D(phi_buffer, xs, ys, zs, xx, yy, zz);
				BUF_IDX3D(phi_buffer, xs, ys, zs, xx, yy, zz) = alpha + (1-alpha) * exp(-beta * grad);
			}
		}
	}
	return phi;
}

double GetImageMinimum(const sitk::Image& image) {
	static sitk::MinimumMaximumImageFilter filter;
	filter.Execute(image);
	return filter.GetMinimum();
}

double GetImageMaximum(const sitk::Image& image) {
	static sitk::MinimumMaximumImageFilter filter;
	filter.Execute(image);
	return filter.GetMaximum();
}

vector<Vec3<int>> ResolvePath(Vec3<int> point, const sitk::Image& distanceMap)
{
	// m_smoothdist 
	vector<uint32_t> size = distanceMap.GetSize();
	int xs = size[0], ys = size[1], zs = size[2];
	const double* distance_buffer = distanceMap.GetBufferAsDouble();

	vector<Vec3<int>> path({ point });

	auto [x, y, z] = point;


	if (BUF_IDX3D(distance_buffer, xs, ys, zs, (int)x, (int)y, (int)z) < 0) return path;

	static vector<double> sqrts = { 0, 1, sqrt(2), sqrt(3) };

	for (int ii = 0; ii < 11111; ++ii) {
		auto [ix, iy, iz] = point;

		if (iz < 2) iz = 2; if (iz >= zs - 2) iz = zs - 3;
		if (ix < 2) ix = 2; if (ix >= xs - 2) ix = xs - 3;
		if (iy < 2) iy = 2; if (iy >= ys - 2) iy = ys - 3;

		
		realnum mmin(0);
		Vec3<int> minimal_direction;
		for (int zo = -1; zo <= 1; ++zo) {
			int zz = iz + zo;
			for (int yo = -1; yo <= 1; ++yo) {
				int yy = iy + yo;
				for (int xo = -1; xo <= 1; ++xo) {
					int xx = ix + xo;
					if (!zo && !yo && !xo) continue;
					double dist = sqrts[abs(xo) + abs(yo) + abs(zo)];
					double d0 = BUF_IDX3D(distance_buffer, xs, ys, zs, (int)ix, (int)iy, (int)iz);
					double d1 = BUF_IDX3D(distance_buffer, xs, ys, zs, xx, yy, zz);
					double val = (d1 - d0) / dist;
					if (val < mmin) {
						mmin = val;
						minimal_direction = { xo,yo,zo };
					}
				}
			}
		}
		point += minimal_direction;
		path.push_back(point);
		if (BUF_IDX3D(distance_buffer, xs, ys, zs, point.x(), point.y(), point.z()) < 1.)
			break;
	}
	return path;
}

double ImageQuantile(const sitk::Image& image, double Q)
{
	int n_pixels = image.GetNumberOfPixels();
	double* vals_buffer = (double*)malloc(n_pixels * sizeof(double));
	memcpy(vals_buffer, image.GetBufferAsDouble(), n_pixels * sizeof(double));
	std::nth_element(vals_buffer, vals_buffer + (int)(n_pixels * Q), vals_buffer + n_pixels);
	double val = vals_buffer[(int)(n_pixels * Q)];
	free(vals_buffer);
	return val;
}

sitk::Image SmartSigmoid(const sitk::Image& image, double qMax, double qMin, double eps)
{
	double q_maxval = ImageQuantile(image, qMax);
	double q_minval = ImageQuantile(image, qMin);
	double alpha = (q_maxval - q_minval) / log((1 / eps - 1) / (1 / 0.95 - 1));
	double beta = log(1 / eps - 1) * alpha + q_minval;
	return sitk::Sigmoid(image, alpha, beta, 1., 0.);
}

std::vector<double> calculateOffsetFromRotation(std::vector<double> rotationMatrix, std::vector<unsigned> imageSize) {
	std::vector<unsigned>  rotated_size = find_rotated_size(imageSize, rotationMatrix);
	sitk::Image sample_img(rotated_size, sitk::sitkInt8);
	sample_img.SetDirection(rotationMatrix);
	std::vector<double> data_center;
	std::transform(imageSize.begin(), imageSize.end(), std::back_inserter(data_center), [](double v) { return (double)v / 2; });

	std::vector<double> sample_center;
	std::transform(rotated_size.begin(), rotated_size.end(), std::back_inserter(sample_center), [](double v) { return v / 2; });
	sample_center = sample_img.TransformContinuousIndexToPhysicalPoint(sample_center);

	std::vector<double> sample_origin;
	std::transform(data_center.begin(), data_center.end(), sample_center.begin(), std::back_inserter(sample_origin), std::minus<double>());
	return sample_origin;
}