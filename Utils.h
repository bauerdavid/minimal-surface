#pragma once
#include <numeric>
#include <algorithm>
#include <sitkImageFileWriter.h>
#include <sitkAdditionalProcedures.h>
#include <sitkImage.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <regex>
#include <unordered_set>
#include <stack>
#include <utility>

//#define SINGLE_THREAD
#ifdef SINGLE_THREAD
#define OMP_PARALLEL_FOR
#define OMP_PARALLEL
#define OMP_FOR
#define OMP_FOR_NOWAIT
#define OMP_CRITICAL(section_name)
#define OMP_ATOMIC
#define OMP_CRITICAL_NO_TAG
#define OMP_PARALLEL_FOR_NUM_THREADS(n_threads)
#define OMP_PARALLEL_NUM_THREADS(n_threads)
#define OMP_FOR_REDUCTION(op, operand)
#define OMP_PARALLEL_FOR_REDUCTION(op, operand)
#else
#define OMP_PARALLEL_FOR __pragma(omp parallel for)
#define OMP_PARALLEL __pragma(omp parallel)
#define OMP_FOR __pragma(omp for)
#define OMP_FOR_NOWAIT __pragma(omp for nowait)
#define OMP_CRITICAL(section_name) __pragma(omp critical(section_name))
#define OMP_ATOMIC __pragma(omp atomic)
#define OMP_CRITICAL_NO_TAG __pragma(omp critical)
#define OMP_PARALLEL_FOR_NUM_THREADS(n_threads) __pragma(omp parallel for num_threads(n_threads))
#define OMP_PARALLEL_NUM_THREADS(n_threads) __pragma(omp parallel num_threads(n_threads))
#define OMP_FOR_REDUCTION(op, operand) __pragma(omp for reduction(op: operand))
#define OMP_PARALLEL_FOR_REDUCTION(op, operand) __pragma(omp parallel for reduction(op: operand))
#endif

#define USE_VECTOR_AS_SET
#define STORE_POINT_AS_INTEGER
#ifdef STORE_POINT_AS_INTEGER
#define POINT3D unsigned long long
#define POINT3D_HASH std::hash<unsigned long long>
#define REMOVABLE_POINT ULLONG_MAX
#else
#define POINT3D Vec3<int>
#define POINT3D_HASH Vec3Hash<int>
#define REMOVABLE_POINT POINT3D(-1, -1, -1)
#endif

#define POINT3D_MAP(val_type) POINT3D, val_type, POINT3D_HASH
#define POINT3D_SET POINT3D, POINT3D_HASH

#define PROFILE_FUNCTIONS

#ifdef PROFILE_FUNCTIONS
class profiler {
	std::chrono::steady_clock::time_point start_time;
	std::string name;
public:
	bool ended_profiling;
	profiler(std::string name, bool push=true);
	~profiler();
	void end_profiling();
	double get_execution_time();
};
struct profile_info {
public:
	size_t n_calls = 0;
	double sum_execution = 0;
	double mean_execution() {
		return sum_execution / n_calls;
	}
	std::unordered_set<std::string> children;
};
class profile_manager {
	static profile_manager& getInstance() {
		static profile_manager instance;
		return instance;
	}
	std::unordered_map<std::string, profile_info> profile_data;
	std::string filename;
	profiler *main_profiler;
	profile_manager() {
		profile_data["__MAIN__"].n_calls = 1;
		main_profiler = new profiler("__MAIN__", false);
		main_profiler->ended_profiling = false;
		call_stack.push("__MAIN__");
	}
	~profile_manager() {
		delete main_profiler;
	}
	std::stack<std::string> call_stack;
public:
	constexpr static int name_w = 100;
	constexpr static int n_calls_w = 7;
	constexpr static int sum_ex_w = 13;
	constexpr static int mean_ex_w = 13;
	static void print_children(std::ostream& os, std::string function, int offset) {
		profile_manager& instance = getInstance();
		auto& info = instance.profile_data[function];
		std::string bullet_point;
		if(offset > 0)
			switch (offset % 3) {
			case 0:
				bullet_point = "- ";
				break;
			case 1:
				bullet_point = "* ";
				break;
			case 2:
				bullet_point = "> ";
				break;
			}
		os << "|"
			<< std::setfill('.') <<  std::setw(name_w) << std::string(offset, ' ') + bullet_point + function << "|"
			<< std::setfill(' ') << std::setw(n_calls_w) << (info.n_calls) << "|"
			<< std::setw(sum_ex_w) << std::setprecision(2) << std::fixed << (info.sum_execution / 1e6) << "|" << std::setw(mean_ex_w)
			<< (info.mean_execution() / 1e6) << "|" << std::endl;

		auto& child_names = info.children;
		for (auto& it = child_names.begin(); it != child_names.end(); it++) {
			std::string child_name = *it;
			if (child_name == "__BLOCK__") continue;
			print_children(os, child_name, offset + 1);
		}
	}

	static void add_func_as_child(std::string parent, std::string child) {
		auto& info = getInstance().profile_data[parent];
		info.children.insert(child);
	}
	static void push_func(std::string func) {
		getInstance().call_stack.push(func);
	}
	static void pop_func() {
		getInstance().call_stack.pop();
	}

	static std::string current_func() {
		return getInstance().call_stack.top();
	}
	static void set_filename(std::string fname) {
		getInstance().filename = fname;
	}
	static void record_info(std::string name, double execution) {
		if (name == "__MAIN__")
			return;
		profile_manager& instance = getInstance();
		profile_info& info = instance.profile_data[name];
		info.n_calls++;
		info.sum_execution += execution;
	}

	static void block_profiling() {
		push_func("__BLOCK__");
	}
	static void unblock_profiling() {
		if (is_blocked())
			pop_func();
	}

	static bool is_blocked() {
		return current_func() == "__BLOCK__";
	}
	static void dump(std::ostream& os) {
		auto& instance = getInstance();
		instance.profile_data["__MAIN__"].sum_execution = instance.main_profiler->get_execution_time();
		os << "|" << std::setfill('=') << std::setw(name_w + n_calls_w + sum_ex_w + mean_ex_w+4) << "|" << std::endl << std::setfill(' ');
		os << std::left;
		os << "|" << std::setw(name_w) << "name"<<"|" << std::setw(n_calls_w) << "#calls"<<"|" << std::setw(sum_ex_w) << "sum time (ms)"<<"|" << std::setw(mean_ex_w) << "avg time (ms)"<<"|" << std::endl;
		os << std::right;
		os << "|" << std::setfill('=') << std::setw(name_w + n_calls_w + sum_ex_w + mean_ex_w+4) << "|" << std::endl << std::setfill(' ');
		os << std::left;
		print_children(os, "__MAIN__", 0);
		os << std::right;
		os << "|" << std::setfill('=') << std::setw(name_w + n_calls_w + sum_ex_w + mean_ex_w+4) << "|" << std::endl << std::setfill(' ');
	}

	static void dump(std::string fname) {
		std::ofstream f_stream;
		f_stream.open(fname);
		dump(f_stream);
		f_stream.close();
	}
};


static inline std::string methodName(const std::string& prettyFunction)
{
	std::string no_template(prettyFunction), prev;
	std::regex template_match("<[^<>]*>");
	do {
		prev = no_template;
		no_template = std::regex_replace(prev, template_match, "");
	} while (no_template != prev);
	size_t lbracket = no_template.rfind("(");
	size_t colons = no_template.substr(0, lbracket).rfind("::");
	size_t begin = no_template.substr(0, colons).rfind(" ") + 1;
	size_t end = lbracket - begin;

	return no_template.substr(begin, end) + "()";
}

#define __METHOD_NAME__ methodName(__FUNCSIG__)
#define _SUB_PROFILE(name) profiler name##_profiler(#name, !profile_manager::is_blocked())
#define _END_SUBPROFILE(name) name##_profiler.end_profiling()
#define _PROFILING profiler _profiler(__METHOD_NAME__, !profile_manager::is_blocked())
#define _END_PROFILING delete _profiler;
#define _DUMP_PROFILE_INFO(output) profile_manager::dump(output);
#define _BLOCK_PROFILING profile_manager::block_profiling()
#define _UNBLOCK_PROFILING profile_manager::unblock_profiling()
#else
#define _SUB_PROFILE(name)
#define _END_SUBPROFILE(name)
#define _PROFILING
#define _END_PROFILING
#define _DUMP_PROFILE_INFO(output)
#define _BLOCK_PROFILING 
#define _UNBLOCK_PROFILING

#endif

#define GET_MACRO(_1, _2, _3, _4, _5, _6, _7, NAME,...) NAME
#define BUF_IDX2D(buffer, xs, ys, xx, yy) buffer[(xx)+(xs)*(yy)]
#define BUF_IDX3D(buffer, xs, ys, zs, xx, yy, zz) buffer[(xx)+(xs)*((yy)+(ys)*(zz))]
//#define BUF_IDX(...) GET_MACRO(__VA_ARGS__, BUF_IDX3D, NOT_DEF, BUF_IDX2D, NOT_DEF, NOT_DEF, NOT_DEF)(__VA_ARGS__)

#include "AreaEikonal.h"
#include "Vec.h"

namespace sitk = itk::simple;

inline POINT3D point_to_representation(int x, int y, int z) {
#ifdef STORE_POINT_AS_INTEGER
	return (POINT3D)x << (2 * 21) | (POINT3D)y << 21 | (POINT3D)z;
#else
	return POINT3D(x, y, z);
#endif
}
template<typename T>
inline std::tuple<T, T, T> representation_to_point(POINT3D representation) {
#ifdef STORE_POINT_AS_INTEGER
	T x = representation >> (2 * 21);
	T y = (representation >> 21) & 0b111111111111111111111;
	T z = representation & 0b111111111111111111111;
	return { x, y, z };
#else
	return { representation.x, representation.y, representation.z };
#endif
}

template<sitk::PixelIDValueEnum pixelIdValueEnum>
struct CType {
	typedef void Type;
};

template<>
struct CType<sitk::sitkFloat32> {
	typedef float Type;
};

template<>
struct CType<sitk::sitkInt32> {
	typedef int Type;
};

template<>
struct CType<sitk::sitkFloat64> {
	typedef double Type;
};

template<typename T>
struct TypeEnumTrait {
public:
	static const sitk::PixelIDValueEnum value = sitk::sitkUnknown;
	static const sitk::InterpolatorEnum interpolator = sitk::sitkLinear;
};

template<>
struct TypeEnumTrait<double> {
public:
	static const sitk::PixelIDValueEnum value = sitk::sitkFloat64;
	static const sitk::InterpolatorEnum interpolator = sitk::sitkLinear;
};

template<>
struct TypeEnumTrait<float> {
public:
	static const sitk::PixelIDValueEnum value = sitk::sitkFloat32;
	static const sitk::InterpolatorEnum interpolator = sitk::sitkBSpline;
};

template<>
struct TypeEnumTrait<int> {
public:
	static const sitk::PixelIDValueEnum value = sitk::sitkInt32;
	static const sitk::InterpolatorEnum interpolator = sitk::sitkNearestNeighbor;
};

template<sitk::PixelIDValueEnum pixelIdValueEnum>
struct PixelManagerTrait {
public:
	/*static void SetPixel(sitk::Image& img, const vector<unsigned int>& idx, CType<pixelIdValueEnum>::Type value);
	static CType<pixelIdValueEnum>::Type GetPixel(sitk::Image& img, const vector<unsigned int>& idx);
	static CType<pixelIdValueEnum>::Type* GetBuffer(sitk::Image& img);
	static void WriteToImage(sitk::Image& img, SVoxImg<SWorkImg<CType<pixelIdValueEnum>::Type>> data);*/
};

template<>
struct PixelManagerTrait<sitk::sitkFloat64> {
public:
	static void SetPixel(sitk::Image& img, const std::vector<unsigned int>& idx, double value) {
		img.SetPixelAsDouble(idx, value);
	}
	static CType<sitk::sitkFloat64>::Type GetPixel(sitk::Image& img, const std::vector<unsigned int>& idx) {
		return img.GetPixelAsDouble(idx);
	}

	static double* GetBuffer(sitk::Image& img) {
		return img.GetBufferAsDouble();
	}

	static const double* GetBuffer(const sitk::Image& img) {
		return img.GetBufferAsDouble();
	}
};

template<>
struct PixelManagerTrait<sitk::sitkFloat32> {
public:
	static void SetPixel(sitk::Image& img, const std::vector<unsigned int>& idx, float value) {
		img.SetPixelAsFloat(idx, value);
	}

	static float GetPixel(sitk::Image& img, const std::vector<unsigned int>& idx) {
		return img.GetPixelAsFloat(idx);
	}

	static float* GetBuffer(sitk::Image& img) {
		return img.GetBufferAsFloat();
	}

	static const float* GetBuffer(const sitk::Image& img) {
		return img.GetBufferAsFloat();
	}
};

template<>
struct PixelManagerTrait<sitk::sitkInt32> {
public:
	static void SetPixel(sitk::Image& img, const std::vector<unsigned int>& idx, int value) {
		img.SetPixelAsInt32(idx, value);
	}
	static int GetPixel(sitk::Image& img, const std::vector<unsigned int>& idx) {
		return img.GetPixelAsInt32(idx);
	}

	static int* GetBuffer(sitk::Image& img) {
		return img.GetBufferAsInt32();
	}

	static const int* GetBuffer(const sitk::Image& img) {
		return img.GetBufferAsInt32();
	}
};

template<>
struct PixelManagerTrait<sitk::sitkUInt8> {
public:
	static void SetPixel(sitk::Image& img, const std::vector<unsigned int>& idx, unsigned char value) {
		img.SetPixelAsUInt8(idx, value);
	}
	static unsigned char GetPixel(sitk::Image& img, const std::vector<unsigned int>& idx) {
		return img.GetPixelAsUInt8(idx);
	}

	static unsigned char* GetBuffer(sitk::Image& img) {
		return img.GetBufferAsUInt8();
	}

	static const unsigned char* GetBuffer(const sitk::Image& img) {
		return img.GetBufferAsUInt8();
	}
};

template<typename T>
std::pair<Vec3<double>, Vec3<double>> best_plane_from_points(const std::vector<std::vector<T>>& c)
{
	// copy coordinates to  matrix in Eigen format
	size_t num_atoms = c.size();
	Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > coord(3, num_atoms);
	for (size_t i = 0; i < num_atoms; ++i) {
		coord.col(i)[0] = c[i][0];
		coord.col(i)[1] = c[i][1];
		coord.col(i)[2] = c[i][2];
	}
	// calculate centroid
	Vec3<double> centroid(coord.row(0).mean(), coord.row(1).mean(), coord.row(2).mean());

	// subtract centroid
	coord.row(0).array() -= centroid.x();
	coord.row(1).array() -= centroid.y();
	coord.row(2).array() -= centroid.z();

	// we only need the left-singular matrix here
	//  http://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
	auto svd = coord.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
	auto U = svd.matrixU().data();
	Vec3<double> meeting_plane_normal(U[6], U[7], U[8]);
	return std::make_pair(centroid, meeting_plane_normal);
}

void rotate(const std::vector<double>& rotation_matrix, const std::vector<double> value, std::vector<double>& dest);

template <typename T> int sgn(T val) {
	return (int)(T(0) < val) - (int)(val < T(0));
}

template<typename T>
void cross_product(const std::vector<T> vec1, const std::vector<T> vec2, std::vector<T>& out) {
	if (vec1.size() != 3 || vec2.size() != 3) {
		return;
	}
	out.push_back(vec1[1] * vec2[2] - vec1[2] * vec2[1]);
	out.push_back(vec1[2] * vec2[0] - vec1[0] * vec2[2]);
	out.push_back(vec1[0] * vec2[1] - vec1[1] * vec2[0]);
}

template<typename T>
std::vector<double> rotation_matrix_from_vectors(std::vector<T>& vec1, std::vector<T>& vec2) {
	std::vector<T> vec1_sqr;
	transform(vec1.begin(), vec1.end(), back_inserter(vec1_sqr), [](T& v) { return v * v; });
	std::vector<T> vec2_sqr;
	transform(vec2.begin(), vec2.end(), back_inserter(vec2_sqr), [](T& v) { return v * v; });

	T sum1 = accumulate(vec1_sqr.begin(), vec1_sqr.end(), 0.0);
	T sum2 = accumulate(vec2_sqr.begin(), vec2_sqr.end(), 0.0);
	transform(vec1.begin(), vec1.end(), vec1.begin(), [sum1](T& c) { return c / sum1; });
	transform(vec2.begin(), vec2.end(), vec2.begin(), [sum2](T& c) { return c / sum2; });
	std::vector<double> v;
	cross_product(vec1, vec2, v);
	double c = inner_product(vec1.begin(), vec1.end(), vec2.begin(), 0.0);
	double norm = sqrt(inner_product(v.begin(), v.end(), v.begin(), 0.0));
	std::vector<double> kmat_plus_eye = std::vector<double>({
		1, -v[2], v[1],
		v[2], 1, -v[0],
		-v[1], v[0], 1 });
	std::vector<double> kmat_2 = std::vector<double>({
		(-v[2] * v[2] - v[1] * v[1]),				 (v[0] * v[1]),				   (v[2] * v[0]),
					   (v[1] * v[0]), (-v[2] * v[2] - v[0] * v[0]),				   (v[2] * v[1]),
					   (v[0] * v[2]),				 (v[1] * v[2]), (-v[0] * v[0] - v[1] * v[1]) });
	transform(kmat_2.begin(), kmat_2.end(), kmat_2.begin(), [c, norm](T& val) {return val * ((1 - c) / (norm * norm)); });
	std::vector<double> rotation_matrix = std::vector<double>(9);
	transform(kmat_plus_eye.begin(), kmat_plus_eye.end(), kmat_2.begin(), rotation_matrix.begin(), plus<double>());
	return rotation_matrix;
}

std::vector<unsigned> find_rotated_size(const std::vector<unsigned int>& original_size, const std::vector<double>& rotation_matrix);

template
<sitk::InterpolatorEnum interpolator = sitk::sitkLinear>
sitk::Image resample_img(sitk::Image src, sitk::Image& dst, const std::vector<double>& rotation_matrix, bool inverse = false, bool useNearestNeighborExtrapolator = false) {
	std::vector<unsigned> data_size = src.GetSize();
	std::vector<unsigned> sample_size = find_rotated_size(data_size, rotation_matrix);
	return resample_img<interpolator>(src, dst, rotation_matrix, sample_size, inverse, useNearestNeighborExtrapolator);
}

template
<sitk::InterpolatorEnum interpolator = sitk::sitkLinear>
sitk::Image resample_img(sitk::Image src, sitk::Image& dst, const std::vector<double>& rotation_matrix, std::vector<unsigned int> sample_size, bool inverse = false, bool useNearestNeighborExtrapolator = false) {
	std::vector<unsigned> data_size = src.GetSize();
	sitk::Image sample_img(sample_size, src.GetPixelID());
	std::vector<double> original_direction = src.GetDirection();
	std::vector<double> original_origin = src.GetOrigin();
	if (inverse)
		src.SetDirection(rotation_matrix);
	else
		sample_img.SetDirection(rotation_matrix);
	std::vector<double> sample_center;
	std::transform(sample_size.begin(), sample_size.end(), std::back_inserter(sample_center), [](double v) { return v / 2; });
	sample_center = sample_img.TransformContinuousIndexToPhysicalPoint(sample_center);

	std::vector<double> data_center;
	std::transform(data_size.begin(), data_size.end(), std::back_inserter(data_center), [](double v) { return (double)v / 2; });
	data_center = src.TransformContinuousIndexToPhysicalPoint(data_center);

	std::vector<double> sample_origin;
	std::transform(data_center.begin(), data_center.end(), sample_center.begin(), std::back_inserter(sample_origin), std::minus<double>());
	sample_img.SetOrigin(sample_origin);
	
	dst = sitk::Resample(src, sample_img, sitk::Transform(), interpolator, -1.0, sitk::sitkUnknown, useNearestNeighborExtrapolator);
	return sample_img;
}




void save_image(std::string filename, const sitk::Image& img);
/*
template<typename T, class TypeEnum = TypeEnumTrait<T>, class PixelManager = PixelManagerTrait<T>>
void save_slice(string filename, SVoxImg <SWorkImg<T>>& data, int xslice) {
	int ys(data.ys), zs(data.zs);
	sitk::Image img(ys, zs, TypeEnum::value);
	T* buffer = PixelManager::GetBuffer(img);
	int xx = xslice;
	for (int zz = 0; zz < zs; zz++) {
		for (int yy = 0; yy < ys; yy++) {
			buffer[yy + ys * zz] = data[zz][yy][xx];
			//PixelManager::SetPixel(img, { (unsigned int)yy, (unsigned int)zz }, data[zz][yy][xx]);
		}
	}
	save_image(filename, img);
}*/

void neg_to_minus1(sitk::Image& img);

template
<sitk::PixelIDValueEnum pixelID=sitk::sitkUInt8, class PixelManager = PixelManagerTrait<pixelID>>
sitk::Image create_6_neighbors_SE() {
	sitk::Image se({ 3, 3, 3 }, pixelID);
	for (int zz = 0; zz < 3; zz++) {
		for (int yy = 0; yy < 3; yy++) {
			for (int xx = 0; xx < 3; xx++) {
				if (xx % 2 + yy % 2 + zz % 2 == 2) {
					PixelManager::SetPixel(se, { (unsigned)xx, (unsigned)yy, (unsigned)zz }, 1);
				}
			}
		}
	}
	return se;
}

template
<sitk::PixelIDValueEnum pixelID = sitk::sitkUInt8, class PixelManager = PixelManagerTrait<pixelID>>
sitk::Image create_19_neighbors_SE() {
	sitk::Image se({ 3, 3, 3 }, pixelID);
	for (int zz = 0; zz < 3; zz++) {
		for (int yy = 0; yy < 3; yy++) {
			for (int xx = 0; xx < 3; xx++) {
				if (xx % 2 + yy % 2 + zz % 2 > 0) {
					PixelManager::SetPixel(se, { (unsigned)xx, (unsigned)yy, (unsigned)zz }, 1);
				}
			}
		}
	}
	return se;
}

sitk::Image CalculatePhi(const sitk::Image& image, double beta = 14., double alpha = 0.01);

template<sitk::PixelIDValueEnum pixelID>
sitk::Image GetImageSlice(const sitk::Image& image, int direction, int slice) {
	vector<unsigned> size = image.GetSize();
	int idxs[3];
	idxs[direction] = slice;
	int loop_direction1 = direction == 0 ? 1 : 0;
	int loop_direction2 = direction == 2 ? 1 : 2;
	sitk::Image slice_image(size[loop_direction1], size[loop_direction2], image.GetPixelID());
	const CType<pixelID>::Type* input_buffer = PixelManagerTrait<pixelID>::GetBuffer(image);
	CType<pixelID>::Type* slice_buffer = PixelManagerTrait<pixelID>::GetBuffer(slice_image);

	for (int i = 0; i < size[loop_direction2]; i++) {
		for (int j = 0; j < size[loop_direction1]; j++) {
			idxs[loop_direction1] = j;
			idxs[loop_direction2] = i;
			BUF_IDX2D(slice_buffer, size[loop_direction1], size[loop_direction2], j, i) = BUF_IDX3D(input_buffer, size[0], size[1], size[2], idxs[0], idxs[1], idxs[2]);
		}
	}
	return slice_image;
}

double GetImageMinimum(const sitk::Image& image);

double GetImageMaximum(const sitk::Image& image);

std::vector<Vec3<int>> ResolvePath(Vec3<int> point, const sitk::Image& distanceMap);