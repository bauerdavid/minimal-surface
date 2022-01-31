#pragma once
#include <utility>
#include "AreaEikonal.h"
#include <numeric>
#include <algorithm>
#include <sitkImageFileWriter.h>
#include <sitkAdditionalProcedures.h>
#include <sitkImage.h>
#include <Eigen/Dense>
#include <Eigen/Core>

#define BUF_IDX(buffer, xs, ys, zs, xx, yy, zz) buffer[(xx) + (xs)*((yy)+(ys)*(zz))]


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

	static void WriteToImage(sitk::Image& img, SVoxImg<SWorkImg<double>>& data) {
		double* buffer = img.GetBufferAsDouble();
		for (int zz = 0; zz < data.zs; zz++) {
			for (int yy = 0; yy < data.ys; yy++) {
				for (int xx = 0; xx < data.xs; xx++) {
					buffer[xx + data.xs * (yy + data.ys * zz)] = data[zz][yy][xx];
				}
			}
		}

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
};

template<typename T>
std::pair<IPoi3<double>, IPoi3<double>> best_plane_from_points(const std::vector<std::vector<T>>& c)
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

void rotate(const std::vector<double> rotation_matrix, const std::vector<double> value, std::vector<double>& dest);

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

void find_rotated_size(std::vector<unsigned int>& original_size, std::vector<double>& rotation_matrix, std::vector<unsigned int>& rotated_size);

template<typename T, sitk::PixelIDValueEnum pixelIdValueEnum = TypeEnumTrait<T>::value, class PixelManager = PixelManagerTrait<pixelIdValueEnum>>
sitk::Image vox_img_2_sitk(SVoxImg<SWorkImg<T>>& data) {
	int xs(data.xs), ys(data.ys), zs(data.zs);
	sitk::Image img({ (unsigned)xs, (unsigned)ys, (unsigned)zs }, pixelIdValueEnum);
	//PixelManager.WriteToImage(img, data);
	CType<pixelIdValueEnum>::Type* buffer = PixelManager::GetBuffer(img);
	for (unsigned int zz = 0; zz < zs; zz++) {
		for (unsigned int yy = 0; yy < ys; yy++) {
			for (unsigned int xx = 0; xx < xs; xx++) {
				//buffer[xx + xs * (yy + ys * zz)] = data[zz][yy][xx];
				PixelManager::SetPixel(img, { xx, yy, zz }, data[zz][yy][xx]);
			}
		}
	}
	return img;
}

template<typename T, sitk::PixelIDValueEnum pixelIdValueEnum = TypeEnumTrait<T>::value, class PixelManager = PixelManagerTrait<pixelIdValueEnum>>
sitk::Image work_img_2_sitk(SWorkImg<T>& data) {
	unsigned int xs(data.xs), ys(data.ys);
	sitk::Image img(xs, ys, pixelIdValueEnum);
	T* buffer = PixelManager::GetBuffer(img);
	for (unsigned int yy = 0; yy < ys; yy++) {
		for (unsigned int xx = 0; xx < xs; xx++) {
			//buffer[xx + xs * yy] = data[yy][xx];
			PixelManager::SetPixel(img, { xx, yy }, data[yy][xx]);
		}
	}
	return img;
}
template
<sitk::InterpolatorEnum interpolator = sitk::sitkLinear>
sitk::Image resample_img(sitk::Image& src, sitk::Image& dst, std::vector<double>& rotation_matrix, std::vector<unsigned int> sample_size, bool inverse = false, bool useNearestNeighborExtrapolator = false) {
	std::vector<unsigned> data_size = src.GetSize();
	if (sample_size.size() == 0)
		find_rotated_size(data_size, rotation_matrix, sample_size);
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

template<typename T, sitk::PixelIDValueEnum pixelIdValueEnum = TypeEnumTrait<T>::value, sitk::InterpolatorEnum interpolatorEnum = TypeEnumTrait<T>::interpolator, class PixelManager = PixelManagerTrait<pixelIdValueEnum>>
sitk::Image resample_vox_img(SVoxImg<SWorkImg<T>>& data, SVoxImg<SWorkImg<T>>& out, std::vector<double>& rotation_matrix, std::vector<unsigned int> sample_size, bool inverse=false, bool useNearestNeighborExtrapolator=false) {
	unsigned int xs(data.xs), ys(data.ys), zs(data.zs);
	std::vector<unsigned int> data_size = { (unsigned int)xs, (unsigned int)ys, (unsigned int)zs };
	sitk::Image data_img = vox_img_2_sitk<T, pixelIdValueEnum, PixelManager>(data);
	sitk::Image resampled;

	sitk::Image sample_image = resample_img<interpolatorEnum>(data_img, resampled, rotation_matrix, sample_size, inverse, useNearestNeighborExtrapolator);
	//double pixVal = resampled.GetPixelAsDouble({ 40, 40, 40 });
	sitk_2_vox_img<T, PixelManager>(resampled, out);
	return sample_image;
}

template<typename T, class PixelManager = PixelManagerTrait<TypeEnumTrait<T>::value>>
void sitk_2_vox_img(sitk::Image& sitk_img, SVoxImg<SWorkImg<T>>& vox_img) {
	vector<unsigned int> sample_size = sitk_img.GetSize();
	vox_img.Set0(sample_size[0], sample_size[1], sample_size[2]);
	for (unsigned int z = 0; z < sample_size[2]; z++) {
		for (unsigned int y = 0; y < sample_size[1]; y++) {
			for (unsigned int x = 0; x < sample_size[0]; x++) {
				T newval = PixelManager::GetPixel(sitk_img, { x, y, z });
				vox_img[z][y][x] = newval;
			}
		}
	}
}


void save_image(std::string filename, sitk::Image& img);
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

template<typename T, typename K=T, sitk::PixelIDValueEnum pixelIdValueEnum = TypeEnumTrait<K>::value, class PixelManager = PixelManagerTrait<pixelIdValueEnum>>
void save_vox_img(std::string filename, SVoxImg < SWorkImg<T>>& data) {
	unsigned int xs(data.xs), ys(data.ys), zs(data.zs);
	sitk::Image img(xs, ys, zs, pixelIdValueEnum);
	K* buffer = PixelManager::GetBuffer(img);
	for (unsigned int zz = 0; zz < zs; zz++) {
		for (unsigned int yy = 0; yy < ys; yy++) {
			for (unsigned int xx = 0; xx < xs; xx++) {
				//buffer[xx + xs * (yy + ys * zz)] = data[zz][yy][xx];
				PixelManager::SetPixel(img, { xx, yy, zz }, data[zz][yy][xx]);
			}
		}
	}
	save_image(filename, img);
}

template<>
inline void save_vox_img<double>
(std::string filename, SVoxImg < SWorkImg<double>>& data) {
	unsigned int xs(data.xs), ys(data.ys), zs(data.zs);
	sitk::Image img(xs, ys, zs, sitk::sitkFloat32);
	float* buffer = img.GetBufferAsFloat();
	for (unsigned int zz = 0; zz < zs; zz++) {
		for (unsigned int yy = 0; yy < ys; yy++) {
			for (unsigned int xx = 0; xx < xs; xx++) {
				//buffer[xx + xs * (yy + ys * zz)] = data[zz][yy][xx];
				PixelManagerTrait<sitk::sitkFloat32>::SetPixel(img, { xx, yy, zz }, data[zz][yy][xx]);
			}
		}
	}
	save_image(filename, img);
}


template<typename T, sitk::PixelIDValueEnum pixelIdValueEnum = TypeEnumTrait<T>::value, class PixelManager = PixelManagerTrait<pixelIdValueEnum>>
void save_work_img(std::string filename, SWorkImg<T>& data) {
	unsigned int xs(data.xs), ys(data.ys);
	sitk::Image img(xs, ys, pixelIdValueEnum);
	for (unsigned int yy = 0; yy < ys; yy++) {
		for (unsigned int xx = 0; xx < xs; xx++) {
			PixelManager::SetPixel(img, { xx, yy }, data[yy][xx]);
		}
	}
	save_image(filename, img);
}

template<>
inline void save_work_img<double>
(std::string filename, SWorkImg<double>& data) {
	unsigned int xs(data.xs), ys(data.ys);
	sitk::Image img(xs, ys, TypeEnumTrait<float>::value);
	for (unsigned int yy = 0; yy < ys; yy++) {
		for (unsigned int xx = 0; xx < xs; xx++) {
			PixelManagerTrait<sitk::sitkFloat32>::SetPixel(img, { xx, yy }, data[yy][xx]);
		}
	}
	save_image(filename, img);
}

template<typename T>
void neg_to_minus1(SVoxImg<SWorkImg<T>>& data) {
	int xs(data.xs), ys(data.ys), zs(data.zs);
	for (int zz = 0; zz < zs; zz++) {
		for (int yy = 0; yy < ys; yy++) {
			for (int xx = 0; xx < xs; xx++) {
				if (data[zz][yy][xx] < 0) {
					data[zz][yy][xx] = -1;
				}
			}
		}
	}
}

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