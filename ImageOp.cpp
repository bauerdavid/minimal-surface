#include "stdafx.h"
#include "ImageOp.h"
#include "math.h"
#include <queue>
#include "Utils.h"
#include <sitkImageOperators.h>

using namespace std;

CImageOp::CImageOp(void)
{
}

CImageOp::~CImageOp(void)
{
}




//-----------------------------------------------------------------------------------------
sitk::Image CImageOp::Blur(const sitk::Image& img) {
	vector<unsigned> size = img.GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	sitk::Image blurred(size, sitk::sitkFloat64);
	const double* sou_buffer = img.GetBufferAsDouble();
	double* des_buffer = blurred.GetBufferAsDouble();
OMP_PARALLEL_FOR
	for (int zz = 1; zz < zs - 1; ++zz) {
		for (int yy = 1; yy < ys - 1; ++yy) {
			for (int xx = 1; xx < xs - 1; ++xx) {
				BUF_IDX3D(des_buffer, xs, ys, zs, xx, yy, zz) = 0.25f * BUF_IDX3D(sou_buffer, xs, ys, zs, xx, yy, zz);
				BUF_IDX3D(des_buffer, xs, ys, zs, xx, yy, zz) += 0.125f * BUF_IDX3D(sou_buffer, xs, ys, zs, xx, yy, zz + 1);
				BUF_IDX3D(des_buffer, xs, ys, zs, xx, yy, zz) += 0.125f * BUF_IDX3D(sou_buffer, xs, ys, zs, xx, yy, zz - 1);
				BUF_IDX3D(des_buffer, xs, ys, zs, xx, yy, zz) += 0.125f * BUF_IDX3D(sou_buffer, xs, ys, zs, xx, yy + 1, zz);
				BUF_IDX3D(des_buffer, xs, ys, zs, xx, yy, zz) += 0.125f * BUF_IDX3D(sou_buffer, xs, ys, zs, xx, yy - 1, zz);
				BUF_IDX3D(des_buffer, xs, ys, zs, xx, yy, zz) += 0.125f * BUF_IDX3D(sou_buffer, xs, ys, zs, xx + 1, yy, zz);
				BUF_IDX3D(des_buffer, xs, ys, zs, xx, yy, zz) += 0.125f * BUF_IDX3D(sou_buffer, xs, ys, zs, xx - 1, yy, zz);
			}
		}
		for (int yy = 1; yy < ys - 1; ++yy) {
			BUF_IDX3D(des_buffer, xs, ys, zs, 0, yy, zz) = BUF_IDX3D(des_buffer, xs, ys, zs, 1, yy, zz);
			BUF_IDX3D(des_buffer, xs, ys, zs, xs - 1, yy, zz) = BUF_IDX3D(des_buffer, xs, ys, zs, xs - 2, yy, zz);
		}
		for (int xx = 1; xx < xs - 1; ++xx) {
			BUF_IDX3D(des_buffer, xs, ys, zs, xx, 0, zz) = BUF_IDX3D(des_buffer, xs, ys, zs, xx, 1, zz);
			BUF_IDX3D(des_buffer, xs, ys, zs, xx, ys - 1, zz) = BUF_IDX3D(des_buffer, xs, ys, zs, xx, ys - 2, zz);
		}
	}
	memcpy(des_buffer, &BUF_IDX3D(des_buffer, xs, ys, zs, 0, 0, 1), xs * ys * sizeof(double));
	memcpy(&BUF_IDX3D(des_buffer, xs, ys, zs, 0, 0, zs - 1), &BUF_IDX3D(des_buffer, xs, ys, zs, 0, 0, zs - 2), xs * ys * sizeof(double));
	return blurred;
}


sitk::Image CImageOp::CreateTestImage()
{
	int xs(XS_), ys(YS_), zs(ZS_);
	sitk::Image test_image = sitk::Image({ (unsigned)xs,(unsigned)ys, (unsigned)zs }, sitk::sitkFloat64);
	double* testimage_buffer = test_image.GetBufferAsDouble();

	realnum s2 = (realnum) 1;


	int cx(xs/2+10), cy(ys/2+10), cz(zs/2);
	realnum Dz = (realnum)zs-cz;
OMP_PARALLEL_FOR
	for (int zz = 0; zz < zs; ++zz) {
		realnum dz2 = (realnum) zz - cz;
		dz2 *= dz2;
		for (int yy = 0; yy < ys; ++yy) {
			realnum dy2 = (realnum) yy - cy;
			dy2 *= dy2;
			for (int xx = 0; xx < xs; ++xx) {
				realnum dx2 = (realnum) xx - cx;
				dx2 *= dx2;
				realnum r2 = dx2/(80*80) + dy2/(50*40) + dz2/(30*30);
				realnum R2 = dx2/(30*30) + dy2/(30*30) + dz2/(50*50);
				if (r2 < s2) {
					BUF_IDX3D(testimage_buffer, xs, ys ,zs, xx, yy, zz) = 0.9f;
					if (zz > cz) BUF_IDX3D(testimage_buffer, xs, ys, zs, xx, yy, zz) -= (dz2)/(Dz*Dz);
				}
				else if (R2 < s2) {
					BUF_IDX3D(testimage_buffer, xs, ys, zs, xx, yy, zz) = 0.9f;
					if (zz < cz) BUF_IDX3D(testimage_buffer, xs, ys, zs, xx, yy, zz) -= 0.25f*(dz2)/(Dz*Dz);
				}
				else BUF_IDX3D(testimage_buffer, xs, ys, zs, xx, yy, zz) = 0.3f;
			}
		}
	}
	sitk::Image kernel = sitk::Image({ 3, 3, 3 }, sitk::sitkFloat64);
	double* kernel_buffer = kernel.GetBufferAsDouble();
	double kernel_values[] = { 0,     0, 0,     0, 0.125,     0, 0,     0, 0,
							  0, 0.125, 0, 0.125,  0.25, 0.125, 0, 0.125, 0,
							  0,     0, 0,     0, 0.125,     0, 0,     0, 0 };
	memcpy(kernel_buffer, kernel_values, 3 * 3 * 3 * sizeof(double));
	test_image = sitk::Convolution(test_image, kernel);
	test_image = sitk::Convolution(test_image, kernel);
	test_image = sitk::Convolution(test_image, kernel);
	test_image = sitk::Convolution(test_image, kernel);

	//CreateTestInput(expcoef);

	return test_image;
}

sitk::Image CImageOp::CreateTestImage2()
{
	int xs(XS_), ys(YS_), zs(ZS_);
	sitk::Image test_image = sitk::Image({ (unsigned)xs, (unsigned)ys , (unsigned)zs }, sitk::sitkFloat64);
	double* image_buffer = test_image.GetBufferAsDouble();

	realnum s2 = (realnum)1;

	int cx(xs / 2 + 0), cy(ys / 2), cz(zs / 2);
	realnum Dz = (realnum)zs - cz;
OMP_PARALLEL_FOR
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				realnum dx2 = (realnum)xx - cx;
				dx2 *= dx2;
				realnum dy2 = (realnum)yy - cy;
				dy2 *= dy2;
				realnum dz2 = (realnum)zz - cz - 0.005 * (xx - cx) * (xx - cx) - 0.01 * (yy - cy) * (yy - cy);
				dz2 *= dz2;

				realnum r2 = dx2 / (80 * 80) + dy2 / (60 * 60) + dz2 / (25 * 25);
				if (r2 < s2) {
					BUF_IDX3D(image_buffer, xs, ys, zs, xx, yy, zz) = 0.9f;
					//if (zz > cz) m_testimage[zz][yy][xx] -= (dz2)/(Dz*Dz);
				}
				else BUF_IDX3D(image_buffer, xs, ys, zs, xx, yy, zz) = 0.3f;

				BUF_IDX3D(image_buffer, xs, ys, zs, xx, yy, zz) += 0.0001f * ((rand() & 0xfff) - 0x7ff);
			}
		}
	}
	sitk::Image kernel = sitk::Image({ 3, 3, 3 }, sitk::sitkFloat64);
	double* kernel_buffer = kernel.GetBufferAsDouble();
	double kernel_values[] = { 0,     0, 0,     0, 0.125,     0, 0,     0, 0,
							  0, 0.125, 0, 0.125,  0.25, 0.125, 0, 0.125, 0,
							  0,     0, 0,     0, 0.125,     0, 0,     0, 0 };
	memcpy(kernel_buffer, kernel_values, 3 * 3 * 3 * sizeof(double));
	test_image = sitk::Convolution(test_image, kernel);
	test_image = sitk::Convolution(test_image, kernel);
	test_image = sitk::Convolution(test_image, kernel);
	test_image = sitk::Convolution(test_image, kernel);

	return test_image;
}
