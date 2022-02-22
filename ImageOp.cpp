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



void CImageOp::GauTest(bool bodd)
{
	vector<unsigned> size = m_testimage.GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);

	sitk::Image &sou = !bodd ? m_testimage:m_testinput;
	sitk::Image &des = !bodd ? m_testinput:m_testimage;
	des = sitk::Image(sou.GetSize(), sou.GetPixelID());
	double* sou_buffer = sou.GetBufferAsDouble();
	double* des_buffer = des.GetBufferAsDouble();
#pragma omp parallel for
	for (int zz = 1; zz < zs-1; ++zz) {
		for (int yy = 1; yy < ys-1; ++yy) {
			for (int xx = 1; xx < xs-1; ++xx) {
				BUF_IDX(des_buffer, xs, ys, zs, xx, yy, zz) = 0.25f * BUF_IDX(sou_buffer, xs, ys, zs, xx, yy, zz);
				BUF_IDX(des_buffer, xs, ys, zs, xx, yy, zz) += 0.125f* BUF_IDX(sou_buffer, xs, ys, zs, xx, yy, zz+1);
				BUF_IDX(des_buffer, xs, ys, zs, xx, yy, zz) += 0.125f* BUF_IDX(sou_buffer, xs, ys, zs, xx, yy, zz-1);
				BUF_IDX(des_buffer, xs, ys, zs, xx, yy, zz) += 0.125f* BUF_IDX(sou_buffer, xs, ys, zs, xx, yy+1, zz);
				BUF_IDX(des_buffer, xs, ys, zs, xx, yy, zz) += 0.125f* BUF_IDX(sou_buffer, xs, ys, zs, xx, yy-1, zz);
				BUF_IDX(des_buffer, xs, ys, zs, xx, yy, zz) += 0.125f* BUF_IDX(sou_buffer, xs, ys, zs, xx+1, yy, zz);
				BUF_IDX(des_buffer, xs, ys, zs, xx, yy, zz) += 0.125f* BUF_IDX(sou_buffer, xs, ys, zs, xx-1, yy, zz);
			}
		}
		for (int yy = 1; yy < ys-1; ++yy) {
			BUF_IDX(des_buffer, xs, ys, zs, 0, yy, zz) = BUF_IDX(des_buffer, xs, ys, zs, 1, yy, zz);
			BUF_IDX(des_buffer, xs, ys, zs, xs-1, yy, zz) = BUF_IDX(des_buffer, xs, ys, zs, xs-2, yy, zz);
		}
		for (int xx = 1; xx < xs-1; ++xx) {
			BUF_IDX(des_buffer, xs, ys, zs, xx, 0, zz) = BUF_IDX(des_buffer, xs, ys, zs, xx, 1, zz);
			BUF_IDX(des_buffer, xs, ys, zs, xx, ys-1, zz) = BUF_IDX(des_buffer, xs, ys, zs, xx, ys-2, zz);
		}
	}
	memcpy(des_buffer, &BUF_IDX(des_buffer, xs, ys, zs, 0, 0, 1), xs * ys * sizeof(double));
	memcpy(&BUF_IDX(des_buffer, xs, ys, zs, 0, 0, zs-1), &BUF_IDX(des_buffer, xs, ys, zs, 0, 0, zs - 2), xs * ys * sizeof(double));

}

sitk::Image& CImageOp::CreateTestImage(int xs, int ys, int zs, int expcoef)
{
	m_testimage = sitk::Image({ (unsigned)xs,(unsigned)ys, (unsigned)zs }, sitk::sitkFloat64);
	m_testinput = sitk::Image({ (unsigned)xs,(unsigned)ys, (unsigned)zs }, sitk::sitkFloat64);
	double* testimage_buffer = m_testimage.GetBufferAsDouble();

	realnum s2 = (realnum) 1;


	int cx(xs/2+10), cy(ys/2+10), cz(zs/2);
	realnum Dz = (realnum)zs-cz;
#pragma omp parallel for
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
					BUF_IDX(testimage_buffer, xs, ys ,zs, xx, yy, zz) = 0.9f;
					if (zz > cz) BUF_IDX(testimage_buffer, xs, ys, zs, xx, yy, zz) -= (dz2)/(Dz*Dz);
				}
				else if (R2 < s2) {
					BUF_IDX(testimage_buffer, xs, ys, zs, xx, yy, zz) = 0.9f;
					if (zz < cz) BUF_IDX(testimage_buffer, xs, ys, zs, xx, yy, zz) -= 0.25f*(dz2)/(Dz*Dz);
				}
				else BUF_IDX(testimage_buffer, xs, ys, zs, xx, yy, zz) = 0.3f;
				//m_testimage[zz][yy][xx] += 0.0001f*((rand()&0xfff)-0x7ff);
			}
		}
	}
	sitk::Image kernel = sitk::Image({ 3, 3, 3 }, sitk::sitkFloat64);
	double* kernel_buffer = kernel.GetBufferAsDouble();
	double kernel_values[] = { 0,     0, 0,     0, 0.125,     0, 0,     0, 0,
							  0, 0.125, 0, 0.125,  0.25, 0.125, 0, 0.125, 0,
							  0,     0, 0,     0, 0.125,     0, 0,     0, 0 };
	memcpy(kernel_buffer, kernel_values, 3 * 3 * 3 * sizeof(double));
	m_testimage = sitk::Convolution(m_testimage, kernel);
	m_testimage = sitk::Convolution(m_testimage, kernel);
	m_testimage = sitk::Convolution(m_testimage, kernel);
	m_testimage = sitk::Convolution(m_testimage, kernel);

	CreateTestInput(expcoef);

	return m_testimage;
}

sitk::Image& CImageOp::CreateTestInput(double expcoef) {
	realnum cf = (realnum)expcoef;
	m_testinput = sitk::Cast(sitk::GradientMagnitude(m_testimage, false), sitk::sitkFloat64);
	vector<unsigned> size = m_testinput.GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	double* testinput_buffer = m_testinput.GetBufferAsDouble();
	//m_testimage.GetGradLen(m_testinput);
	//m_imagegrad = m_testinput;

#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				////m_testinput[zz][yy][xx] = (0.01+m_testinput[zz][yy][xx]); //1.0/(0.000001+m_testinput[zz][yy][xx])
				//m_testinput[zz][yy][xx] = 1.0/(0.000001+0.999999*exp(-cf*m_testinput[zz][yy][xx]));
				BUF_IDX(testinput_buffer, xs, ys, zs, xx, yy, zz) = 1.0 / (0.01 + 0.99 * exp(-cf * BUF_IDX(testinput_buffer, xs, ys, zs, xx, yy, zz)));
			}
		}
	}
}


void CImageOp::GetXTestDisTrans()
{
	_PROFILING;
	int xs(m_loc.xs), ys(m_loc.ys);

	std::unordered_set<unsigned long>::const_iterator pe = m_bound.end();
	std::unordered_set<unsigned long>::const_iterator pb = m_bound.begin();

	realnum bandthick = (realnum)m_bandthick;

	for (int yy = 0; yy < ys; ++yy) {
		for (int xx = 0; xx < xs; ++xx) {
			std::unordered_set<unsigned long>::const_iterator pp ;
			int max_dist(1 << 30); // big number
			for (pp = pb; pp != pe; ++pp) {
				unsigned long da = *pp;
				int cx = (da & 0xffff);
				int cy = (da >> 16);
				cx -= xx; 
				if (cx < 0) cx *= -1;
				if (cx > m_bandthick) continue;
				cy -= yy;
				if (cy < 0) cy *= -1;
				if (cy > m_bandthick) continue;
				int d2 = cx * cx + cy * cy;
				if (d2 < max_dist) max_dist = d2;
			}
			realnum dd = sqrt((realnum)max_dist);
			if (dd <= bandthick) m_loc[yy][xx] = m_loc[yy][xx] > 0.5 ? dd:-dd; // actually 1.0 (that is > 0.5)
			//else m_loc[yy][xx] = _NO_BOU_;
			//else m_loc[yy][xx] = 0;
			else m_loc[yy][xx] = m_loc[yy][xx] > 0.5 ? bandthick : -bandthick; // -bandthick try this instead

			//if (!max_dist) m_loc[yy][xx] = 0; // temporal (test)
		}
	}

	//m_loc *= 0.5 / bandthick; m_loc += 0.5; // temporal (test)
}

void CImageOp::XTestfill(short x, short y)
{
	_PROFILING;
	std::queue<unsigned long> x_queue;
	std::queue<unsigned long> y_queue;
	x_queue.push(x);
	y_queue.push(y);
	std::unordered_set<unsigned long>::const_iterator pe = m_bound.end();
	short xs = m_loc.xs, ys = m_loc.ys;
	while (!x_queue.empty()) {
		short x = x_queue.front();
		short y = y_queue.front();
		x_queue.pop();
		y_queue.pop();
		std::unordered_set<unsigned long>::const_iterator pp = m_bound.find((y << 16) + x);
		if (pp == pe && m_loc[y][x] > 0.5) {
			m_loc[y][x] = 0;
			if (x < xs - 1) {
				x_queue.push(x + 1);
				y_queue.push(y);
			}
			if (y < ys - 1) {
				x_queue.push(x);
				y_queue.push(y + 1);
			}
			if (x > 0) {
				x_queue.push(x - 1);
				y_queue.push(y);
			}
			if (y > 0) {
				x_queue.push(x );
				y_queue.push(y - 1);
			}
		}
	}
}

void CImageOp::GetXTestBound(int ix, std::vector<CVec3>& out) // contour on plane MAIN 1
{
	vector<unsigned> size = m_testimage.GetSize();
	double* image_buffer = m_testimage.GetBufferAsDouble();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	out.clear();
	m_bound.clear();
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			int xx = ix;
			if (BUF_IDX(image_buffer, xs, ys, zs, xx, yy, zz) >= 0.5f) {
				if (BUF_IDX(image_buffer, xs, ys, zs, xx, yy, zz+1) < 0.5f || BUF_IDX(image_buffer, xs, ys, zs, xx, yy, zz-1) < 0.5f
				 || BUF_IDX(image_buffer, xs, ys, zs, xx, yy+1, zz) < 0.5f || BUF_IDX(image_buffer, xs, ys, zs, xx, yy-1, zz) < 0.5f) {
					 CVec3 in((realnum)xx,(realnum)yy,(realnum)zz);
					 out.push_back(in); // to show paths
					 m_bound.emplace((zz<<16)+yy);
				}
			}
		}
	}

	m_loc.Set(ys,zs,1);
	XTestfill(0, 0);
	GetXTestDisTrans();
}

void CImageOp::GetPlaneDistMap(int distmap_ys, int distmap_zs, std::unordered_set<unsigned long>& boundset)
{
	_PROFILING;
	m_bound = boundset;
	m_loc.Set(distmap_ys, distmap_zs, 1);
	XTestfill(0, 0);
	GetXTestDisTrans();
}


sitk::Image& CImageOp::CImageOp::GetIniMap(int xs, int ys, int zs, int ix) // initial map MAIN 2
{
	_PROFILING;
	m_inimap = sitk::Image({ (unsigned int)xs, (unsigned int)ys, (unsigned int)zs }, sitk::sitkFloat64) + _NO_BOU_;
	double* inimap_buffer = m_inimap.GetBufferAsDouble();
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			BUF_IDX(inimap_buffer, xs, ys, zs, ix, yy, zz) = m_loc[zz][yy];
		}
	}

	return m_inimap; // all _NO_BOU_ except the true distance map around the boundary curve
}

//-----------------------------------------------------------------------------------------


sitk::Image& CImageOp::CreateTestImage2(int xs, int ys, int zs, int expcoef)
{
	m_testimage = sitk::Image({ (unsigned)xs, (unsigned)ys , (unsigned)zs }, sitk::sitkFloat64);
	m_testinput = sitk::Image({ (unsigned)xs, (unsigned)ys , (unsigned)zs }, sitk::sitkFloat64);
	double* image_buffer = m_testimage.GetBufferAsDouble();
	double* input_buffer = m_testinput.GetBufferAsDouble();

	realnum s2 = (realnum) 1;
	realnum cf = (realnum) expcoef;


	int cx(xs/2+0), cy(ys/2), cz(zs/2);
	realnum Dz = (realnum)zs-cz;
#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				realnum dx2 = (realnum) xx - cx;
				dx2 *= dx2;
				realnum dy2 = (realnum) yy - cy;
				dy2 *= dy2;
				realnum dz2 = (realnum) zz - cz - 0.005*(xx - cx)*(xx - cx) - 0.01*(yy - cy)*(yy - cy);
				dz2 *= dz2;

				realnum r2 = dx2/(80*80) + dy2/(60*60) + dz2/(25*25);
				if (r2 < s2) {
					BUF_IDX(image_buffer, xs, ys, zs, xx, yy, zz) = 0.9f;
					//if (zz > cz) m_testimage[zz][yy][xx] -= (dz2)/(Dz*Dz);
				}
				else BUF_IDX(image_buffer, xs, ys, zs, xx, yy, zz) = 0.3f;
				
				BUF_IDX(image_buffer, xs, ys, zs, xx, yy, zz) += 0.0001f*((rand()&0xfff)-0x7ff);
			}
		}
	}

	GauTest(false);
	GauTest(true);
	GauTest(false);
	GauTest(true);

	m_testinput = sitk::GradientMagnitude(m_testimage, false);

	//m_testimage.GetGradLen(m_testinput);
	//m_imagegrad = m_testinput;

#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				////m_testinput[zz][yy][xx] = (0.01+m_testinput[zz][yy][xx]); //1.0/(0.000001+m_testinput[zz][yy][xx])
				//m_testinput[zz][yy][xx] = 1.0/(0.000001+0.999999*exp(-cf*m_testinput[zz][yy][xx]));
				BUF_IDX(input_buffer, xs, ys, zs, xx, yy, zz) = 1.0/(0.01+0.99*exp(-cf* BUF_IDX(input_buffer, xs, ys, zs, xx, yy, zz)));
			}
		}
	}
	//m_testimage = m_testinput;

	return m_testimage;
}

