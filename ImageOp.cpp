#include "StdAfx.h"
#include "ImageOp.h"
#include "math.h"

CImageOp::CImageOp(void)
{
}

CImageOp::~CImageOp(void)
{
}

void CImageOp::GaussFilterImage(int nGauss)
{
	int xs = m_aux.xs, ys = m_aux.ys;
	if (xs != m_loc.xs || ys != m_loc.ys) return;

	int gg = nGauss*2-1;
	for (int kk = 0; kk < gg; ++kk) {
		SWorkImg<realnum> &sou = (kk&1) ? m_aux:m_loc;
		SWorkImg<realnum> &des = (kk&1) ? m_loc:m_aux;
		for (int ii = 1; ii < ys-1; ++ii) {
			for (int jj = 1; jj < xs-1; ++jj) {
				des[ii][jj] = 0.25f*sou[ii][jj];
				des[ii][jj] += 0.125f*sou[ii+1][jj];
				des[ii][jj] += 0.125f*sou[ii-1][jj];
				des[ii][jj] += 0.125f*sou[ii][jj+1];
				des[ii][jj] += 0.125f*sou[ii][jj-1];
				des[ii][jj] += 0.0625f*sou[ii+1][jj+1];
				des[ii][jj] += 0.0625f*sou[ii-1][jj-1];
				des[ii][jj] += 0.0625f*sou[ii-1][jj+1];
				des[ii][jj] += 0.0625f*sou[ii+1][jj-1];
			}
		}
	}
}

void CImageOp::CalcOneComponent(realnum r, SWorkImg<realnum> &comp)
{
	int xs = comp.xs, ys = comp.ys;

	int ir = (int)(r+1.0);
	for (int yy = 0; yy < ys; ++yy) {
		for (int xx = 0; xx < xs; ++xx) {

			m_aux[yy][xx] = 0;
			for (int qq = -ir; qq <= ir; ++qq) {
				int iy = yy+qq;
				if (iy < 0) iy = 0; else if (iy >= ys) iy = ys-1;
				for (int pp = -ir; pp <= ir; ++pp) {
					int ix = xx+pp;
					if (ix < 0) ix = 0; else if (ix >= xs) ix = xs-1;
					if ((realnum)(pp*pp+qq*qq) <= r*r) {
						m_aux[yy][xx] += comp[iy][ix];
					}
				}
			}

		}
	}

	m_aux *= 1.0/r;
	comp = m_aux;

}

void CImageOp::CalculateCircularGradientFlow(realnum r, SWorkImg<realnum> &img)
{
	m_bOk = false;
	int xs = img.xs, ys = img.ys;
	if (!xs || !ys) return;
	img.GetImgHesse(m_gxx,m_gxy,m_gyy);
	m_aux.Set(xs,ys);
	if (xs != m_gxx.xs || ys != m_gxx.ys) return;
	if (xs != m_gxy.xs || ys != m_gxy.ys) return;
	if (xs != m_gyy.xs || ys != m_gyy.ys) return;
	if (xs != m_aux.xs || ys != m_aux.ys) return;

	CalcOneComponent(r,m_gxx);
	CalcOneComponent(r,m_gxy);
	CalcOneComponent(r,m_gyy);
		
	m_bOk = true;
}

int CImageOp::GetDirectionalFlow(realnum r, CVec2 dir, SWorkImg<realnum> &img, SWorkImg<realnum> &out, bool bForced, int expcoef, bool bReverse, int nG)
{
	int xs = img.xs, ys = img.ys;
	if (m_loc.xs != xs || m_loc.ys != ys || bForced) {
		m_loc = img;
		if (!bReverse) m_loc *= -1;
		m_aux = m_loc;
		GaussFilterImage(nG);

		CalculateCircularGradientFlow(r,m_loc);
		if (!m_bOk) return 0;
	}

	{	// orientation score
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				m_aux[yy][xx] = dir.x*dir.x*m_gxx[yy][xx]+2*dir.x*dir.y*m_gxy[yy][xx]+dir.y*dir.y*m_gyy[yy][xx];
				m_aux[yy][xx] -= 0.25*2; // NEW
				if (m_aux[yy][xx] < 0) m_aux[yy][xx] = 0;
			}
		}
	}

	realnum maxv(-1); 
	for (int yy = 0; yy < ys; ++yy) {
		for (int xx = 0; xx < xs; ++xx) {
			if (m_aux[yy][xx] > maxv) maxv = m_aux[yy][xx];
		}
	}

	out = m_aux;
	if (xs != out.xs || ys != out.ys) return 0;
	if (maxv <= 1e-22) return 0;
	out *= 1.0/maxv;

	
	{	// data coefficient
		realnum cf = (realnum) expcoef;
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				out[yy][xx] = 0.001+0.999*exp(-cf*out[yy][xx]);
				//if (out[yy][xx] < 1e-6) out[yy][xx] = 1e-6; //!!!
			}
		}
	}

	return 1;
}


//-----------------------------------------------------------------------------------------



void CImageOp::GauTest(bool bodd)
{
	int xs(m_testimage.xs), ys(m_testimage.ys), zs(m_testimage.zs);

	SVoxImg<SWorkImg<realnum>> &sou = !bodd ? m_testimage:m_testinput;
	SVoxImg<SWorkImg<realnum>> &des = !bodd ? m_testinput:m_testimage;
#pragma omp parallel for
	for (int zz = 1; zz < zs-1; ++zz) {
		for (int yy = 1; yy < ys-1; ++yy) {
			for (int xx = 1; xx < xs-1; ++xx) {
				des[zz][yy][xx] = 0.25f*sou[zz][yy][xx];
				des[zz][yy][xx] += 0.125f*sou[zz+1][yy][xx];
				des[zz][yy][xx] += 0.125f*sou[zz-1][yy][xx];
				des[zz][yy][xx] += 0.125f*sou[zz][yy+1][xx];
				des[zz][yy][xx] += 0.125f*sou[zz][yy-1][xx];
				des[zz][yy][xx] += 0.125f*sou[zz][yy][xx+1];
				des[zz][yy][xx] += 0.125f*sou[zz][yy][xx-1];
			}
		}
		for (int yy = 1; yy < ys-1; ++yy) {
			des[zz][yy][0] = des[zz][yy][1];
			des[zz][yy][xs-1] = des[zz][yy][xs-2];
		}
		for (int xx = 1; xx < xs-1; ++xx) {
			des[zz][0][xx] = des[zz][1][xx];
			des[zz][ys-1][xx] = des[zz][ys-2][xx];
		}
	}
	des[0] = des[1];
	des[zs-1] = des[zs-2];
}

SVoxImg<SWorkImg<realnum>>& CImageOp::CreateTestImage(int xs, int ys, int zs, int expcoef)
{
	m_testimage.Set0(xs,ys,zs);
	m_testinput.Set0(xs,ys,zs);

	realnum s2 = (realnum) 1;
	realnum cf = (realnum) expcoef;


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
					m_testimage[zz][yy][xx] = 0.9f;
					if (zz > cz) m_testimage[zz][yy][xx] -= (dz2)/(Dz*Dz);
				}
				else if (R2 < s2) {
					m_testimage[zz][yy][xx] = 0.9f;
					if (zz < cz) m_testimage[zz][yy][xx] -= 0.25f*(dz2)/(Dz*Dz);
				}
				else m_testimage[zz][yy][xx] = 0.3f;
				//m_testimage[zz][yy][xx] += 0.0001f*((rand()&0xfff)-0x7ff);
			}
		}
	}

	GauTest(false);
	GauTest(true);
	GauTest(false);
	GauTest(true);


	m_testimage.GetGradLen(m_testinput);
	//m_imagegrad = m_testinput;

#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				////m_testinput[zz][yy][xx] = (0.01+m_testinput[zz][yy][xx]); //1.0/(0.000001+m_testinput[zz][yy][xx])
				//m_testinput[zz][yy][xx] = 1.0/(0.000001+0.999999*exp(-cf*m_testinput[zz][yy][xx]));
				m_testinput[zz][yy][xx] = 1.0/(0.01+0.99*exp(-cf*m_testinput[zz][yy][xx]));
			}
		}
	}
	//m_testimage = m_testinput;

	return m_testimage;
}

SWorkImg<realnum>& CImageOp::GetXImageSlice(int ix)
{
	int xs(m_testimage.xs), ys(m_testimage.ys), zs(m_testimage.zs);
	m_loc.Set(ys,zs);
	for (int zz = 0; zz < zs; ++zz)
		for (int yy = 0; yy < ys; ++yy)
			m_loc[zz][yy] = m_testimage[zz][yy][ix];

	return m_loc;
}

void CImageOp::GetXTestDisTrans()
{
	int xs(m_loc.xs), ys(m_loc.ys);

	/*m_gxx = m_loc;
	m_gyy.Set(xs, ys, 0.0);

	for (int ii = 0; ii < m_bandthick; ++ii) {
		SWorkImg<realnum>& sou = !(ii & 1) ? m_gxx : m_gyy;
		SWorkImg<realnum>& des = (ii & 1) ? m_gxx : m_gyy;
		for (int yy = 0 + 1; yy < ys - 1; ++yy) {
			for (int xx = 0 + 1; xx < xs - 1; ++xx) {
				des[yy][xx] = m_loc[yy][xx];
				realnum minn = sou[yy + 1][xx];
				if (sou[yy - 1][xx] < minn) minn = sou[yy - 1][xx];
				if (sou[yy][xx - 1] < minn) minn = sou[yy][xx - 1];
				if (sou[yy][xx + 1] < minn) minn = sou[yy][xx + 1];
				des[yy][xx] += minn;
			}
		}
	}

	m_loc = m_gyy;*/

	std::unordered_set<unsigned long>::const_iterator pe = m_bound.end();
	std::unordered_set<unsigned long>::const_iterator pb = m_bound.begin();

	realnum bandthick = (realnum)m_bandthick;

	for (int yy = 0; yy < ys; ++yy) {
		for (int xx = 0; xx < xs; ++xx) {
			std::unordered_set<unsigned long>::const_iterator pp ;
			int d2m(1 << 30);
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
				if (d2 < d2m) d2m = d2;
			}
			realnum dd = sqrt((realnum)d2m);
			if (dd <= bandthick) m_loc[yy][xx] = m_loc[yy][xx] > 0.5 ? dd:-dd; // actually 1.0 (that is > 0.5)
			//else m_loc[yy][xx] = _NO_BOU_;
			//else m_loc[yy][xx] = 0;
			else m_loc[yy][xx] = m_loc[yy][xx] > 0.5 ? bandthick : -bandthick; // -bandthick try this instead

			//if (!d2m) m_loc[yy][xx] = 0; // temporal (test)
		}
	}

	//m_loc *= 0.5 / bandthick; m_loc += 0.5; // temporal (test)
}

void CImageOp::XTestfill(short x, short y)
{
	std::unordered_set<unsigned long>::const_iterator pe = m_bound.end();
	std::unordered_set<unsigned long>::const_iterator pp = m_bound.find((y<<16)+x);
	short xs = m_loc.xs, ys = m_loc.ys;
	if (pp == pe && m_loc[y][x] > 0.5) {
		m_loc[y][x] = 0;
		if (x < xs-1) XTestfill(x + 1, y);
		if (y < ys-1) XTestfill(x, y + 1);
		if (x > 0) XTestfill(x - 1, y);
		if (y > 0) XTestfill(x, y - 1);
	}

}


void CImageOp::GetXTestIntern()
{
	int xs(m_loc.xs), ys(m_loc.ys);
	m_aux.Set(xs, ys);

	std::unordered_set<unsigned long>::const_iterator pp;
	std::unordered_set<unsigned long>::const_iterator pe = m_bound.end();

	for (int zz = 0; zz < ys; ++zz) {
		int out(2);
		for (int yy = 0; yy < xs; ++yy) {
			pp = m_bound.find((zz << 16) + yy);
			if (pp != pe) { // contour point
				if (out == 2) out = 1; // 1st boundary
				if (out == -2) out = -1;
			}
			else {
				if (out == 1) out = -2;
				if (out == -1) out = 2;
			}
			if (out == 1 || out == -2)
				m_loc[zz][yy] = 1;
			else 
				m_loc[zz][yy] = 0;
		}
	}

	/*m_aux = m_loc;
	for (int zz = 0+1; zz < ys-1; ++zz) {
		for (int yy = 0; yy < xs; ++yy) {
			if (m_aux[zz][yy] < 0.5) {
				if (m_aux[zz + 1][yy] > 0.5 && m_aux[zz - 1][yy] > 0.5)
					m_loc[zz][yy] = 1;
			}
			else {
				if (m_aux[zz + 1][yy] < 0.5 && m_aux[zz - 1][yy] < 0.5)
					m_loc[zz][yy] = 0;
			}
		}
	}
	for (int zz = 0; zz < ys; ++zz)
		m_loc[zz][0] = m_loc[zz][xs-1] = 0;
	for (int yy = 0; yy < xs; ++yy)
		m_loc[0][yy] = m_loc[ys-1][yy] = 0;*/


	//GetXTestDisTrans();

}

void CImageOp::GetXTestBound(int ix, std::vector<CVec3>& out) // contour on plane MAIN 1
{
	int xs(m_testimage.xs), ys(m_testimage.ys), zs(m_testimage.zs);
	out.clear();
	m_bound.clear();
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			int xx = ix;
			if (m_testimage[zz][yy][xx] >= 0.5f) {
				if (m_testimage[zz+1][yy][xx] < 0.5f || m_testimage[zz-1][yy][xx] < 0.5f
				 || m_testimage[zz][yy+1][xx] < 0.5f || m_testimage[zz][yy-1][xx] < 0.5f) {
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

void CImageOp::GetPlaneDistMap(std::unordered_set<unsigned long>& boundset)
{
	int xs(m_testimage.xs), ys(m_testimage.ys), zs(m_testimage.zs);
	m_bound = boundset;

	m_loc.Set(ys, zs, 1);
	XTestfill(0, 0);
	GetXTestDisTrans();
}


SVoxImg<SWorkImg<realnum>>& CImageOp::CImageOp::GetIniMap(int ix) // initial map MAIN 2
{
	m_inimap.Set(m_testimage.xs, m_testimage.ys, m_testimage.zs);
	for (int zz = 0; zz < m_inimap.zs; ++zz) {
		m_inimap[zz] = _NO_BOU_;
	}
	for (int zz = 0; zz < m_inimap.zs; ++zz) {
		for (int yy = 0; yy < m_inimap.ys; ++yy) {
			m_inimap[zz][yy][ix] = m_loc[zz][yy];
		}
	}

	return m_inimap; // all _NO_BOU_ except the true distance map around the boundary curve
}

//-----------------------------------------------------------------------------------------


SVoxImg<SWorkImg<realnum>>& CImageOp::CreateTestImage2(int xs, int ys, int zs, int expcoef)
{
	m_testimage.Set0(xs,ys,zs);
	m_testinput.Set0(xs,ys,zs);

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
					m_testimage[zz][yy][xx] = 0.9f;
					//if (zz > cz) m_testimage[zz][yy][xx] -= (dz2)/(Dz*Dz);
				}
				else m_testimage[zz][yy][xx] = 0.3f;
				
				m_testimage[zz][yy][xx] += 0.0001f*((rand()&0xfff)-0x7ff);
			}
		}
	}

	GauTest(false);
	GauTest(true);
	GauTest(false);
	GauTest(true);


	m_testimage.GetGradLen(m_testinput);
	//m_imagegrad = m_testinput;

#pragma omp parallel for
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				////m_testinput[zz][yy][xx] = (0.01+m_testinput[zz][yy][xx]); //1.0/(0.000001+m_testinput[zz][yy][xx])
				//m_testinput[zz][yy][xx] = 1.0/(0.000001+0.999999*exp(-cf*m_testinput[zz][yy][xx]));
				m_testinput[zz][yy][xx] = 1.0/(0.01+0.99*exp(-cf*m_testinput[zz][yy][xx]));
			}
		}
	}
	//m_testimage = m_testinput;

	return m_testimage;
}

