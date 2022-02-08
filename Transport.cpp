#include "stdafx.h"
#include "Transport.h"
#include "SimpleITK.h"
#include "Utils.h"
namespace sitk = itk::simple;
using namespace std;


void CTransport::TrInit(sitk::Image &distmap, SVoxImg<SWorkImg<realnum>>& inimap, realnum maxdistance)
{
	/*
		called with:
			distmap = m_liftedEikonal.m_phasefield.m_distance[0]
			inimap = passi
			maxdistance = m_liftedEikonal.currentdistance[0]
	*/
	_PROFILING;
	vector<unsigned> size = distmap.GetSize();
	SVoxImg<SWorkImg<double>> distmap_vox;
	double* distmap_buffer = distmap.GetBufferAsDouble();
	int xs = size[0], ys = size[1], zs = size[2];
	if (inimap.xs != xs) return;
	if (inimap.ys != ys) return;
	if (inimap.zs != zs) return;

	m_isboundary = sitk::Image(size, sitk::sitkInt32);
	int* isboundary_buffer = m_isboundary.GetBufferAsInt32();
	m_transportfunction[0] = sitk::Image(size, sitk::sitkFloat64);
	double* trf0_buffer = m_transportfunction[0].GetBufferAsDouble();
	m_transportfunction[1] = sitk::Image(size, sitk::sitkFloat64);
	double* trf1_buffer = m_transportfunction[1].GetBufferAsDouble();

	realnum mmin(_NO_BOU_*2); 
	for (int zz = 0; zz < zs; ++zz)
		for (int yy = 0; yy < ys; ++yy)
			for (int xx = 0; xx < xs; ++xx) {
				if (inimap[zz][yy][xx] < mmin) mmin = inimap[zz][yy][xx]; // minimal value by plane distance function
				if (distmap_buffer[xx + xs * (yy + ys * zz)] < -0.5) // not calculated
					distmap_buffer[xx + xs * (yy + ys * zz)] = 1.2 * maxdistance; // new
				if (!xx || xx == xs - 1 || !yy || yy == ys - 1 || !zz || zz == zs - 1) // external borders
					distmap_buffer[xx + xs * (yy + ys * zz)] = 1.2 * maxdistance; // new
			}
	m_gx = sitk::Derivative(distmap, 0, 1, false);
	m_gy = sitk::Derivative(distmap, 1, 1, false);
	m_gz = sitk::Derivative(distmap, 2, 1, false);

	m_min = mmin;
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				if (inimap[zz][yy][xx] < _NO_BOU_) { // _NO_BOU_ = not boundary, < _NO_BOU_ = calculated boundary
					BUF_IDX(trf1_buffer, xs, ys, zs, xx, yy, zz) = BUF_IDX(trf0_buffer, xs, ys, zs, xx, yy, zz) = inimap[zz][yy][xx]; // passed arrival plane
					BUF_IDX(isboundary_buffer, xs, ys, zs, xx, yy, zz) = 1; // 1: read only boundary
				}
				else {
					if (!xx || xx == xs - 1 || !yy || yy == ys - 1 || !zz || zz == zs - 1) { // external borders
						BUF_IDX(trf1_buffer, xs, ys, zs, xx, yy, zz) = BUF_IDX(trf0_buffer, xs, ys, zs, xx, yy, zz) = mmin;
						BUF_IDX(isboundary_buffer, xs, ys, zs, xx, yy, zz) = 1; // 1: read only boundary
					}//TODO asda
					else {
						BUF_IDX(trf1_buffer, xs, ys, zs, xx, yy, zz) = BUF_IDX(trf0_buffer, xs, ys, zs, xx, yy, zz) = 0;// -mmin; // initial transport value, is 0 better?
						BUF_IDX(isboundary_buffer, xs, ys, zs, xx, yy, zz) = 0; // at points where calculation is needed

					}
				}
			}
		}
	}
	m_active = 1;
}

void CTransport::TrIterate(int bev)
{

	_PROFILING;
	int targi = bev ? 0 : 1;
	sitk::Image& trf = m_transportfunction[bev];// m_transportfunction[bev];
	sitk::Image& trft = m_transportfunction[targi];// m_transportfunction[targi];
	double* trf_buffer = trf.GetBufferAsDouble();
	double* trft_buffer = trft.GetBufferAsDouble();
	sitk::Image& bound = m_isboundary;
	int* bound_buffer = bound.GetBufferAsInt32();
	vector<uint32_t> size = trf.GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	static const realnum alpha = 0.003;// small value

	double* gx_buffer = m_gx.GetBufferAsDouble();
	double* gy_buffer = m_gy.GetBufferAsDouble();
	double* gz_buffer = m_gz.GetBufferAsDouble();

	#pragma omp parallel for
	for (int zyx = 0; zyx < zs * ys * xs; zyx++) {
		int zz = zyx / (ys * xs);
		if (zz > 0 && zz < zs - 1)
		{
			int yy = (zyx / xs) % ys;
			if (yy > 0 && yy < ys - 1)
			{
				int xx = zyx % xs;
				if (xx > 0 && xx < xs - 1)
				{
					if (BUF_IDX(bound_buffer, xs, ys, zs, xx, yy, zz) != 0) {
						continue; // boundary to completely omit (no transport from uninitialized distance map)
					}

					realnum val;
					int x, y, z;
					realnum gx(BUF_IDX(gx_buffer, xs, ys, zs, xx, yy, zz)), gy(BUF_IDX(gy_buffer, xs, ys, zs, xx, yy, zz)), gz(BUF_IDX(gz_buffer, xs, ys, zs, xx, yy, zz));
					x = gx > 0 ? xx + 1 : xx - 1;
					y = gy > 0 ? yy + 1 : yy - 1;
					z = gz > 0 ? zz + 1 : zz - 1;
					if (gx < 0) gx *= -1;
					if (gy < 0) gy *= -1;
					if (gz < 0) gz *= -1;
					val = gx * BUF_IDX(trf_buffer, xs, ys, zs, x, yy, zz) + gy * BUF_IDX(trf_buffer, xs, ys, zs, xx, y, zz) + gz * BUF_IDX(trf_buffer, xs, ys, zs, xx, yy, z);
					val /= gx + gy + gz + alpha;// zero: at points where calculation is needed, 1 read only boundary points
					BUF_IDX(trft_buffer, xs, ys, zs, xx, yy, zz) = val;  // update wherever calculation is needed
				}
			}
		}
	}

}

void CTransport::TrControl(int nIter)
{
	_PROFILING;
	if (!m_active) return;

	for (int ii = 0; ii < nIter; ++ii)
		TrIterate(ii&1);
}

void CTransport::GetDispSlice(int along, int at, SDisImg& r) // better colorization!
{
	if (!m_active) return;
	realnum zlim(1e-22);
	int bcol(0);
	sitk::Image& trf = m_transportfunction[0];
	double* trf_buffer = trf.GetBufferAsDouble();
	int* bound_buffer = m_isboundary.GetBufferAsInt32();
	vector<uint32_t> size = trf.GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);

	if (along == Talox) {
		if (ys != r.xs || zs != r.ys) {
			r.Clean();
			r.dat = new unsigned long[ys * zs];
			if (!r.dat) return;
			r.xs = ys; r.ys = zs;
		}
		
		int xx(at);
		for (int zz = 0; zz < zs; ++zz) {
			for (int yy = 0; yy < ys; ++yy) {
				if (BUF_IDX(bound_buffer, xs, ys, zs, xx, yy, zz) == -1) // not initialized
					r[zz][yy] = (0xff << 16) + (0xff << 8) + 0xff;
				else { //if (m_isboundary[zz][yy][xx] == 1) // boundary or internal
					if (BUF_IDX(trf_buffer, xs, ys, zs, xx, yy, zz) < -zlim) // -m_min: maximal value
						r[zz][yy] = (bcol << 16) + (bcol << 0) + (int)256*(127 + 0xff * ((0.5* BUF_IDX(trf_buffer, xs, ys, zs, xx, yy, zz) / m_min)));
					else if (BUF_IDX(trf_buffer, xs, ys, zs, xx, yy, zz) > zlim)
						r[zz][yy] = ((int)(127 -0xff * (0.5 * BUF_IDX(trf_buffer, xs, ys, zs, xx, yy, zz) / m_min)) << 16) + (bcol << 8) + bcol;
					else 
						r[zz][yy] = (0xff << 16) + (0xff << 8) + 0xff;
				}

			}
		}
		
	}
	else if (along == Taloy) {
		if (xs != r.xs || zs != r.ys) {
			r.Clean();
			r.dat = new unsigned long[xs * zs];
			if (!r.dat) return;
			r.xs = xs; r.ys = zs;
		}

		int yy(at);
		for (int zz = 0; zz < zs; ++zz) {
			for (int xx = 0; xx < xs; ++xx) {
				if (BUF_IDX(bound_buffer, xs, ys, zs, xx, yy, zz) == -1) // not initialized
					r[zz][xx] = (0xff << 16) + (0xff << 8) + 0xff;
				else { //if (m_isboundary[zz][yy][xx] == 1) // boundary or internal
					if (BUF_IDX(trf_buffer, xs, ys, zs, xx, yy, zz) < -zlim) // -m_min: maximal value
						r[zz][xx] = (bcol << 16) + (bcol << 0) + (int)256*(127 + 0xff * ((0.5 * BUF_IDX(trf_buffer, xs, ys, zs, xx, yy, zz) / m_min)));
					else if (BUF_IDX(trf_buffer, xs, ys, zs, xx, yy, zz) > zlim)
						r[zz][xx] = ((int)(127 -0xff * (0.5 * BUF_IDX(trf_buffer, xs, ys, zs, xx, yy, zz) / m_min)) << 16) + (bcol << 8) + bcol;
					else
						r[zz][xx] = (0xff << 16) + (0xff << 8) + 0xff;
				}
			}
		}

	}
	else if (along == Taloz) {
		if (xs != r.xs || ys != r.ys) {
			r.Clean();
			r.dat = new unsigned long[xs * ys];
			if (!r.dat) return;
			r.xs = xs; r.ys = ys;
		}

		int zz(at);
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				if (BUF_IDX(bound_buffer, xs, ys, zs, xx, yy, zz) == -1) // not initialized
					r[yy][xx] = (0xff << 16) + (0xff << 8) + 0xff;
				else { //if (m_isboundary[zz][yy][xx] == 1) // boundary or internal
					if (BUF_IDX(trf_buffer, xs, ys, zs, xx, yy, zz) < -zlim) // -m_min: maximal value
						r[yy][xx] = (bcol << 16) + (bcol << 0) + (int)256*(127 + 0xff * ((0.5 * BUF_IDX(trf_buffer, xs, ys, zs, xx, yy, zz) / m_min)));
					else if (BUF_IDX(trf_buffer, xs, ys, zs, xx, yy, zz) > zlim)
						r[yy][xx] = ((int)(127 -0xff * (0.5 * BUF_IDX(trf_buffer, xs, ys, zs, xx, yy, zz) / m_min)) << 16) + (bcol << 8) + bcol;
					else
						r[yy][xx] = (0xff << 16) + (0xff << 8) + 0xff;
				}
			}
		}

	}


}

