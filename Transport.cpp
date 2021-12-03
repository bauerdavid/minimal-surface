#include "stdafx.h"
#include "Transport.h"


void CTransport::TrInit(SVoxImg<SWorkImg<realnum>> &distmap, SVoxImg<SWorkImg<realnum>>& inimap, realnum maxdistance)
{
	/*
		called with:
			distmap = m_liftedEikonal.m_phasefield.m_distance[0]
			inimap = passi
			maxdistance = m_liftedEikonal.currentdistance[0]
	*/

	int xs = distmap.xs, ys = distmap.ys, zs = distmap.zs;
	if (inimap.xs != xs) return;
	if (inimap.ys != ys) return;
	if (inimap.zs != zs) return;
	m_gx.Set(xs, ys, zs);
	m_gy.Set(xs, ys, zs);
	m_gz.Set(xs, ys, zs);

	m_isboundary.Set(xs, ys, zs);
	m_transportfunction[0].Set0(xs, ys, zs);
	m_transportfunction[1].Set0(xs, ys, zs);

	realnum mmin(_NO_BOU_*2); 
	for (int zz = 0; zz < zs; ++zz)
		for (int yy = 0; yy < ys; ++yy)
			for (int xx = 0; xx < xs; ++xx) {
				if (inimap[zz][yy][xx] < mmin) mmin = inimap[zz][yy][xx]; // minimal value by plane distance function
				if (distmap[zz][yy][xx] < -0.5) // not calculated
					distmap[zz][yy][xx] = 1.2 * maxdistance; // new
				if (!xx || xx == xs - 1 || !yy || yy == ys - 1 || !zz || zz == zs - 1) // external borders
					distmap[zz][yy][xx] = 1.2 * maxdistance; // new
			}

	distmap.GetGrad(m_gx, m_gy, m_gz);

	m_min = mmin;
	for (int zz = 0; zz < zs; ++zz) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				if (inimap[zz][yy][xx] < _NO_BOU_) { // _NO_BOU_ = not boundary, < _NO_BOU_ = calculated boundary
					m_transportfunction[1][zz][yy][xx] = m_transportfunction[0][zz][yy][xx] = inimap[zz][yy][xx]; // passed arrival plane
					m_isboundary[zz][yy][xx] = 1; // 1: read only boundary
				}
				else {
					if (!xx || xx == xs - 1 || !yy || yy == ys - 1 || !zz || zz == zs - 1) { // external borders
						m_transportfunction[1][zz][yy][xx] = m_transportfunction[0][zz][yy][xx] = mmin;
						m_isboundary[zz][yy][xx] = 1; // 1: read only boundary
					}//TODO asda
					else {
						m_transportfunction[1][zz][yy][xx] = m_transportfunction[0][zz][yy][xx] = 0;// -mmin; // initial transport value, is 0 better?
						m_isboundary[zz][yy][xx] = 0; // at points where calculation is needed
						/* // new
						if (distmap[zz][yy][xx] >= 0.0 
						      && (distmap[zz + 1][yy][xx] < -0.5 || distmap[zz - 1][yy][xx] < -0.5
							   || distmap[zz][yy + 1][xx] < -0.5 || distmap[zz][yy - 1][xx] < -0.5
							   || distmap[zz][yy][xx - 1] < -0.5 || distmap[zz][yy][xx + 1] < -0.5)
							) { // no transport from uninitialized distance map
							m_gx[zz][yy][xx] = m_gy[zz][yy][xx] = m_gz[zz][yy][xx] = 0; 
							}
						*/
					}
				}
			}
		}
	}
	m_active = 1;
}

void CTransport::TrIterate(int bev)
{

	int targi = bev ? 0 : 1;
	SVoxImg<SWorkImg<realnum>>& trf = m_transportfunction[bev];// m_transportfunction[bev];
	SVoxImg<SWorkImg<realnum>>& trft = m_transportfunction[targi];// m_transportfunction[targi];

	SVoxImg<SWorkImg<int>>& bound = m_isboundary;

	int xs(trf.xs), ys(trf.ys), zs(trf.zs);
	realnum alpha = 0.003;// small value

	#pragma omp parallel for
	for (int zz = 1; zz < zs-1; ++zz)
		for (int yy = 1; yy < ys-1; ++yy)
			for (int xx = 1; xx < xs-1; ++xx) {
				if (bound[zz][yy][xx] != 0) {
					continue; // boundary to completely omit (no transport from uninitialized distance map)
				}
				
				realnum val;
				int x, y, z;
				realnum gx(m_gx[zz][yy][xx]), gy(m_gy[zz][yy][xx]), gz(m_gz[zz][yy][xx]);
				x = gx > 0 ? xx + 1 : xx - 1;
				y = gy > 0 ? yy + 1 : yy - 1;
				z = gz > 0 ? zz + 1 : zz - 1;
				if (gx < 0) gx *= -1;
				if (gy < 0) gy *= -1;
				if (gz < 0) gz *= -1;
				val = gx * trf[zz][yy][x] + gy * trf[zz][y][xx] + gz * trf[zz][yy][xx];
				val /= gx + gy + gz + alpha;
				if (!bound[zz][yy][xx]) // zero: at points where calculation is needed, 1 read only boundary points
					trft[zz][yy][xx] = val;  // update wherever calculation is needed
			}

}

void CTransport::TrControl(int nIter)
{
	if (!m_active) return;

	for (int ii = 0; ii < nIter; ++ii)
		TrIterate(ii&1);
}

void CTransport::GetDispSlice(int along, int at, SDisImg& r) // better colorization!
{
	if (!m_active) return;
	realnum zlim(1e-22);
	int bcol(0);
	SVoxImg<SWorkImg<realnum>>& trf = m_transportfunction[0];
	int xs(trf.xs), ys(trf.ys), zs(trf.zs);

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
				if (m_isboundary[zz][yy][xx] == -1) // not initialized
					r[zz][yy] = (0xff << 16) + (0xff << 8) + 0xff;
				else { //if (m_isboundary[zz][yy][xx] == 1) // boundary or internal
					if (trf[zz][yy][xx] < -zlim) // -m_min: maximal value
						r[zz][yy] = (bcol << 16) + (bcol << 0) + (int)256*(127 + 0xff * ((0.5*trf[zz][yy][xx] / m_min)));
					else if (trf[zz][yy][xx] > zlim)
						r[zz][yy] = ((int)(127 -0xff * (0.5 * trf[zz][yy][xx] / m_min)) << 16) + (bcol << 8) + bcol;
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
				if (m_isboundary[zz][yy][xx] == -1) // not initialized
					r[zz][xx] = (0xff << 16) + (0xff << 8) + 0xff;
				else { //if (m_isboundary[zz][yy][xx] == 1) // boundary or internal
					if (trf[zz][yy][xx] < -zlim) // -m_min: maximal value
						r[zz][xx] = (bcol << 16) + (bcol << 0) + (int)256*(127 + 0xff * ((0.5 * trf[zz][yy][xx] / m_min)));
					else if (trf[zz][yy][xx] > zlim)
						r[zz][xx] = ((int)(127 -0xff * (0.5 * trf[zz][yy][xx] / m_min)) << 16) + (bcol << 8) + bcol;
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
				if (m_isboundary[zz][yy][xx] == -1) // not initialized
					r[yy][xx] = (0xff << 16) + (0xff << 8) + 0xff;
				else { //if (m_isboundary[zz][yy][xx] == 1) // boundary or internal
					if (trf[zz][yy][xx] < -zlim) // -m_min: maximal value
						r[yy][xx] = (bcol << 16) + (bcol << 0) + (int)256*(127 + 0xff * ((0.5 * trf[zz][yy][xx] / m_min)));
					else if (trf[zz][yy][xx] > zlim)
						r[yy][xx] = ((int)(127 -0xff * (0.5 * trf[zz][yy][xx] / m_min)) << 16) + (bcol << 8) + bcol;
					else
						r[yy][xx] = (0xff << 16) + (0xff << 8) + 0xff;
				}
			}
		}

	}


}

