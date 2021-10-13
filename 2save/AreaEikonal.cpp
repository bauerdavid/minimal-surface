#include "StdAfx.h"
#include "AreaEikonal.h"

CAreaEikonal::CAreaEikonal(void)
{
	m_valid = 0;
}

CAreaEikonal::~CAreaEikonal(void)
{
}

// Image prep&check

void CAreaEikonal::GetIntImage(SWorkImg<realnum> &work, SWorkImg<realnum> &vx, SWorkImg<realnum> &vy)
{
	if (!m_valid) return;

	int xs = m_work.xs, ys = m_work.ys;
	vx.Set(xs,ys);
	vy.Set(xs,ys);

	for (int yy = 0; yy < ys; ++yy) {
		for (int xx = 0; xx < xs; ++xx) {
			if (!xx) vx[yy][xx] = work[yy][xx];
			else vx[yy][xx] = vx[yy][xx-1]+work[yy][xx];
			if (!yy) vy[yy][xx] = work[yy][xx];
			else vy[yy][xx] = vy[yy-1][xx]+work[yy][xx];
		}
	}

	m_valid = 2;
}

void CAreaEikonal::GetDivergence(SWorkImg<realnum> &vx, SWorkImg<realnum> &vy, SWorkImg<realnum> &div)
{
	if (m_valid < 2) return;

	div.GetDiv(vx,vy);
	div *= 0.5;
}

// Phasefield stuff

realnum g_w = 3;

void CAreaEikonal::PhaseInit(IPoi reginit, IPoi arrival)
{
	if (!m_valid) return;
	if (!m_grays && !m_color) return;
	int spacex, spacey;	
	if (m_grays) { spacex = m_work.xs; spacey = m_work.ys; }
	else { spacex = m_workr.xs; spacey = m_workr.ys; }
	//g_cyc = 0;

	m_phasefield.m_field.Set(spacex,spacey,0.0);
	m_phasefield.m_velo.Set(spacex,spacey,0.0);
	m_phasefield.m_gradlen.Set(spacex,spacey,0.0);
	m_phasefield.m_n[0].Set(spacex,spacey,0.0);
	m_phasefield.m_n[1].Set(spacex,spacey,0.0);
	m_phasefield.m_aux.Set(spacex,spacey,0.0);

	m_phasefield.m_distance.Set(spacex,spacey,0.0);
	m_phasefield.m_distgrad[0].Set(spacex,spacey,0.0);
	m_phasefield.m_distgrad[1].Set(spacex,spacey,0.0);

	m_minpath.clear();

	SWorkImg<realnum> &field = m_phasefield.m_field;
	SWorkImg<realnum> &distance = m_phasefield.m_distance;
	for (int yy = 0; yy < field.ys; ++yy) {
		for (int xx = 0; xx < field.xs; ++xx) {
			field[yy][xx] = -1.0;
			distance[yy][xx] = -1.0;
		}
	}

	int hs = (int)(1.5f*g_w/2);
	m_currentdistance = 0;
	for (int yy = reginit.y-hs; yy < reginit.y+hs; ++yy) {
		for (int xx = reginit.x-hs; xx < reginit.x+hs; ++xx) {
			int dx = xx-reginit.x, dy = yy-reginit.y;
			realnum dd = (realnum)(dx*dx+dy*dy);
			if ((int)dd < hs*hs/4) {
				field[yy][xx] = 1.0f;
				dd /= 10000;
				distance[yy][xx] = dd;
				if (m_currentdistance < dd) m_currentdistance = dd;
			}
		}
	}

	for (int ii = 0; ii < 50; ++ii)
		RegularizePhaseField(m_phasefield.m_field,m_phasefield.m_velo);

	m_reference.x = (realnum)reginit.x;
	m_reference.y = (realnum)reginit.y;

	m_distancefrom = reginit;
	m_distanceto = arrival;
	m_resolvready = 0; // arrival bumping

	m_inittype = arrival.x > 0 ? 2:1;
}

void CAreaEikonal::RegularizePhaseField(SWorkImg<realnum> &field, SWorkImg<realnum> &velo)
{
	//++g_cyc;
	int xs = field.xs, ys = field.ys;
	realnum fac = 21.0f/(g_w*g_w), dfac = -(g_w*g_w)/16.0f;

	SWorkImg<realnum> &aux = m_phasefield.m_aux;

	for (int ii = 0; ii < 1; ++ii) { // 1 - 25

		aux.GetLaplace(field);
		velo.GetLaplace(aux);
		velo *= dfac;

		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				velo[yy][xx] -= aux[yy][xx];
			}
		}

		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				realnum fmm = field[yy][xx];
				velo[yy][xx] -= fac*fmm*fmm*fmm;
				velo[yy][xx] += fac*fmm;
			}
		}

		realnum totalfieldweight = 0.00075f;
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				field[yy][xx] += totalfieldweight*velo[yy][xx];
			}
		}

	}
	for (int yy = 0; yy < ys; ++yy) {
		for (int xx = 0; xx < xs; ++xx) {
			velo[yy][xx] = 0;	// zero out here
		}
	}

}

void CAreaEikonal::CalculateFundQuant()
{
	SWorkImg<realnum> &field = m_phasefield.m_field;
	SWorkImg<realnum> &velo = m_phasefield.m_velo;
	SWorkImg<realnum> &gradlen = m_phasefield.m_gradlen;
	SWorkImg<realnum> &nx = m_phasefield.m_n[0];
	SWorkImg<realnum> &ny = m_phasefield.m_n[1];
	int xs = field.xs, ys = field.ys;

	// basic quantities: unit normal components & gradient length
	//nx.GetGradX(field); ny.GetGradY(field);
	field.GetImgCGrad(nx,ny,false);
	for (int yy = 0; yy < ys; ++yy) {
		for (int xx = 0; xx < xs; ++xx) {
			realnum gx = nx[yy][xx], gy = ny[yy][xx];
			realnum len = sqrt(gx*gx+gy*gy);
			if (len < 0.01f) {
				gradlen[yy][xx] = 0; nx[yy][xx] = 0; ny[yy][xx] = 0;
			}
			else {
				realnum ilen = 1.0f/len;
				gradlen[yy][xx] = len; 
				nx[yy][xx] *= ilen; ny[yy][xx] *= ilen;
			}
		}
	}

}

void CAreaEikonal::Iterate()
{
	static int call = 33;
	SWorkImg<realnum> &field = m_phasefield.m_field;
	SWorkImg<realnum> &velo = m_phasefield.m_velo;
	SWorkImg<realnum> &gradlen = m_phasefield.m_gradlen;
	SWorkImg<realnum> &nx = m_phasefield.m_n[0];
	SWorkImg<realnum> &ny = m_phasefield.m_n[1];
	SWorkImg<realnum> &distance = m_phasefield.m_distance;
	int xs = field.xs, ys = field.ys;

	////////////////////////////////////////////////////////////////////
	// fundamental quantities
	////////////////////////////////////////////////////////////////////
	
	CalculateFundQuant(); // for object-level evolution

	////////////////////////////////////////////////////////////////////
	// evolution logic
	////////////////////////////////////////////////////////////////////
	
	{	// Inhomogeneous test

	for (int yy = 1; yy < ys-1; ++yy) {
		for (int xx = 1; xx < xs-1; ++xx) {

			if (field[yy][xx] < 0.95f) { 
				if (distance[yy+1][xx] >= 0 || distance[yy-1][xx] >= 0 //|| distance[yy+1][xx+1] >= 0 || distance[yy-1][xx-1] >= 0
				 || distance[yy][xx+1] >= 0 || distance[yy][xx-1] >= 0) { //|| distance[yy+1][xx-1] >= 0 || distance[yy-1][xx+1] >= 0	

					 // G = f(r)I -> G^-1 = [1/f(r)] I
					realnum eikon(0);

					eikon += -m_work[yy][xx]+1.1f; // path in light
					/*eikon += m_work[yy][xx]+0.15f; // path in dark
					eikon += 3.5f;*/

					if (eikon < 1e-6f) eikon = 1e-6f; // cause: inverse, must be positive (causal)
					eikon *= eikon; // I^4
					eikon = 1.0f/eikon;
					
					velo[yy][xx] = gradlen[yy][xx]*eikon; // |n| = 1
				
				}
			}
			// loop
		}
	}

	}


	////////////////////////////////////////////////////////////////////
	// conditioning (speed limitation)
	////////////////////////////////////////////////////////////////////
//speedlim:
	realnum minv = 1e11f, maxv = 0;

	for (int yy = 0; yy < ys; ++yy) {
		for (int xx = 0; xx < xs; ++xx) {
			realnum fvelo = fabs(velo[yy][xx]);
			if (fvelo > maxv) maxv = fvelo;
			if (fvelo < minv) minv = fvelo;
		}
	}
	maxv += 1e-22f;
	realnum mdatspedmax = 2+0.5f   -2; // set max speed
	maxv = mdatspedmax/maxv;
	m_currentdistance += maxv;
	for (int yy = 0; yy < ys; ++yy) {
		for (int xx = 0; xx < xs; ++xx) {
			velo[yy][xx] *= maxv;
		}
	}
//update:
	for (int yy = 0+1; yy < ys-1; ++yy) {
		for (int xx = 0; xx < xs; ++xx) {
			field[yy][xx] += velo[yy][xx];
			if (field[yy][xx] > 0 && distance[yy][xx] < -0.5f) 
				distance[yy][xx] = m_currentdistance;
		}
	}


	////////////////////////////////////////////////////////////////////
	// phase field regularization with velocity zeroed out
	////////////////////////////////////////////////////////////////////

	static int cl=0; if (!(cl&7)) RegularizePhaseField(field,velo); ++cl;
	if (m_inittype == 2) { if (field[m_distanceto.y][m_distanceto.x] > 0.9f) ++m_resolvready; }

}

void CAreaEikonal::ResolvePath(realnum x,realnum y, bool bClear)
{
	SWorkImg<realnum> &distance = m_phasefield.m_distance;
	SWorkImg<realnum> &gradx = m_phasefield.m_distgrad[0];
	SWorkImg<realnum> &grady = m_phasefield.m_distgrad[1];
	int xs = distance.xs, ys = distance.ys;

	if (bClear) m_minpath.clear();
	int ix((int)x), iy((int)y);

	if (ix < 2 || ix > xs-3) return;
	if (iy < 2 || iy > ys-3) return;
	if (distance[iy][ix] < 0) return;

	CVec2 path(x,y);
	m_minpath.push_back(path);
	distance.GetImgCGrad(gradx,grady,false);

	for (int ii = 0; ii < 11111; ++ii) {
		ix = (int)(path.x+0.5f); iy = (int)(path.y+0.5f);
		CVec2 dir(gradx[iy][ix],grady[iy][ix]);
		dir.Norm(); dir *= 0.64f; path -= dir;
		m_minpath.push_back(path);
		
		dir = path; dir -= m_reference;
		if (dir.x*dir.x+dir.y*dir.y < 0.5f*1.5f) {
			m_minpath.push_back(m_reference);
			break;
		}
	}

}
