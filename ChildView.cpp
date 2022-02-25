
// ChildView.cpp : implementation of the CChildView class
//

#include "stdafx.h"
#include "MinArea.h"
#include "ChildView.h"
#include "Utils.h"
#include "AreaEikonal.h"
#ifdef _DEBUG
#define new DEBUG_NEW
#endif
#include <iostream>
#include <fstream>
#include <sitkImage.h>
#include <sitkAdditionalProcedures.h>
#include <omp.h>
#include "FMSigned.h"
#define DISTANCE_ITERATION 1
#define DISTANCE_MEAN_ITERATION 2
#define PLANE_PHASEFIELD_ITERATION 3
#define TRANSPORT_FUNCTION_ITERATION 4
#define DONE_ITERATION 5

#define EXPCOEF_RANGE 100
using namespace std;
namespace sitk = itk::simple;
// CChildView

CChildView::CChildView()
{
	m_pControl = 0;
	m_threadactivated = 0;
	m_bispoint = 0;
	m_zsee = 0;
	m_ysee = 0;
	m_xsee = 0;

	m_bsee = false;
	m_curvsee = 0;
	m_btransportview = false;

	m_valid = 0;
	m_expfac = EXP_COEF_DEF;//+10
}

CChildView::~CChildView()
{
}


BEGIN_MESSAGE_MAP(CChildView, CWnd)
	ON_WM_PAINT()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_RBUTTONUP()
END_MESSAGE_MAP()



// CChildView message handlers

BOOL CChildView::PreCreateWindow(CREATESTRUCT& cs) 
{
	if (!CWnd::PreCreateWindow(cs))
		return FALSE;

	cs.dwExStyle |= WS_EX_CLIENTEDGE;
	cs.style &= ~WS_BORDER;
	cs.lpszClass = AfxRegisterWndClass(CS_HREDRAW|CS_VREDRAW|CS_DBLCLKS, 
		::LoadCursor(NULL, IDC_ARROW), reinterpret_cast<HBRUSH>(COLOR_WINDOW+1), NULL);

	if (!m_pControl) {
		m_pControl = new CControlDlg;
		m_pControl->Create(IDD_DIALOG1);
		m_pControl->ShowWindow(SW_SHOW);
		m_pControl->m_pView = this;
	}

	return TRUE;
}

void CChildView::DispImage(CPaintDC &dc, SDisImg &disp, int offx, int offy)
{
	int xs = disp.xs;
	int ys = disp.ys;
	void *ptr = disp.dat;
	BITMAPINFOHEADER bmi;
	memset(&bmi,0,sizeof(bmi));
	bmi.biWidth = xs;
	bmi.biHeight = -ys;
	bmi.biBitCount = 32;
	bmi.biCompression = BI_RGB;
	bmi.biPlanes = 1;
	bmi.biSize = sizeof(BITMAPINFOHEADER);

	int q = ::StretchDIBits(dc.m_hDC,offx,offy,xs,ys,0,0,xs,ys,ptr,(BITMAPINFO*)&bmi,DIB_RGB_COLORS,SRCCOPY);
	if (q == GDI_ERROR) { DWORD err = GetLastError(); }
}

CPen g_bluepen(PS_SOLID,1,0xff0000);
CPen g_greenpen(PS_SOLID,1,0x00ff00);
CPen g_redpen(PS_SOLID,1,0xff);
CPen g_yellowpen(PS_SOLID,2,0xffff);
CPen g_orangepen(PS_SOLID,2,0x7fff);

unsigned GetDepthColor(realnum val)
{
	int ival((int)val);
	while (ival > 1199) ival -= 1200;
	if (ival < 0) ival = 0;

	unsigned red(0), green(0), blue(0);
	if  (ival < 200) {
		blue += (int)(ival/2) +155;
	} // bluemax
	else if  (ival < 400) {
		blue += 200 +55;
		green += (ival-200) +55;
	} // bluemax + greenmax
	else if  (ival < 500) {
		green += 200 +55;
		blue += (499-ival)*2 +55;
	} // greenmax
	else if  (ival < 700) {
		green += 200 +55;
		red += ival-499 +55;
	} // greenmax + redmax
	else if  (ival < 900) {
		red += 200 +55;
		green += 899-ival +55;
	} // redmax
	else if  (ival < 1000) {
		red += 200 +55;
		blue += (ival-899)*2 +55;
	} // redmax + bluemax
	else {
		red += 1199-ival +55;
		blue += 1199-ival +55;
	}
	
	return red+(green<<8)+(blue<<16);
}

void CChildView::Z2Dir(int z, int &dx, int &dy)
{
	realnum a = z*2*PI/ZS_;
	dx = (int)(30*cos(a));
	dy = (int)(30*sin(a));
}


void CChildView::OnPaint() 
{
	CPaintDC dc(this); // device context for painting
			
	// TODO: Add your message handler code here
	int offx = 0, offy = 0;
	if (m_disp.xs) {
		DispImage(dc,m_disp,offx,offy);
		offx += m_disp.xs+1;
		if (m_dispd1.xs) DispImage(dc,m_dispd1,offx,offy);
		offx += m_dispd1.xs+1;
		if (m_dispd2.xs) {
			DispImage(dc,m_dispd2,offx,offy);
			std::vector<POINT3D> &bound = m_liftedEikonal.m_boundcontour;
			int si = bound.size();
			if (si) {
				for (int ii = 0; ii < si; ++ii) {
					auto [x, y, z] = representation_to_point<int>(bound[ii]);
					//if(m_xsee == x)
						dc.SetPixelV(offx+y,z,0xffff);
				}
			}
		}
		offx += m_dispd2.xs+1;
		if (m_dispd3.xs) DispImage(dc,m_dispd3,offx,offy); // no dispd3
		
		CPen *oldpen = dc.SelectObject(&g_bluepen);
		dc.MoveTo(m_disp.xs,m_zsee); dc.LineTo(offx,m_zsee);
		dc.SelectObject(&g_greenpen);
		dc.MoveTo(0,m_ysee); dc.LineTo(m_disp.xs,m_ysee);
		dc.MoveTo(m_ysee+offx-m_dispd2.xs,0); dc.LineTo(m_ysee+offx-m_dispd2.xs,m_dispd2.ys);
		dc.SelectObject(&g_redpen);
		dc.MoveTo(m_xsee,0); dc.LineTo(m_xsee,m_disp.ys);
		dc.MoveTo(m_xsee+m_disp.xs+1,0); dc.LineTo(m_xsee+m_disp.xs+1,m_dispd1.ys);
		dc.SelectObject(oldpen);

		//offx = 0; offy = m_disp.ys-1; if (m_dislic.xs) DispImage(dc,m_dislic,offx,offy);
	}	

	vector<unsigned> size = m_liftedEikonal.m_phasefield.m_field[0].GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);

	if (m_threadactivated)  {
		//char txt[100];
		//sprintf_s(txt, 99, "%d %d         ", g_ign, g_ian);
		//dc.TextOutA(300, 10, txt);
		for (int ii = 0; ii < 2; ++ii) {
			double* field_buffer = m_liftedEikonal.m_phasefield.m_field[ii].GetBufferAsDouble();
			int* meeting_plane_buffer = m_liftedEikonal.m_phasefield.meeting_plane_positions.GetBufferAsInt32();

			int zsee = 0;
			if (m_pControl) {
				zsee = m_zsee;
			}

			int pasi = m_liftedEikonal.m_minpath[ii][00].size();

			if (m_threadactivated != PLANE_PHASEFIELD_ITERATION) {
				if (m_disp.xs && !m_btransportview) {
					for (int yy = 0; yy < ys; ++yy) {
						for (int xx = 0; xx < xs; ++xx) {
							if (m_threadactivated == DISTANCE_ITERATION || m_threadactivated == DONE_ITERATION) {
								if (abs(m_liftedEikonal.m_phasefield.plane_normal.x*xx+ m_liftedEikonal.m_phasefield.plane_normal.y*yy+ m_liftedEikonal.m_phasefield.plane_normal.z*zsee+ m_liftedEikonal.m_phasefield.plane_offset) < 0.5) {
									dc.SetPixelV(xx, yy, 0xff00ff);
								}
								else if (meeting_plane_buffer[xx + xs * (yy + ys * zsee)] > 0) {
									dc.SetPixelV(xx, yy, 0xff0000);
								}
								else if (field_buffer[xx + xs * (yy + ys * zsee)] > -0.9f && field_buffer[xx + xs * (yy + ys * zsee)] < 0.3f)
									dc.SetPixelV(xx, yy, !ii ? 0xff : 0xff00);
							}
						}
					}
				}
			}

			if (m_dispd1.xs) for (int zz = 0; zz < zs; ++zz) {
				for (int xx = 0; xx < xs; ++xx) {
					if (m_threadactivated  == DISTANCE_ITERATION || m_threadactivated == DONE_ITERATION) {
						if (abs(m_liftedEikonal.m_phasefield.plane_normal.x * xx + m_liftedEikonal.m_phasefield.plane_normal.y * m_ysee + m_liftedEikonal.m_phasefield.plane_normal.z * zz + m_liftedEikonal.m_phasefield.plane_offset) < 0.5) {
							dc.SetPixelV(m_disp.xs + 1 + xx, zz, 0xff00ff);
						}
						else if (meeting_plane_buffer[xx + xs * (m_ysee + ys * zz)] > 0) {
							dc.SetPixelV(m_disp.xs + 1 + xx, zz, 0xff0000);
						}
						else if (field_buffer[xx + xs * (m_ysee + ys * zz)] > -0.9f && field_buffer[xx + xs * (m_ysee + ys * zz)] < 0.3f)
							dc.SetPixelV(m_disp.xs+1+xx,zz,!ii?0xff:0xff00);
					}
					/*else if (!pasi) {
						realnum tocc(dismap[zsee][yy][xx]); // ARRIVA!!
						if (tocc >= 0) 
							dc.SetPixelV(xx,yy,GetDepthColor(tocc/16));
					}*/
				}
			}
			if (m_dispd2.xs) { // xslice
				for (int zz = 0; zz < zs; ++zz) {
					for (int yy = 0; yy < ys; ++yy) {
						if (m_threadactivated == DISTANCE_ITERATION || m_threadactivated == DONE_ITERATION) {
							if (abs(m_liftedEikonal.m_phasefield.plane_normal.x * m_xsee + m_liftedEikonal.m_phasefield.plane_normal.y * yy + m_liftedEikonal.m_phasefield.plane_normal.z * zz + m_liftedEikonal.m_phasefield.plane_offset) < 0.5) {
								dc.SetPixelV(2 * (m_disp.xs + 1) + yy, zz, 0xff00ff);
							}
							else if (meeting_plane_buffer[m_xsee + xs * (yy + ys * zz)] > 0) {
								dc.SetPixelV(2 * (m_disp.xs + 1) + yy, zz, 0xff0000);
							}
							else if (field_buffer[m_xsee + xs * (yy + ys * zz)] > -0.9f && field_buffer[m_xsee + xs * (yy + ys * zz)] < 0.3f)
								dc.SetPixelV(2 * (m_disp.xs + 1) + yy, zz, !ii ? 0xff : 0xff00);
						}
						if (m_threadactivated == PLANE_PHASEFIELD_ITERATION) {
							SWorkImg<realnum>& field2D = m_liftedEikonal.m_inicountourCalculator.m_field2D;
							if (field2D[zz][yy] > 0.5f)
								dc.SetPixelV(2 * (m_disp.xs + 1) + yy, zz, 0xff);
							else if(field2D[zz][yy] < -0.5f)
								dc.SetPixelV(2 * (m_disp.xs + 1) + yy, zz, 0xff0000);
						}
						/*else if (!pasi) {
							realnum tocc(dismap[zsee][yy][xx]); // ARRIVA!!
							if (tocc >= 0)
								dc.SetPixelV(xx,yy,GetDepthColor(tocc/16));
						}*/
					}
				}
			}
			if (pasi) {
				CPen* oldpen = dc.SelectObject(&g_yellowpen);

				int8_t* path_x_buffer = m_path_image_x.GetBufferAsInt8();
				int8_t* path_y_buffer = m_path_image_y.GetBufferAsInt8();
				int8_t* path_z_buffer = m_path_image_z.GetBufferAsInt8();
				if (BUF_IDX(path_z_buffer, xs, ys, zs, 0, 0, m_zsee) == -100) {
					for (int yy = 0; yy < ys; yy++) {
						for (int xx = 0; xx < xs; xx++) {
							BUF_IDX(path_z_buffer, xs, ys, zs, xx, yy, m_zsee) = 0;
						}
					}
					for (int jj = 0; jj < MAXMINPATH; ++jj) {
						//CPen *oldpen = dc.SelectObject(&g_yellowpen);
						std::vector<CVec3>& minpath = m_liftedEikonal.m_minpath[ii][jj];
						pasi = minpath.size();
						if (!pasi) break;
						for (int ii = 0; ii < pasi; ++ii) {
							if(BUF_IDX(path_z_buffer, xs, ys, zs, (int)minpath[ii].x, (int)minpath[ii].y, m_zsee) < 1)
								BUF_IDX(path_z_buffer, xs, ys, zs, (int)minpath[ii].x, (int)minpath[ii].y, m_zsee) = minpath[ii].z < m_zsee ? -1 : 1;
						}
					}
				}
				for (int yy = 0; yy < ys; yy++) {
					for (int xx = 0; xx < xs; xx++) {
						int8_t val = BUF_IDX(path_z_buffer, xs, ys, zs, xx, yy, m_zsee);
						if (val) {
							int color = val < 0 ? 0x7fff : 0xffff;
							dc.SetPixelV(xx, yy, color);
						}
					}
				}

				if (BUF_IDX(path_y_buffer, xs, ys, zs, 0, m_ysee, 0) == -100) {
					for (int zz = 0; zz < zs; zz++) {
						for (int xx = 0; xx < xs; xx++) {
							BUF_IDX(path_y_buffer, xs, ys, zs, xx, m_ysee, zz) = 0;
						}
					}
					for (int jj = 0; jj < MAXMINPATH; ++jj) {
						//CPen *oldpen = dc.SelectObject(&g_yellowpen);
						std::vector<CVec3>& minpath = m_liftedEikonal.m_minpath[ii][jj];
						pasi = minpath.size();
						if (!pasi) break;
						for (int ii = 0; ii < pasi; ++ii) {
							if (BUF_IDX(path_y_buffer, xs, ys, zs, (int)minpath[ii].x, m_ysee, (int)minpath[ii].z) < 1)
								BUF_IDX(path_y_buffer, xs, ys, zs, (int)minpath[ii].x, m_ysee, (int)minpath[ii].z) = minpath[ii].y < m_ysee ? -1 : 1;
						}
					}
				}
				for (int zz = 0; zz < zs; zz++) {
					for (int xx = 0; xx < xs; xx++) {
						int8_t val = BUF_IDX(path_y_buffer, xs, ys, zs, xx, m_ysee, zz);
						if (val) {
							int color = val < 0 ? 0x7fff : 0xffff;
							dc.SetPixelV(xx + m_disp.xs + 1, zz, color);
						}
					}
				}

				if (BUF_IDX(path_x_buffer, xs, ys, zs, m_xsee, 0, 0) == -100) {
					for (int zz = 0; zz < zs; zz++) {
						for (int yy = 0; yy < ys; yy++) {
							BUF_IDX(path_x_buffer, xs, ys, zs, m_xsee, yy, zz) = 0;
						}
					}
					for (int jj = 0; jj < MAXMINPATH; ++jj) {
						//CPen *oldpen = dc.SelectObject(&g_yellowpen);
						std::vector<CVec3>& minpath = m_liftedEikonal.m_minpath[ii][jj];
						pasi = minpath.size();
						if (!pasi) break;
						for (int ii = 0; ii < pasi; ++ii) {
							if (BUF_IDX(path_x_buffer, xs, ys, zs, m_xsee, (int)minpath[ii].y, (int)minpath[ii].z) < 1)
								BUF_IDX(path_x_buffer, xs, ys, zs, m_xsee, (int)minpath[ii].y, (int)minpath[ii].z) = minpath[ii].x < m_xsee ? -1 : 1;
						}
					}
				}
				for (int zz = 0; zz < zs; zz++) {
					for (int yy = 0; yy < ys; yy++) {
						int8_t val = BUF_IDX(path_x_buffer, xs, ys, zs, m_xsee, yy, zz);
						if (val) {
							int color = val < 0 ? 0x7fff : 0xffff;
							dc.SetPixelV(yy + 2 * (m_disp.xs + 1), zz, color);
						}
					}
				}
				
				
				

					//dc.MoveTo((int)minpath[0].x, (int)minpath[0].y);
					

					//dc.MoveTo((int)minpath[0].x + m_disp.xs + 1, (int)minpath[0].z);
					
					//dc.MoveTo((int)minpath[0].y + 2 * (m_disp.xs + 1), (int)minpath[0].z);

				dc.SelectObject(oldpen);
			}
		}


		{
			double* field_buffer = m_liftedEikonal.m_phasefield.m_field[0].GetBufferAsDouble();

			int xx(0), yy(0), zz(0);
			dc.MoveTo(xx,300-(int)(10*field_buffer[xx + xs * (m_ysee + ys * m_zsee)]));
			for (xx = 1; xx < xs; ++xx) dc.LineTo(xx*5,300-(int)(10*field_buffer[xx + xs * (m_ysee + ys * m_zsee)]));
			dc.MoveTo(yy,350-(int)(10*field_buffer[m_xsee + xs * (yy + ys * m_zsee)]));
			for (yy = 1; yy < ys; ++yy) dc.LineTo(yy*5,350-(int)(10*field_buffer[m_xsee + xs * (yy + ys * m_zsee)]));
			dc.MoveTo(zz,400-(int)(10*field_buffer[m_xsee + xs * (m_ysee + ys * zz)]));
			for (zz = 1; zz < zs; ++zz) dc.LineTo(zz*5,400-(int)(10*field_buffer[m_xsee + xs * (m_ysee + ys * zz)]));

		}

	}
			

	if (m_disp.xs) {
		if (m_bispoint > 0) {
			if (abs(m_zsee-m_start_point.z) < 3) dc.FillSolidRect(m_start_point.x-1, m_start_point.y-1,3,3,0xff);
			if (abs(m_ysee- m_start_point.y) < 3) dc.FillSolidRect(m_start_point.x+m_disp.xs, m_start_point.z-1,3,3,0xff);
			if (abs(m_xsee- m_start_point.x) < 3) dc.FillSolidRect(m_start_point.y+2*m_disp.xs+1, m_start_point.z -1,3,3,0xff);
		}
		if (m_bispoint > 1) {
			if (abs(m_zsee- m_end_point.z) < 3) dc.FillSolidRect(m_end_point.x-1, m_end_point.y-1,3,3,0xffff00);
			if (abs(m_ysee- m_end_point.y) < 3) dc.FillSolidRect(m_end_point.x+m_disp.xs, m_end_point.z-1,3,3,0xffff00);
			if (abs(m_xsee- m_end_point.x) < 3) dc.FillSolidRect(m_end_point.y+2*m_disp.xs+1, m_end_point.z-1,3,3,0xffff00);
		}
	}


	// Do not call CWnd::OnPaint() for painting messages
}

UINT BackgroundThread(LPVOID params)
{
	omp_set_nested(1);
	omp_set_dynamic(1);
	CChildView* view = (CChildView*)params;
	int cyc = 0;
	extern bool g_modeswitch;
	bool dumped = false;
	double start = omp_get_wtime();
	ofstream f_stream;
	f_stream.open("Y:/BIOMAG/shortest path/" INFO_FILENAME);
#ifdef DEBUG_STATES
	save_image("Y:/BIOMAG/shortest path/interm_imgs/ph0_data.tif", view->m_liftedEikonal.m_phasefield.m_data);
	save_image("Y:/BIOMAG/shortest path/interm_imgs/ph0_testimage.tif", view->m_imageOp.GetTestImage());
	save_image("Y:/BIOMAG/shortest path/interm_imgs/ph0_dist0.tif", view->m_liftedEikonal.m_phasefield.m_distance[0]);
	save_image("Y:/BIOMAG/shortest path/interm_imgs/ph0_dist1.tif", view->m_liftedEikonal.m_phasefield.m_distance[1]);
	save_image("Y:/BIOMAG/shortest path/interm_imgs/ph0_field0.tif", view->m_liftedEikonal.m_phasefield.m_field[0]);
	save_image("Y:/BIOMAG/shortest path/interm_imgs/ph0_field1.tif", view->m_liftedEikonal.m_phasefield.m_field[1]);
#endif
	while (view->m_threadactivated) {

		if (view->m_threadactivated == DISTANCE_ITERATION) { // 3D
			if (!view->m_liftedEikonal.m_phasefield.m_bdone)
			{
				view->m_liftedEikonal.m_phasefield.Iterate(g_modeswitch);
				if (view->m_liftedEikonal.m_phasefield.m_bdone) {
#ifdef DEBUG_STATES
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph1_flow_idx.tif", view->m_liftedEikonal.m_phasefield.m_flow_idx);
#endif

					view->m_liftedEikonal.m_phasefield.ExtractMeetingPlane();
					view->m_liftedEikonal.m_phasefield.CombineDistance();
					//view->m_liftedEikonal.rotation_matrix = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
					view->m_liftedEikonal.m_rotated_phasefield.Initialize(view->m_liftedEikonal.m_phasefield, view->m_liftedEikonal.m_phasefield.rotation_matrix);
					view->m_liftedEikonal.m_phasefield.SmoothDistances();

					vector<double> plane_center_physical = { (view->m_liftedEikonal.m_phasefield.plane_center.x), (view->m_liftedEikonal.m_phasefield.plane_center.y), (view->m_liftedEikonal.m_phasefield.plane_center.z) };
					vector<double> plane_center_transformed;
					plane_center_transformed = view->m_liftedEikonal.m_rotated_phasefield.m_sample_image.TransformPhysicalPointToContinuousIndex(plane_center_physical);
					view->m_liftedEikonal.m_rotated_phasefield.m_plane_slice = plane_center_transformed[0];
#ifndef ITERATE_ROTATED
					view->m_liftedEikonal.m_rotated_phasefield.CalculateAlignedCombinedDistance(view->m_start_point.x, view->m_end_point.x);
					view->m_liftedEikonal.m_rotated_phasefield.SmoothDistances();
#endif
					view->m_liftedEikonal.m_rotated_phasefield.m_bdone = false;
					
					// calculate data center
					vector<unsigned int> size = view->m_liftedEikonal.m_phasefield.m_distance[0].GetSize();
					double xs(size[0]);
					double ys(size[1]);
					double zs(size[2]);
					vector<unsigned int> rot_size = view->m_liftedEikonal.m_rotated_phasefield.m_distance[0].GetSize();

#ifdef DEBUG_STATES
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph1_combined_distance.tif", view->m_liftedEikonal.m_phasefield.m_combined_distance);
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph1_distance_0.tif", view->m_liftedEikonal.m_phasefield.m_distance[0]);
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph1_distance_1.tif", view->m_liftedEikonal.m_phasefield.m_distance[1]);
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph1_r_combined_distance.tif", view->m_liftedEikonal.m_rotated_phasefield.m_combined_distance);
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph1_r_distance0.tif", view->m_liftedEikonal.m_rotated_phasefield.m_distance[0]);
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph1_r_distance1.tif", view->m_liftedEikonal.m_rotated_phasefield.m_distance[1]);
#endif

					f_stream << "size (x, y, z):";
					f_stream << " " << rot_size[0];
					f_stream << " " << rot_size[1];
					f_stream << " " << rot_size[2];
					f_stream << endl;
					f_stream << "plane center: " << plane_center_transformed[0] << endl;
#ifndef ITERATE_ROTATED
					view->m_threadactivated = DISTANCE_MEAN_ITERATION;
#endif
					view->Invalidate(FALSE);
					
				}
			}
			else if (!view->m_liftedEikonal.m_rotated_phasefield.m_bdone) {
				view->m_liftedEikonal.m_rotated_phasefield.Iterate(g_modeswitch);
				if (view->m_liftedEikonal.m_rotated_phasefield.m_bdone) {
					//view->m_liftedEikonal.m_rotated_phasefield.CombineDistance();
					vector<double> start_point({ view->m_start_point.x, view->m_start_point.y, view->m_start_point.z });
					vector<double> end_point({ view->m_end_point.x, view->m_end_point.y, view->m_end_point.z });
					start_point = view->m_liftedEikonal.m_phasefield.m_sample_image.TransformContinuousIndexToPhysicalPoint(start_point);
					start_point = view->m_liftedEikonal.m_rotated_phasefield.m_sample_image.TransformPhysicalPointToContinuousIndex(start_point);
					end_point = view->m_liftedEikonal.m_phasefield.m_sample_image.TransformContinuousIndexToPhysicalPoint(end_point);
					end_point = view->m_liftedEikonal.m_rotated_phasefield.m_sample_image.TransformPhysicalPointToContinuousIndex(end_point);
					view->m_liftedEikonal.m_rotated_phasefield.CalculateAlignedCombinedDistance(start_point[0], end_point[0]);
					view->m_liftedEikonal.m_rotated_phasefield.SmoothDistances();
					/*sitk::MinimumMaximumImageFilter minmax;
					minmax.Execute(view->m_liftedEikonal.m_rotated_phasefield.m_distance[0]);
					double maxval = minmax.GetMaximum();
					minmax.Execute(view->m_liftedEikonal.m_rotated_phasefield.m_distance[1]);
					maxval = max(maxval, minmax.GetMaximum());
					view->m_liftedEikonal.m_rotated_phasefield.m_distance[0] /= maxval;
					view->m_liftedEikonal.m_rotated_phasefield.m_distance[1] /= maxval;*/
					view->m_threadactivated = DISTANCE_MEAN_ITERATION;
#ifdef DEBUG_STATES
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph1.5_r_combined_distance.tif", view->m_liftedEikonal.m_rotated_phasefield.m_combined_distance);
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph1.5_r_distance0.tif", view->m_liftedEikonal.m_rotated_phasefield.m_distance[0]);
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph1.5_r_distance1.tif", view->m_liftedEikonal.m_rotated_phasefield.m_distance[1]);
#endif
					view->Invalidate(FALSE);
				}
			}

			if (++cyc > 7) { view->Invalidate(FALSE); cyc = 0; }

			

		}
		else if (view->m_threadactivated == DISTANCE_MEAN_ITERATION) {
			view->m_liftedEikonal.m_inicountourCalculator.GetDistancemean(view->m_liftedEikonal.m_rotated_phasefield.m_combined_distance, view->m_liftedEikonal.m_rotated_phasefield.m_plane_slice);
#ifdef DEBUG_STATES
			save_image("Y:/BIOMAG/shortest path/interm_imgs/ph2_r_combined_distance.tif", view->m_liftedEikonal.m_rotated_phasefield.m_combined_distance);
			save_image("Y:/BIOMAG/shortest path/interm_imgs/ph2_r_distance_0.tif", view->m_liftedEikonal.m_rotated_phasefield.m_distance[0]);
			save_image("Y:/BIOMAG/shortest path/interm_imgs/ph2_r_distance_1.tif", view->m_liftedEikonal.m_rotated_phasefield.m_distance[1]);
			save_image("Y:/BIOMAG/shortest path/interm_imgs/ph2_r_data.tif", view->m_liftedEikonal.m_rotated_phasefield.m_data);
			save_image("Y:/BIOMAG/shortest path/interm_imgs/ph2_r_flow_idx.tif", view->m_liftedEikonal.m_rotated_phasefield.m_flow_idx);

			Sleep(300);
			save_work_img("Y:/BIOMAG/shortest path/interm_imgs/ph2_slice_gx.tif",view->m_liftedEikonal.m_inicountourCalculator.m_gx2D);
			save_work_img("Y:/BIOMAG/shortest path/interm_imgs/ph2_slice_gy.tif", view->m_liftedEikonal.m_inicountourCalculator.m_gy2D);
			save_work_img("Y:/BIOMAG/shortest path/interm_imgs/ph2_slice_distance0.tif", view->m_liftedEikonal.m_inicountourCalculator.m_distance2D[0]);
			save_work_img("Y:/BIOMAG/shortest path/interm_imgs/ph2_slice_distance1.tif", view->m_liftedEikonal.m_inicountourCalculator.m_distance2D[1]);
#endif
			if (view->m_liftedEikonal.m_inicountourCalculator.m_bInited)
				view->m_threadactivated = PLANE_PHASEFIELD_ITERATION;
		}
		else if (view->m_threadactivated == PLANE_PHASEFIELD_ITERATION) { // plane
			if (view->m_liftedEikonal.m_inicountourCalculator.m_nect) {
				view->m_liftedEikonal.m_inicountourCalculator.Iterate();
				
			}
			else {
				//view->m_imageOp.GetXTestBound(view->m_xsee, view->m_liftedEikonal.m_boundcontour);
				vector<unsigned> size = view->m_liftedEikonal.m_rotated_phasefield.m_distance[0].GetSize();
				int xs(size[0]), ys(size[1]), zs(size[2]);
				unordered_set<unsigned long>& bound = view->m_liftedEikonal.m_inicountourCalculator.RetrieveBound();
				view->m_liftedEikonal.m_boundcontour.clear();
				for (auto& it = bound.begin(); it != bound.end(); it++) {
					int z = *it >> 16;
					int y = *it & 0b1111111111111111;
					vector<double> rotated;
					rotated = view->m_liftedEikonal.m_rotated_phasefield.m_sample_image.TransformContinuousIndexToPhysicalPoint({ view->m_liftedEikonal.m_rotated_phasefield.m_plane_slice, (double)y, (double)z });
					rotated = view->m_liftedEikonal.m_phasefield.m_sample_image.TransformPhysicalPointToContinuousIndex(rotated);
					//rotate(view->m_liftedEikonal.m_phasefield.rotation_matrix, { view->m_liftedEikonal.m_rotated_phasefield.m_plane_slice, (double)y, (double)z }, rotated);
					view->m_liftedEikonal.m_boundcontour.push_back(point_to_representation(round(rotated[0]), round(rotated[1]), round(rotated[2])));
				}
				view->m_imageOp.GetPlaneDistMap(ys, zs, bound);
				if (xs > 0) {
					sitk::Image& passi = view->m_imageOp.GetIniMap(xs, ys, zs, view->m_liftedEikonal.m_rotated_phasefield.m_plane_slice);
					realnum maxdist = view->m_liftedEikonal.m_rotated_phasefield.m_currentdistance;
					view->m_transport.TrInit(view->m_liftedEikonal.m_rotated_phasefield.m_combined_distance, passi, maxdist);
#ifdef DEBUG_STATES
					save_work_img("Y:/BIOMAG/shortest path/interm_imgs/ph3_slice_distance.tif", view->m_imageOp.m_loc);
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph3_transport_gx.tif", view->m_transport.m_gx);
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph3_transport_gy.tif", view->m_transport.m_gy);
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph3_transport_gz.tif", view->m_transport.m_gz);
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph3_transport_init0.tif", view->m_transport.m_transportfunction[0]);
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph3_transport_init1.tif", view->m_transport.m_transportfunction[1]);
					save_image("Y:/BIOMAG/shortest path/interm_imgs/ph3_transport_bound.tif", view->m_transport.m_isboundary);
					save_work_img("Y:/BIOMAG/shortest path/interm_imgs/ph3_slice_field.tif", view->m_liftedEikonal.m_inicountourCalculator.m_field2D);
#endif

				}
				view->m_threadactivated = TRANSPORT_FUNCTION_ITERATION;

			}
			if (++cyc > 7) { view->Invalidate(FALSE); cyc = 0; }
		}
		else if (view->m_threadactivated == TRANSPORT_FUNCTION_ITERATION) {
			view->m_transport.TrControl(100000);
			vector<unsigned> size = view->m_liftedEikonal.m_phasefield.m_data.GetSize();
			unsigned int xs(size[0]), ys(size[1]), zs(size[2]);
			vector<unsigned int> sample_size = { xs, ys, zs };
#ifdef DEBUG_STATES
			save_image("Y:/BIOMAG/shortest path/interm_imgs/ph4_transport0.tif", view->m_transport.m_transportfunction[0]);
			save_image("Y:/BIOMAG/shortest path/interm_imgs/ph4_bound.tif", view->m_transport.m_isboundary);
#endif
			resample_img(view->m_transport.m_transportfunction[0], view->m_transport.m_transportfunction[0], view->m_liftedEikonal.m_phasefield.rotation_matrix, sample_size, true);
			//neg_to_minus1(view->m_transport.m_transportfunction[0]);
			resample_img<sitk::sitkNearestNeighbor>(view->m_transport.m_isboundary, view->m_transport.m_isboundary, view->m_liftedEikonal.m_phasefield.rotation_matrix, sample_size, true);
#ifdef DEBUG_STATES
			save_image("Y:/BIOMAG/shortest path/interm_imgs/ph4_unrot_transportfn.tif", view->m_transport.m_transportfunction[0]);
			save_image("Y:/BIOMAG/shortest path/interm_imgs/ph4_unrot_bound.tif", view->m_transport.m_isboundary);
#endif
			view->m_transport.GetDispSlice(Talox, view->m_xsee, view->m_dispd2); // TaloX - enum
			view->m_transport.GetDispSlice(Taloy, view->m_ysee, view->m_dispd1);
			view->m_transport.GetDispSlice(Taloz, view->m_zsee, view->m_disp);
			view->Invalidate();
			//view->m_prevthreadactivated = 2;
			Sleep(300);
			view->m_threadactivated = DONE_ITERATION;
		}
		else if (view->m_threadactivated == DONE_ITERATION) {
			if (!dumped) {
				double length = omp_get_wtime() - start;
				f_stream << "execution time: " << length << endl;
				f_stream.close();
				_DUMP_PROFILE_INFO("Y:/BIOMAG/shortest path/profiling.txt");
				dumped = true;
			}
			Sleep(300);
		}

	}

	return 0;
}

void CChildView::InitThread()
{
	if (m_threadactivated) return;
	if (!m_bispoint) {
		m_start_point.x = 40;
		m_end_point.x = 201;
		m_start_point.y = m_end_point.y = 87;
		m_start_point.z = m_end_point.z = 54;
		m_bispoint = 2;
	}
	
	
	if (m_bispoint <= 1) {
		m_end_point.x = -100;
		m_end_point.y = -100;
	}
	m_imageOp.CreateTestInput(m_expfac);
	sitk::Image& data = m_imageOp.GetTestInput();
	vector<uint32_t> size = data.GetSize();
	if (m_bispoint == 1) {
		m_end_point.z = m_xsee;
	}
	//sitk::Image data_im = vox_img_2_sitk(data);
	m_liftedEikonal.m_phasefield.Initialize(data, m_start_point, m_end_point);


	CWinThread* thread = AfxBeginThread(BackgroundThread,this,THREAD_PRIORITY_NORMAL,0,CREATE_SUSPENDED,0);
	if(thread == 0) return;
	m_threadactivated = DISTANCE_ITERATION;
	thread->m_bAutoDelete = true;
	thread->ResumeThread();
	Invalidate();

}

void CChildView::PauseThread(int threadstat) // 1-(re-)start 2-suspended
{
	if (!m_threadactivated) return;
	if (m_threadactivated < DONE_ITERATION) {
		m_prevthreadactivated = m_threadactivated;
		if (m_prevthreadactivated == 3) m_prevthreadactivated = 2;
		m_threadactivated = DONE_ITERATION;
	}
	else {
		m_threadactivated = m_prevthreadactivated;
	}
	Invalidate();
}

void CChildView::StopThread()
{
	m_threadactivated = 0;
	_DUMP_PROFILE_INFO("Y:/BIOMAG/shortest path/profiling.txt")
	Invalidate();
}

void CChildView::GetAllPathX()
{
	std::vector<POINT3D>& bound = m_liftedEikonal.m_boundcontour;
	int si = bound.size();
	if (si > MAXMINPATH) si = MAXMINPATH;
	for (int jj = 0; jj < MAXMINPATH; ++jj) {
		m_liftedEikonal.m_minpath[0][jj].clear();
		//m_liftedEikonal.m_minpath[1][jj].clear();
	}
	int& mj = m_liftedEikonal.m_j[0];
	mj = 0;
	for (int jj = 0; jj < si; ++jj) {
		auto [x, y, z] = representation_to_point<int>(bound[jj]);
		m_liftedEikonal.ResolvePath(x + 1, y, z, true, 0, mj);
		m_liftedEikonal.ResolvePath(x - 1, y, z, true, 0, mj);
		if (m_liftedEikonal.m_minpath[0][mj].size()) ++mj;

		Invalidate();
	}


}

void CChildView::OnLButtonDown(UINT nFlags, CPoint point) {
	pressed = true;
	CWnd::OnLButtonDown(nFlags, point);
}


void CChildView::OnLButtonUp(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	if (!pressed) return;

		int8_t* x_path_buffer = m_path_image_x.GetBufferAsInt8();
		int8_t* y_path_buffer = m_path_image_y.GetBufferAsInt8();
		int8_t* z_path_buffer = m_path_image_z.GetBufferAsInt8();
		vector<unsigned> size = m_path_image_x.GetSize();
		memset(x_path_buffer, -100, size[0] * size[1] * size[2]);
		memset(y_path_buffer, -100, size[0] * size[1] * size[2]);
		memset(z_path_buffer, -100, size[0] * size[1] * size[2]);
	
	if (m_threadactivated == DONE_ITERATION) {
		if (point.x < m_disp.xs) {
			m_liftedEikonal.ResolvePath(point.x, point.y, m_zsee, true, 0, m_liftedEikonal.m_j[0]);
			if (m_liftedEikonal.m_minpath[0][m_liftedEikonal.m_j[0]].size()) ++m_liftedEikonal.m_j[0];
			m_liftedEikonal.ResolvePath(point.x, point.y, m_zsee, true, 1, m_liftedEikonal.m_j[1]);
			if (m_liftedEikonal.m_minpath[1][m_liftedEikonal.m_j[1]].size()) ++m_liftedEikonal.m_j[1];
		}
		else if (point.x < 2 * m_disp.xs + 1) {
			m_liftedEikonal.ResolvePath(point.x - m_disp.xs - 1, m_ysee, point.y, true, 0, m_liftedEikonal.m_j[0]);
			if (m_liftedEikonal.m_minpath[0][m_liftedEikonal.m_j[0]].size()) ++m_liftedEikonal.m_j[0];
			m_liftedEikonal.ResolvePath(point.x - m_disp.xs - 1, m_ysee, point.y, true, 1, m_liftedEikonal.m_j[1]);
			if (m_liftedEikonal.m_minpath[1][m_liftedEikonal.m_j[1]].size()) ++m_liftedEikonal.m_j[1];
		}
		else if(point.x < m_disp.xs+m_dispd1.xs+m_dispd2.xs) {
			m_liftedEikonal.ResolvePath(m_xsee, point.x - 2 * (m_disp.xs + 1), point.y, true, 0, m_liftedEikonal.m_j[0]);
			if (m_liftedEikonal.m_minpath[0][m_liftedEikonal.m_j[0]].size()) ++m_liftedEikonal.m_j[0];
			m_liftedEikonal.ResolvePath(m_xsee, point.x - 2 * (m_disp.xs + 1), point.y, true, 1, m_liftedEikonal.m_j[1]);
			if (m_liftedEikonal.m_minpath[1][m_liftedEikonal.m_j[1]].size()) ++m_liftedEikonal.m_j[1];
		}
		Invalidate();
	}
	else if (!m_threadactivated) {
		vector<uint32_t> size = m_imageOp.m_testimage.GetSize();
		int xs(size[0]), ys(size[1]);
		if (m_disp.xs && point.x < xs && point.y < ys) {
			if (!m_bispoint) {
				m_start_point = CVec3(point.x, point.y, m_zsee);
				++m_bispoint;
			}
			else {
				m_end_point = CVec3(point.x, point.y, m_zsee);
				++m_bispoint;
			}
		}
		Invalidate();
	}

	CWnd::OnLButtonUp(nFlags, point);
	pressed = false;
}

void CChildView::OnRButtonUp(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default
	//m_bispoint = 0;
	Invalidate();

	CWnd::OnRButtonUp(nFlags, point);
}



// C:\Users\jmolnar\Documents\Visual Studio 2008\Projects\MinArea\MinArea\ChildView.cpp : implementation file
//

#include "stdafx.h"
#include "MinArea.h"
#include "ChildView.h"


// CControlDlg dialog

IMPLEMENT_DYNAMIC(CControlDlg, CDialog)

CControlDlg::CControlDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CControlDlg::IDD, pParent)
{
	m_pView = 0;
	m_startstate = 0;
	m_flowray = 3;
	m_bini = false;
	m_bforced = true;
}

CControlDlg::~CControlDlg()
{
}

void CControlDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_BUTTON4, m_startstop);
	DDX_Control(pDX, IDC_EDIT1, m_cdist);
	DDX_Control(pDX, IDC_SLIDER1, m_zlevel);
	DDX_Control(pDX, IDC_EDIT2, m_ezlevel);
	DDX_Control(pDX, IDC_SLIDER2, m_cstartdir);
	DDX_Control(pDX, IDC_SLIDER3, m_cstopdir);
	DDX_Control(pDX, IDC_EDIT3, m_cflowray);
	DDX_Control(pDX, IDC_SLIDER4, m_cdatafac);
	DDX_Control(pDX, IDC_CHECK2, m_ccursee);
	DDX_Control(pDX, IDC_BUTTON8, m_cmsw);
}


BEGIN_MESSAGE_MAP(CControlDlg, CDialog)
	ON_BN_CLICKED(IDC_BUTTON1, &CControlDlg::OnBnClickedOpen)
	ON_BN_CLICKED(IDC_BUTTON2, &CControlDlg::OnBnClickedIntImage1)
	ON_BN_CLICKED(IDC_BUTTON3, &CControlDlg::OnBnClickedGetBC)
	ON_BN_CLICKED(IDC_BUTTON4, &CControlDlg::OnBnClickedStartStop)
	ON_BN_CLICKED(IDC_BUTTON5, &CControlDlg::OnBnClickedSuspend)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER1, &CControlDlg::OnNMCustomdrawSlider1)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER2, &CControlDlg::OnNMCustomdrawSlider2)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER3, &CControlDlg::OnNMCustomdrawSlider3)
	ON_EN_CHANGE(IDC_EDIT3, &CControlDlg::OnEnChangeEdit3)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER4, &CControlDlg::OnNMCustomdrawSlider4)
	ON_EN_CHANGE(IDC_EDIT2, &CControlDlg::OnEnChangeEdit2)
	ON_BN_CLICKED(IDC_CHECK1, &CControlDlg::OnBnClickedCheckSee)
	ON_EN_CHANGE(IDC_EDIT1, &CControlDlg::OnEnChangeEdit1)
	ON_BN_CLICKED(IDC_BUTTON6, &CControlDlg::OnBnClickedMinPath)
	ON_BN_CLICKED(IDC_CHECK2, &CControlDlg::OnBnClickedCSee)
	ON_BN_CLICKED(IDC_BUTTON7, &CControlDlg::OnBnClickedIntImage2)
	ON_BN_CLICKED(IDC_BUTTON8, &CControlDlg::OnBnClickedWWOCorr)
	ON_BN_CLICKED(IDC_CHECK3, &CControlDlg::OnBnClickedSwitchView)
	ON_BN_CLICKED(IDC_BUTTON10, &CControlDlg::OnBnClickedTestBuild)
END_MESSAGE_MAP()


// CControlDlg message handlers

void CControlDlg::OnBnClickedOpen()
{
	// TODO: Add your control notification handler code here

	CFileDialog cfd(TRUE);
	cfd.m_ofn.Flags |= OFN_ALLOWMULTISELECT;
	char strFiles[(MAX_PATH+1)*10+1];
	strFiles[0] = 0;
	cfd.m_ofn.lpstrFile = strFiles;
	cfd.m_ofn.nMaxFile  = (MAX_PATH+1)*10+1;

	if (cfd.DoModal() == IDOK) {
		POSITION pos = cfd.GetStartPosition();
		if(!pos) return;
		bool b1st = true;
		CString cfn;
		int ser = -1;
		do {
			cfn = cfd.GetNextPathName(pos); 
			sitk::Image ci = sitk::ReadImage(string((LPCTSTR)cfn), sitk::sitkUnknown);
			auto type = ci.GetPixelID();
			ci = sitk::Cast(ci, sitk::sitkFloat64);
			m_pView->m_imageOp.m_testimage = sitk::RescaleIntensity(ci, 0, 1);
			vector<uint32_t> size = m_pView->m_imageOp.m_testimage.GetSize();
			int xs = size[0], ys = size[1], zs = size[2];
			//sitk_2_vox_img(ci, m_pView->m_imageOp.m_testimage);
			m_pView->m_imageOp.GauTest(false);
			m_pView->m_imageOp.GauTest(true);
			m_pView->m_imageOp.GauTest(false);
			m_pView->m_imageOp.GauTest(true);
			m_pView->m_imageOp.CreateTestInput(m_pView->m_expfac);
			m_pView->m_path_image_x = sitk::Image(size, sitk::sitkInt8)-100;
			m_pView->m_path_image_y = sitk::Image(size, sitk::sitkInt8)-100;
			m_pView->m_path_image_z = sitk::Image(size, sitk::sitkInt8)-100;
			//CImage ci; ci.Load(LPCTSTR(cfn));

			if (!xs || !ys || !zs) return;
			double* im_buffer = m_pView->m_imageOp.m_testimage.GetBufferAsDouble();
			m_pView->m_work.Set(xs, ys);
			for (int yy = 0; yy < ys; yy++) {
				for (int xx = 0; xx < xs; xx++) {
					m_pView->m_work[yy][xx] = BUF_IDX(im_buffer, xs, ys, zs, xx, yy, m_pView->m_zsee);
				}
			}
			//m_pView->m_work = m_pView->m_imageOp.m_testimage[m_pView->m_zsee];
			m_pView->m_valid = 1;
			m_pView->m_grays = true;
			m_pView->m_color = false;
			m_zlevel.SetRange(0, zs - 1, 1);
			m_cstartdir.SetRange(0, ys - 1, 1);
			m_cstopdir.SetRange(0, xs - 1, 1);
			m_bini = true;

			m_cdatafac.SetRange(1, EXPCOEF_RANGE, 1);
			m_cdatafac.SetPos(EXP_COEF_DEF);
			char txt[22];
			sprintf_s(txt, 20, "%d", m_flowray);
			m_cflowray.SetWindowTextA(txt);
			m_pView->m_work.GetDispImg(m_pView->m_disp);
			OnNMCustomdrawSlider2(0, 0);
			OnNMCustomdrawSlider3(0, 0);
			m_pView->Invalidate();

		} while(pos);
		// some init ok here

	}


}

void CControlDlg::OnBnClickedIntImage1() // synth. image 1
{
	// TODO: Add your control notification handler code here
	if (m_bini) return;
	m_pView->m_imageOp.CreateTestImage(XS_,YS_,ZS_);
	{
		vector<uint32_t> size = m_pView->m_imageOp.m_testimage.GetSize();
		m_pView->m_path_image_x = sitk::Image(size, sitk::sitkInt8) - 100;
		m_pView->m_path_image_y = sitk::Image(size, sitk::sitkInt8) - 100;
		m_pView->m_path_image_z = sitk::Image(size, sitk::sitkInt8) - 100;
		int xs = size[0], ys = size[1], zs = size[2]; 
		double* im_buffer = m_pView->m_imageOp.m_testimage.GetBufferAsDouble();
		m_pView->m_work.Set(xs, ys);
		for (int yy = 0; yy < ys; yy++) {
			for (int xx = 0; xx < xs; xx++) {
				m_pView->m_work[yy][xx] = BUF_IDX(im_buffer, xs, ys, zs, xx, yy, m_pView->m_zsee);
			}
		}
		m_pView->m_valid = 1;
		m_pView->m_grays = true;
		m_pView->m_color = false;
	}
	{
		m_zlevel.SetRange(0,ZS_-1,1); // from z
		m_cstartdir.SetRange(0,YS_-1,1); // from y
		m_cstopdir.SetRange(0,XS_-1,1); // from x
		m_bini = true;

		m_cdatafac.SetRange(1, EXPCOEF_RANGE,1); m_cdatafac.SetPos(EXP_COEF_DEF);
		char txt[22]; sprintf_s(txt,20,"%d",m_flowray);	m_cflowray.SetWindowTextA(txt);
	}
	m_pView->m_work.GetDispImg(m_pView->m_disp);
	OnNMCustomdrawSlider2(0,0);
	OnNMCustomdrawSlider3(0,0);
	m_pView->Invalidate();

}

void CControlDlg::OnBnClickedGetBC() // contour on plane
{
	// TODO: Add your control notification handler code here
	//m_pView->m_imageOp.GetXTestBound(m_pView->m_xsee/* - 1*/, m_pView->m_liftedEikonal.m_boundcontour);
	//if (m_pView->m_liftedEikonal.m_phasefield.m_distance[0].xs > 0) {
	//	int xs(m_pView->m_liftedEikonal.m_phasefield.m_distance[0].xs), ys(m_pView->m_liftedEikonal.m_phasefield.m_distance[0].ys), zs(m_pView->m_liftedEikonal.m_phasefield.m_distance[0].zs);
	//	SVoxImg<SWorkImg<realnum>> & passi = m_pView->m_imageOp.GetIniMap(xs, ys, zs, m_pView->m_xsee/* - 1*/);
	//	//m_pView->m_transport.TrInit(m_pView->m_liftedEikonal.m_phasefield.m_smoothdist[0],passi, m_pView->m_liftedEikonal.m_currentdistance[0]);
	//	m_pView->m_transport.TrInit(m_pView->m_liftedEikonal.m_phasefield.m_distance[0], passi, m_pView->m_liftedEikonal.m_phasefield.m_currentdistance);

		//--m_pView->m_xsee;
		/*
		m_pView->m_transport.GetDispSlice(Talox, m_pView->m_xsee, m_pView->m_dispd2); // TaloX - enum
		m_pView->m_transport.GetDispSlice(Taloy, m_pView->m_ysee, m_pView->m_dispd1);
		m_pView->m_transport.GetDispSlice(Taloz, m_pView->m_zsee, m_pView->m_disp);*/
	//}

	//m_pView->Invalidate();
}

void CControlDlg::OnBnClickedStartStop()
{
	// TODO: Add your control notification handler code here
	if (!m_startstate) {
		m_pView->InitThread();
		m_startstate = 1;
		m_startstop.SetWindowTextA("Stop");
	}
	else {
		m_pView->StopThread();
		m_startstate = 0;
		m_startstop.SetWindowTextA("Start");
	}
}


void CControlDlg::OnBnClickedSuspend() // pause/continue
{
	// TODO: Add your control notification handler code here

	if (!m_startstate) return;

	if (m_startstate == 1) {
		m_startstate = 2; // pause
	}
	else {
		m_startstate = 1;
	}

	m_pView->PauseThread(m_startstate);
}

////////////////////////////////////////////////////////
// view from z - y - x
////////////////////////////////////////////////////////

void CControlDlg::OnNMCustomdrawSlider1(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	// TODO: Add your control notification handler code here
	vector<uint32_t> size = m_pView->m_imageOp.m_testimage.GetSize();
	if (!m_bini) goto ret;
		m_pView->m_zsee = m_zlevel.GetPos();
		char txt[22];
		sprintf_s(txt,20,"%d",m_pView->m_zsee);
		m_ezlevel.SetWindowTextA(txt);

		int xs = size[0], ys = size[1], zs = size[2];
		double* im_buffer = m_pView->m_imageOp.m_testimage.GetBufferAsDouble();
		m_pView->m_work.Set(xs, ys);
		for (int yy = 0; yy < ys; yy++) {
			for (int xx = 0; xx < xs; xx++) {
				m_pView->m_work[yy][xx] = BUF_IDX(im_buffer, xs, ys, zs, xx, yy, m_pView->m_zsee);
			}
		}

	if (!m_pView->m_btransportview) {
		m_pView->m_work.GetDispImg(m_pView->m_disp);
	}
	else {
		m_pView->m_transport.GetDispSlice(Taloz, m_pView->m_zsee, m_pView->m_disp);
	}
		
	m_pView->Invalidate();
ret:
	*pResult = 0;
}


void CControlDlg::OnNMCustomdrawSlider2(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	// TODO: Add your control notification handler code here
	if (!m_bini) goto ret;
	sitk::Image &testimg = m_pView->m_imageOp.m_testimage;
	double* testimg_buffer = testimg.GetBufferAsDouble();

	if (!m_pView->m_btransportview) {
		vector<uint32_t> size = m_pView->m_imageOp.m_testimage.GetSize();
		int xs = size[0], ys = size[1], zs = size[2];
		SWorkImg<realnum>& yslice = m_pView->m_intey;
		yslice.Set(xs, zs);
		int ylev = m_cstartdir.GetPos();

		for (int zz = 0; zz < zs; ++zz) {
			for (int xx = 0; xx < xs; ++xx) {
				yslice[zz][xx] = BUF_IDX(testimg_buffer, xs, ys, zs, xx, ylev, zz);
			}
		}

		yslice.GetDispImg(m_pView->m_dispd1);
		m_pView->m_ysee = m_cstartdir.GetPos();
	}
	else {
		m_pView->m_ysee = m_cstartdir.GetPos();
		m_pView->m_transport.GetDispSlice(Taloy, m_pView->m_ysee, m_pView->m_dispd1);
	}

	m_pView->Invalidate();
ret:
	if (pResult) *pResult = 0;
}

void CControlDlg::OnNMCustomdrawSlider3(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	// TODO: Add your control notification handler code here
	if (!m_bini) goto ret;


	sitk::Image &testimg = m_pView->m_imageOp.m_testimage;
	double* testimg_buffer = testimg.GetBufferAsDouble();
	if (!m_pView->m_btransportview) {
		vector<uint32_t> size = m_pView->m_imageOp.m_testimage.GetSize();
		int xs = size[0], ys = size[1], zs = size[2]; 
		SWorkImg<realnum>& xslice = m_pView->m_intex;
		xslice.Set(ys, zs);
		int xlev = m_cstopdir.GetPos();

		for (int zz = 0; zz < zs; ++zz) {
			for (int yy = 0; yy < ys; ++yy) {
				xslice[zz][yy] = BUF_IDX(testimg_buffer, xs, ys, zs, xlev, yy, zz);
			}
		}

		xslice.GetDispImg(m_pView->m_dispd2);
		m_pView->m_xsee = m_cstopdir.GetPos();
	}
	else {
		m_pView->m_xsee = m_cstopdir.GetPos();
		m_pView->m_transport.GetDispSlice(Talox, m_pView->m_xsee, m_pView->m_dispd2);
	}


	m_pView->Invalidate();
ret:
	if (pResult) *pResult = 0;
}

void CControlDlg::OnEnChangeEdit3()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
	char txt[22];
	m_cflowray.GetWindowTextA(txt,20);
	m_flowray = atoi(txt);
	m_bforced = true;
	if (m_flowray < 1) {
		m_flowray = 1;
		m_cflowray.SetWindowTextA("1");
	}
	//if (m_bini && !m_startstate) OnBnClickedIntImage1();
}

void CControlDlg::OnNMCustomdrawSlider4(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	// TODO: Add your control notification handler code here
	if (m_bini) {
		m_pView->m_expfac = m_cdatafac.GetPos();
		m_bforced = true;
	}

	*pResult = 0;
}

void CControlDlg::OnEnChangeEdit2()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
}

void CControlDlg::OnBnClickedCheckSee()
{
	// TODO: Add your control notification handler code here
	m_pView->m_bsee = !m_pView->m_bsee;
	m_pView->Invalidate();
}

void CControlDlg::OnEnChangeEdit1()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialog::OnInitDialog()
	// function and call CRichEditCtrl().SetEventMask()
	// with the ENM_CHANGE flag ORed into the mask.

	// TODO:  Add your control notification handler code here
}

void CControlDlg::OnBnClickedMinPath() // get minimal paths
{
	// TODO: Add your control notification handler code here
	m_pView->GetAllPathX();
}

void CControlDlg::OnBnClickedCSee()
{
	// TODO: Add your control notification handler code here
	m_pView->m_curvsee = m_ccursee.GetCheck();
	m_pView->Invalidate();
}

void CControlDlg::OnBnClickedIntImage2() // synth. image 2
{
	// TODO: Add your control notification handler code here
	if (m_bini) return;
	m_pView->m_imageOp.CreateTestImage2(XS_,YS_,ZS_);
	{
		vector<uint32_t> size = m_pView->m_imageOp.m_testimage.GetSize();
		int xs = size[0], ys = size[1], zs = size[2];
		double* im_buffer = m_pView->m_imageOp.m_testimage.GetBufferAsDouble();
		m_pView->m_work.Set(xs, ys);
		for (int yy = 0; yy < ys; yy++) {
			for (int xx = 0; xx < xs; xx++) {
				m_pView->m_work[yy][xx] = BUF_IDX(im_buffer, xs, ys, zs, xx, yy, m_pView->m_zsee);
			}
		}
		m_pView->m_valid = 1;
		m_pView->m_grays = true;
		m_pView->m_color = false;
	}
	{
		m_zlevel.SetRange(0,ZS_-1,1); // from z
		m_cstartdir.SetRange(0,YS_-1,1); // from y
		m_cstopdir.SetRange(0,XS_-1,1); // from x
		m_bini = true;

		m_cdatafac.SetRange(1, EXPCOEF_RANGE,1); m_cdatafac.SetPos(EXP_COEF_DEF);
		char txt[22]; sprintf_s(txt,20,"%d",m_flowray);	m_cflowray.SetWindowTextA(txt);
	}
	m_pView->m_work.GetDispImg(m_pView->m_disp);
	OnNMCustomdrawSlider2(0,0);
	OnNMCustomdrawSlider3(0,0);
	m_pView->Invalidate();
}

void CControlDlg::OnBnClickedWWOCorr() // correction using mean surf metric
{
	// TODO: Add your control notification handler code here
	extern bool g_modeswitch;
	g_modeswitch = !g_modeswitch;
	if (g_modeswitch)
		m_cmsw.SetWindowTextA("With corr");
	else
		m_cmsw.SetWindowTextA("W/O corr");
}


void CControlDlg::OnBnClickedSwitchView()
{
	// TODO: Add your control notification handler code here
	m_pView->m_btransportview = !m_pView->m_btransportview;
	if (m_pView->m_btransportview) {

		m_pView->m_transport.GetDispSlice(Talox, m_pView->m_xsee, m_pView->m_dispd2);
		m_pView->m_transport.GetDispSlice(Taloy, m_pView->m_ysee, m_pView->m_dispd1);
		m_pView->m_transport.GetDispSlice(Taloz, m_pView->m_zsee, m_pView->m_disp);

	}

	m_pView->Invalidate();
}


void CControlDlg::OnBnClickedTestBuild()
{
	// TODO: Add your control notification handler code here
	m_pView->m_transport.TrControl(333);
	if (m_pView->m_btransportview) {

		m_pView->m_transport.GetDispSlice(Talox, m_pView->m_xsee, m_pView->m_dispd2);
		m_pView->m_transport.GetDispSlice(Taloy, m_pView->m_ysee, m_pView->m_dispd1);
		m_pView->m_transport.GetDispSlice(Taloz, m_pView->m_zsee, m_pView->m_disp);

	}
	m_pView->Invalidate();
}
