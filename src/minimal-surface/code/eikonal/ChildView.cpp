
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
#include "FMSigned.h"
#include <thread>
#ifdef DEBUG_STATES
#include "Debug.h"
#endif

#define EXPCOEF_RANGE 100
using namespace std;
namespace sitk = itk::simple;
// CChildView

CChildView::CChildView()
{
	m_pControl = 0;
	mAlgorithmStage = 0;
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
			std::vector<POINT3D> &bound = mBoundaryContour;
			int si = (int)bound.size();
			if (si) {
				for (int ii = 0; ii < si; ++ii) {
					auto [x, y, z] = representation_to_point<int>(bound[ii]);
					//if(m_xsee == x)
						dc.SetPixelV(offx+y,z,0xffff);
				}
			}
		}
		offx += m_dispd2.xs+1;
		
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

	vector<unsigned> size = mInputImage.GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);

	if (mAlgorithmStage)  {
		for (int ii = 0; ii < 2; ++ii) {
			if (mPhaseFieldMap[ii] == NULL) return;
			const double* field_buffer = mPhaseFieldMap[ii]->GetBufferAsDouble();
			const int* meeting_points_buffer = mMeetingPointsMap->GetBufferAsInt32();
			int zsee = 0;
			if (m_pControl) {
				zsee = m_zsee;
			}

			if (mAlgorithmStage == DISTANCE_ITERATION || mAlgorithmStage == DONE_ITERATION) {
				if (m_disp.xs && !m_btransportview) {
					for (int yy = 0; yy < ys; ++yy) {
						for (int xx = 0; xx < xs; ++xx) {
							Vec3<double> pixel(xx, yy, zsee);
							if (abs((mPlaneNormal*pixel).Sum() + mPlaneOffset) < 0.5) {
								dc.SetPixelV(xx, yy, 0xff00ff);
							}
							else if (meeting_points_buffer != NULL && meeting_points_buffer[xx + xs * (yy + ys * zsee)] > 0) {
								dc.SetPixelV(xx, yy, 0xff0000);
							}
							else if (field_buffer != NULL && field_buffer[xx + xs * (yy + ys * zsee)] > -0.9f && field_buffer[xx + xs * (yy + ys * zsee)] < 0.3f)
								dc.SetPixelV(xx, yy, !ii ? 0xff : 0xff00);
						}
					}
				}

				if (m_dispd1.xs) 
					for (int zz = 0; zz < zs; ++zz) {
						for (int xx = 0; xx < xs; ++xx) {
							Vec3<double> pixel(xx, m_ysee, zz);
							if (abs((mPlaneNormal * pixel).Sum() + mPlaneOffset) < 0.5) {
								dc.SetPixelV(m_disp.xs + 1 + xx, zz, 0xff00ff);
							}
							else if (meeting_points_buffer != NULL && meeting_points_buffer[xx + xs * (m_ysee + ys * zz)] > 0) {
								dc.SetPixelV(m_disp.xs + 1 + xx, zz, 0xff0000);
							}
							else if (field_buffer != NULL && field_buffer[xx + xs * (m_ysee + ys * zz)] > -0.9f && field_buffer[xx + xs * (m_ysee + ys * zz)] < 0.3f)
								dc.SetPixelV(m_disp.xs + 1 + xx, zz, !ii ? 0xff : 0xff00);
						
						}
					}
				if (m_dispd2.xs) { // xslice
					for (int zz = 0; zz < zs; ++zz) {
						for (int yy = 0; yy < ys; ++yy) {
							Vec3<double> pixel(m_xsee, yy, zz);
							if (abs((mPlaneNormal * pixel).Sum() + mPlaneOffset) < 0.5) {
								dc.SetPixelV(2 * (m_disp.xs + 1) + yy, zz, 0xff00ff);
							}
							else if (meeting_points_buffer != NULL && meeting_points_buffer[m_xsee + xs * (yy + ys * zz)] > 0) {
								dc.SetPixelV(2 * (m_disp.xs + 1) + yy, zz, 0xff0000);
							}
							else if (field_buffer != NULL && field_buffer[m_xsee + xs * (yy + ys * zz)] > -0.9f && field_buffer[m_xsee + xs * (yy + ys * zz)] < 0.3f)
								dc.SetPixelV(2 * (m_disp.xs + 1) + yy, zz, !ii ? 0xff : 0xff00);
						}
					}
				}
			}

			
			else if (mPhaseField2DMap != NULL && mAlgorithmStage == PLANE_PHASEFIELD_ITERATION && m_dispd2.xs) { // xslice
				const double* field2d_buffer = mEstimator.GetInitialContourCalculator().GetPhaseField().GetBufferAsDouble();
				vector<unsigned> size = mEstimator.GetInitialContourCalculator().GetPhaseField().GetSize();
				int ph_xs(size[0]), ph_ys(size[1]);
					for (int zz = 0; zz < zs; ++zz) {
						for (int yy = 0; yy < ys; ++yy) {
							if (BUF_IDX2D(field2d_buffer, ph_xs, ph_ys, yy, zz) > 0.5f)
								dc.SetPixelV(2 * (m_disp.xs + 1) + yy, zz, 0xff);
							else if (BUF_IDX2D(field2d_buffer, ph_xs, ph_ys, yy, zz) < -0.5f)
								dc.SetPixelV(2 * (m_disp.xs + 1) + yy, zz, 0xff0000);
						}
					}
			}
			if (!mEstimator.mMinimalPaths[ii].empty()) {
				CPen* oldpen = dc.SelectObject(&g_yellowpen);

				int8_t* path_x_buffer = mPathImageX.GetBufferAsInt8();
				int8_t* path_y_buffer = mPathImageY.GetBufferAsInt8();
				int8_t* path_z_buffer = mPathImageZ.GetBufferAsInt8();
				if (BUF_IDX3D(path_z_buffer, xs, ys, zs, 0, 0, m_zsee) == -100) {
					for (int yy = 0; yy < ys; yy++) {
						for (int xx = 0; xx < xs; xx++) {
							BUF_IDX3D(path_z_buffer, xs, ys, zs, xx, yy, m_zsee) = 0;
						}
					}
					for(int ii=0; ii<2; ii++)
						for (int jj = 0; jj < mEstimator.mMinimalPaths[ii].size(); ++jj) {
							std::vector<Vec3<int>>& minpath = mEstimator.mMinimalPaths[ii][jj];
							int point_count = (int)minpath.size();
							if (!point_count) break;
							for (int p = 0; p < point_count; ++p) {
								if(BUF_IDX3D(path_z_buffer, xs, ys, zs, (int)minpath[p].x(), (int)minpath[p].y(), m_zsee) < 1)
									BUF_IDX3D(path_z_buffer, xs, ys, zs, (int)minpath[p].x(), (int)minpath[p].y(), m_zsee) = minpath[p].z() < m_zsee ? -1 : 1;
							}
						}
				}
				for (int yy = 0; yy < ys; yy++) {
					for (int xx = 0; xx < xs; xx++) {
						int8_t val = BUF_IDX3D(path_z_buffer, xs, ys, zs, xx, yy, m_zsee);
						if (val) {
							int color = val < 0 ? 0x7fff : 0xffff;
							dc.SetPixelV(xx, yy, color);
						}
					}
				}

				if (BUF_IDX3D(path_y_buffer, xs, ys, zs, 0, m_ysee, 0) == -100) {
					for (int zz = 0; zz < zs; zz++) {
						for (int xx = 0; xx < xs; xx++) {
							BUF_IDX3D(path_y_buffer, xs, ys, zs, xx, m_ysee, zz) = 0;
						}
					}
					for (int ii = 0; ii < 2; ii++)
						for (int jj = 0; jj < mEstimator.mMinimalPaths[ii].size(); ++jj) {
							std::vector<Vec3<int>>& minpath = mEstimator.mMinimalPaths[ii][jj];
							int point_count = (int)minpath.size();
							if (!point_count) break;
							for (int p = 0; p < point_count; ++p) {
								if (BUF_IDX3D(path_y_buffer, xs, ys, zs, (int)minpath[p].x(), m_ysee, (int)minpath[p].z()) < 1)
									BUF_IDX3D(path_y_buffer, xs, ys, zs, (int)minpath[p].x(), m_ysee, (int)minpath[p].z()) = minpath[p].y() < m_ysee ? -1 : 1;
							}
						}
				}
				for (int zz = 0; zz < zs; zz++) {
					for (int xx = 0; xx < xs; xx++) {
						int8_t val = BUF_IDX3D(path_y_buffer, xs, ys, zs, xx, m_ysee, zz);
						if (val) {
							int color = val < 0 ? 0x7fff : 0xffff;
							dc.SetPixelV(xx + m_disp.xs + 1, zz, color);
						}
					}
				}

				if (BUF_IDX3D(path_x_buffer, xs, ys, zs, m_xsee, 0, 0) == -100) {
					for (int zz = 0; zz < zs; zz++) {
						for (int yy = 0; yy < ys; yy++) {
							BUF_IDX3D(path_x_buffer, xs, ys, zs, m_xsee, yy, zz) = 0;
						}
					}
					for (int ii = 0; ii < 2; ii++)
						for (int jj = 0; jj < mEstimator.mMinimalPaths[ii].size(); ++jj) {
							std::vector<Vec3<int>>& minpath = mEstimator.mMinimalPaths[ii][jj];
							int point_count = (int)minpath.size();
							if (!point_count) break;
							for (int p = 0; p < point_count; ++p) {
								if (BUF_IDX3D(path_x_buffer, xs, ys, zs, m_xsee, (int)minpath[p].y(), (int)minpath[p].z()) < 1)
									BUF_IDX3D(path_x_buffer, xs, ys, zs, m_xsee, (int)minpath[p].y(), (int)minpath[p].z()) = minpath[p].x() < m_xsee ? -1 : 1;
							}
						}
				}
				for (int zz = 0; zz < zs; zz++) {
					for (int yy = 0; yy < ys; yy++) {
						int8_t val = BUF_IDX3D(path_x_buffer, xs, ys, zs, m_xsee, yy, zz);
						if (val) {
							int color = val < 0 ? 0x7fff : 0xffff;
							dc.SetPixelV(yy + 2 * (m_disp.xs + 1), zz, color);
						}
					}
				}
				dc.SelectObject(oldpen);
			}
		}


		/*{
			const double* field_buffer = mEstimator.mAreaEikonal.GetPhaseFieldMap(0).GetBufferAsDouble();

			int xx(0), yy(0), zz(0);
			dc.MoveTo(xx,300-(int)(10*field_buffer[xx + xs * (m_ysee + ys * m_zsee)]));
			for (xx = 1; xx < xs; ++xx) dc.LineTo(xx*5,300-(int)(10*field_buffer[xx + xs * (m_ysee + ys * m_zsee)]));
			dc.MoveTo(yy,350-(int)(10*field_buffer[m_xsee + xs * (yy + ys * m_zsee)]));
			for (yy = 1; yy < ys; ++yy) dc.LineTo(yy*5,350-(int)(10*field_buffer[m_xsee + xs * (yy + ys * m_zsee)]));
			dc.MoveTo(zz,400-(int)(10*field_buffer[m_xsee + xs * (m_ysee + ys * zz)]));
			for (zz = 1; zz < zs; ++zz) dc.LineTo(zz*5,400-(int)(10*field_buffer[m_xsee + xs * (m_ysee + ys * zz)]));

		}*/

	}
			

	if (m_disp.xs) {
		if (m_bispoint > 0) {
			if (abs(m_zsee-mStartPoint.z()) < 3) dc.FillSolidRect((int)mStartPoint.x() -1, (int)mStartPoint.y() -1,3,3,0xff);
			if (abs(m_ysee- mStartPoint.y()) < 3) dc.FillSolidRect((int)mStartPoint.x() +m_disp.xs, (int)mStartPoint.z() -1,3,3,0xff);
			if (abs(m_xsee- mStartPoint.x()) < 3) dc.FillSolidRect((int)mStartPoint.y() +2*m_disp.xs+1, (int)mStartPoint.z() -1,3,3,0xff);
		}
		if (m_bispoint > 1) {
			if (abs(m_zsee- mEndPoint.z()) < 3) dc.FillSolidRect((int)mEndPoint.x() -1, (int)mEndPoint.y() -1,3,3,0xffff00);
			if (abs(m_ysee- mEndPoint.y()) < 3) dc.FillSolidRect((int)mEndPoint.x() +m_disp.xs, (int)mEndPoint.z() -1,3,3,0xffff00);
			if (abs(m_xsee- mEndPoint.x()) < 3) dc.FillSolidRect((int)mEndPoint.y() +2*m_disp.xs+1, (int)mEndPoint.z() -1,3,3,0xffff00);
		}
	}


	// Do not call CWnd::OnPaint() for painting messages
}

UINT BackgroundThread(LPVOID params)
{
	CChildView* view = (CChildView*)params;
	view->mEstimator.Calculate(view->mInputImage, view->mStartPoint, view->mEndPoint);

	GetDispSliceFromTransportFunction(view->mEstimator.GetTransportFunctionCalculator(), Talox, view->m_xsee, view->m_dispd2); // TaloX - enum
	GetDispSliceFromTransportFunction(view->mEstimator.GetTransportFunctionCalculator(), Taloy, view->m_ysee, view->m_dispd1);
	GetDispSliceFromTransportFunction(view->mEstimator.GetTransportFunctionCalculator(), Taloz, view->m_zsee, view->m_disp);
	view->Invalidate();
	_DUMP_PROFILE_INFO(debug_profiling_folder + "profiling.txt");
	return 0;
}

void CChildView::InitThread()
{
	if (mAlgorithmStage) return;
	if (!m_bispoint) {
		mStartPoint.x() = 40;
		mEndPoint.x() = 201;
		mStartPoint.y() = mEndPoint.y() = 87;
		mStartPoint.z() = mEndPoint.z() = 54;
		m_bispoint = 2;
	}
	
	
	if (m_bispoint <= 1) {
		mEndPoint.x() = -100;
		mEndPoint.y() = -100;
	}
	if (m_bispoint == 1) {
		mEndPoint.z() = m_xsee;
	}
	AreaEikonalES::iter_callback_type iter_cb = [this](int iteration) {
		if (iteration % 7 == 0) {
			std::thread invalidate_thread([this]() {Sleep(50); Invalidate(); });
			invalidate_thread.detach();
		}
	};
#ifdef DEBUG_STATES
	//save initial maps
	save_image(debug_states_folder + "ph0_testimage.tif", mInputImage);
	EventSource::data_callback_type dist_init_debug_cb = [](const sitk::Image& map, int idx) {
		save_image(debug_states_folder + "ph0_dist" + std::to_string(idx)+".tif", map);
	};
	mEstimator.HookStageDataInitializedEvent(AreaEikonalStage, dist_init_debug_cb, Distance);

	EventSource::data_callback_type phase_init_debug_cb = [](const sitk::Image& map, int idx) {
		save_image(debug_states_folder + "ph0_field" + std::to_string(idx) + ".tif", map);
	};
	mEstimator.HookStageDataInitializedEvent(AreaEikonalStage, phase_init_debug_cb, PhaseField);

	EventSource::data_callback_type phi_init_debug_cb = [](const sitk::Image& map, int) {
		save_image(debug_states_folder + "ph0_phi.tif", map);
	};
	mEstimator.HookCalculatedPhiMapEvent(phi_init_debug_cb);

	const sitk::Image* final_distmap[2];
	EventSource::data_callback_type distmap_update_debug_cb = [&final_distmap](const sitk::Image& map, int idx) {
		final_distmap[idx] = &map;
	};
	EventSource::basic_callback_type finished_eikonal_debug_cb = [&final_distmap]() {
		save_image(debug_states_folder + "ph1_distance_0.tif", *final_distmap[0]);
		//save_image(debug_states_folder + "ph1_distance_1.tif", *final_distmap[1]);
	};
	mEstimator.HookStageUpdatedEvent(AreaEikonalStage, distmap_update_debug_cb, Distance);
	mEstimator.HookStageFinishedEvent(AreaEikonalStage, finished_eikonal_debug_cb);

	EventSource::data_callback_type init_rot_distmap_debug_cb = [](const sitk::Image& map, int idx) {
		save_image(debug_states_folder + "ph1_r_distance_" + to_string(idx) + ".tif", map);
	};
	mEstimator.HookStageDataInitializedEvent(RotatedAreaEikonalStage, init_rot_distmap_debug_cb, Distance);

	EventSource::data_callback_type rot_distmap_update_debug_cb = [&final_distmap](const sitk::Image& map, int idx) {
		final_distmap[idx] = &map;
	};
	EventSource::basic_callback_type rot_finished_eikonal_debug_cb = [&final_distmap]() {
		save_image(debug_states_folder + "ph1.5_r_distance_0.tif", *final_distmap[0]);
		save_image(debug_states_folder + "ph1.5_r_distance_1.tif", *final_distmap[1]);
	};
	mEstimator.HookStageUpdatedEvent(RotatedAreaEikonalStage, rot_distmap_update_debug_cb, Distance);
	mEstimator.HookStageFinishedEvent(RotatedAreaEikonalStage, rot_finished_eikonal_debug_cb);

	/*EventSource::data_callback_type grad_x_init_debug_cb = [](const sitk::Image& map, int) {
		save_image(debug_states_folder + "ph2_slice_gx.tif", map);
	};
	mEstimator.mInitialContourCalculator.HookInitializeGradXEvent(grad_x_init_debug_cb);

	EventSource::data_callback_type grad_y_init_debug_cb = [](const sitk::Image& map, int) {
		save_image(debug_states_folder + "ph2_slice_gy.tif", map);
	};
	mEstimator.mInitialContourCalculator.HookInitializeGradYEvent(grad_y_init_debug_cb);*/

	const sitk::Image* final_transport;
	EventSource::data_callback_type transport_updated_debug_cb = [&final_transport](const sitk::Image& map, int) {
		final_transport = &map;
	};
	mEstimator.HookStageUpdatedEvent(TransportFunctionStage, transport_updated_debug_cb);
	EventSource::basic_callback_type transport_finished_debug_cb = [&final_transport]() {
		save_image(debug_states_folder + "ph4_transport0.tif", *final_transport);
	};
	mEstimator.HookStageFinishedEvent(TransportFunctionStage, transport_finished_debug_cb);

	EventSource::basic_callback_type calculation_finished_debug_cb = [this]() {
		save_image(debug_states_folder + "ph4_unrot_transportfn.tif", mEstimator.GetTransportFunctionCalculator().GetTransportFunction(0));
		save_image(debug_states_folder + "ph3_transport_gx.tif", mEstimator.GetTransportFunctionCalculator().GetGradX());
		save_image(debug_states_folder + "ph3_transport_gy.tif", mEstimator.GetTransportFunctionCalculator().GetGradY());
		save_image(debug_states_folder + "ph3_transport_gz.tif", mEstimator.GetTransportFunctionCalculator().GetGradZ()); };
#endif
	EventSource::iter_callback_type change_state_cb = [this](int iteration) {
		mAlgorithmStage = iteration;
		switch (iteration) {
		case DISTANCE_ITERATION:
			break;
		case DISTANCE_MEAN_ITERATION:
			break;
		case PLANE_PHASEFIELD_ITERATION:
			break;
		case TRANSPORT_FUNCTION_ITERATION:
#ifdef DEBUG_STATES
			save_image(debug_states_folder + "ph3_slice_field.tif", mEstimator.GetInitialContourCalculator().GetPhaseField());
#endif
			
			break;
		}
		Invalidate(FALSE);
	};

	mEstimator.HookIterationEvent(change_state_cb);

	EventSource::basic_callback_type finished_cb = [this]() {
		mAlgorithmStage = DONE_ITERATION;
	};

	EventSource::data_callback_type field_init_cb = [this](const sitk::Image& map, int idx) {
		if(mPhaseFieldMap[idx] != &map)
			mPhaseFieldMap[idx] = &map;
	};
	EventSource::data_callback_type meeting_points_map_init_cb = [this](const sitk::Image& meetPointsMap, int) {
		if(mMeetingPointsMap != &meetPointsMap)
			mMeetingPointsMap = &meetPointsMap;
	};
	AreaEikonalES::plane_update_callback_type meeting_plane_update_cb = [this](Vec3<double> center, Vec3<double> normal) {
		mPlaneNormal = normal;
		mPlaneOffset = -(center * normal).Sum();
	};
	mEstimator.HookStageIterationEvent(AreaEikonalStage, iter_cb);
	mEstimator.HookStageDataInitializedEvent(AreaEikonalStage, meeting_points_map_init_cb, MeetingPoints);
	mEstimator.HookStageDataInitializedEvent(AreaEikonalStage, field_init_cb, PhaseField);
	mEstimator.HookUpdateMeetingPlaneEvent(meeting_plane_update_cb);

	EventSource::data_callback_type field_2d_update_cb = [this](const sitk::Image& field2D, int) {
		if(mPhaseField2DMap != &field2D)
			mPhaseField2DMap = &field2D;
	};

	mEstimator.HookStageIterationEvent(PlanePhaseFieldStage, iter_cb);
	mEstimator.HookStageDataInitializedEvent(PlanePhaseFieldStage, field_2d_update_cb);
	mEstimator.HookStageUpdatedEvent(PlanePhaseFieldStage, field_2d_update_cb);

	//mEstimator.mAreaEikonal.Initialize(m_imageOp.m_testimage, mStartPoint, mEndPoint, m_expfac, 0.01);
	CWinThread* thread = AfxBeginThread(BackgroundThread,this,THREAD_PRIORITY_NORMAL,0,CREATE_SUSPENDED,0);
	if(thread == 0) return;
	mAlgorithmStage = DISTANCE_ITERATION;
	thread->m_bAutoDelete = true;
	thread->ResumeThread();
	//Invalidate();

}

void CChildView::PauseThread(int threadstat) // 1-(re-)start 2-suspended
{
	if (!mAlgorithmStage) return;
	if (mAlgorithmStage < DONE_ITERATION) {
		mAlgorithmStage = DONE_ITERATION;
	}
	Invalidate();
}

void CChildView::StopThread()
{
	mAlgorithmStage = 0;
	_DUMP_PROFILE_INFO(debug_profiling_folder + "profiling.txt")
	Invalidate();
}

void CChildView::GetAllPathX()
{
	vector<double> plane_center_physical = mPlaneCenter;
	vector<double> plane_center_transformed = mEstimator.GetRotatedAreaEikonal().TransformPhysicalPointToContinuousIndex(plane_center_physical);
	double plane_slice = plane_center_transformed[0];
	vector<unsigned> size = mEstimator.GetRotatedAreaEikonal().GetDistanceMap(0).GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	unordered_set<unsigned long> boundary_points = mEstimator.GetInitialContourCalculator().RetrieveBound();
	mBoundaryContour.clear();
	for (auto& it = boundary_points.begin(); it != boundary_points.end(); it++) {
		int z = *it >> 16;
		int y = *it & 0b1111111111111111;
		vector<double> rotated;
		rotated = mEstimator.GetRotatedAreaEikonal().TransformContinuousIndexToPhysicalPoint({plane_slice, (double)y, (double)z});
		rotated = mEstimator.GetAreaEikonal().TransformPhysicalPointToContinuousIndex(rotated);
		mBoundaryContour.push_back(point_to_representation((int)round(rotated[0]), (int)round(rotated[1]), (int)round(rotated[2])));
	}
	std::vector<POINT3D>& bound = mBoundaryContour;
	int si = (int)bound.size();
	if (si > MAXMINPATH) si = MAXMINPATH;
	for (int jj = 0; jj < MAXMINPATH; ++jj) {
		mEstimator.mMinimalPaths[0].clear();
		mEstimator.mMinimalPaths[1].clear();
	}

	for (int jj = 0; jj < si; ++jj) {
		auto [x, y, z] = representation_to_point<int>(bound[jj]);
		if (z < 2 || z >= zs - 2) continue;
		if (x < 2 || x >= xs - 2) continue;
		if (y < 2 || y >= ys - 2) continue;
		mEstimator.mMinimalPaths[0].push_back(mEstimator.GetAreaEikonal().ResolvePath({x, y, z}, 0));
		mEstimator.mMinimalPaths[1].push_back(mEstimator.GetAreaEikonal().ResolvePath({x, y, z}, 1));

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

		int8_t* x_path_buffer = mPathImageX.GetBufferAsInt8();
		int8_t* y_path_buffer = mPathImageY.GetBufferAsInt8();
		int8_t* z_path_buffer = mPathImageZ.GetBufferAsInt8();
		vector<unsigned> size = mPathImageX.GetSize();
		memset(x_path_buffer, -100, size[0] * size[1] * size[2]);
		memset(y_path_buffer, -100, size[0] * size[1] * size[2]);
		memset(z_path_buffer, -100, size[0] * size[1] * size[2]);
	
	if (mAlgorithmStage == DONE_ITERATION) {
		if (point.x < m_disp.xs) {
			mEstimator.mMinimalPaths[0].push_back(mEstimator.GetAreaEikonal().ResolvePath({point.x, point.y, m_zsee}, 0));
			mEstimator.mMinimalPaths[1].push_back(mEstimator.GetAreaEikonal().ResolvePath({ point.x, point.y, m_zsee }, 1));
		}
		else if (point.x < 2 * m_disp.xs + 1) {
			mEstimator.mMinimalPaths[0].push_back(mEstimator.GetAreaEikonal().ResolvePath({ point.x - m_disp.xs - 1, m_ysee, point.y }, 0));
			mEstimator.mMinimalPaths[1].push_back(mEstimator.GetAreaEikonal().ResolvePath({ point.x - m_disp.xs - 1, m_ysee, point.y }, 1));
		}
		else if(point.x < m_disp.xs+m_dispd1.xs+m_dispd2.xs) {
			mEstimator.mMinimalPaths[0].push_back(mEstimator.GetAreaEikonal().ResolvePath({ m_xsee, point.x - 2 * (m_disp.xs + 1), point.y }, 0));
			mEstimator.mMinimalPaths[1].push_back(mEstimator.GetAreaEikonal().ResolvePath({ m_xsee, point.x - 2 * (m_disp.xs + 1), point.y }, 1));
		}
		Invalidate();
	}
	else if (!mAlgorithmStage) {
		vector<uint32_t> size = mInputImage.GetSize();
		int xs(size[0]), ys(size[1]);
		if (m_disp.xs && point.x < xs && point.y < ys) {
			if (!m_bispoint) {
				mStartPoint = Vec3<double>(point.x, point.y, m_zsee);
				++m_bispoint;
			}
			else {
				mEndPoint = Vec3<double>(point.x, point.y, m_zsee);
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
			m_pView->mInputImage = sitk::RescaleIntensity(ci, 0, 1);
			vector<uint32_t> size = m_pView->mInputImage.GetSize();
			int xs = size[0], ys = size[1], zs = size[2];
			//sitk_2_vox_img(ci, m_pView->m_imageOp.m_testimage);
			m_pView->mInputImage = m_pView->m_imageOp.Blur(m_pView->mInputImage);
			m_pView->mInputImage = m_pView->m_imageOp.Blur(m_pView->mInputImage);
			m_pView->mInputImage = m_pView->m_imageOp.Blur(m_pView->mInputImage);
			m_pView->mInputImage = m_pView->m_imageOp.Blur(m_pView->mInputImage);
			m_pView->mPathImageX = sitk::Image(size, sitk::sitkInt8)-100;
			m_pView->mPathImageY = sitk::Image(size, sitk::sitkInt8)-100;
			m_pView->mPathImageZ = sitk::Image(size, sitk::sitkInt8)-100;
			//CImage ci; ci.Load(LPCTSTR(cfn));

			if (!xs || !ys || !zs) return;
			const double* im_buffer = m_pView->mInputImage.GetBufferAsDouble();
			m_pView->m_work.Set(xs, ys);
			for (int yy = 0; yy < ys; yy++) {
				for (int xx = 0; xx < xs; xx++) {
					m_pView->m_work[yy][xx] = BUF_IDX3D(im_buffer, xs, ys, zs, xx, yy, m_pView->m_zsee);
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
	m_pView->mInputImage = m_pView->m_imageOp.CreateTestImage();
	
	vector<uint32_t> size = m_pView->mInputImage.GetSize();
	m_pView->mPathImageX = sitk::Image(size, sitk::sitkInt8) - 100;
	m_pView->mPathImageY = sitk::Image(size, sitk::sitkInt8) - 100;
	m_pView->mPathImageZ = sitk::Image(size, sitk::sitkInt8) - 100;
	int xs = size[0], ys = size[1], zs = size[2]; 
	double* im_buffer = m_pView->mInputImage.GetBufferAsDouble();
	m_pView->m_work.Set(xs, ys);
	for (int yy = 0; yy < ys; yy++) {
		for (int xx = 0; xx < xs; xx++) {
			m_pView->m_work[yy][xx] = BUF_IDX3D(im_buffer, xs, ys, zs, xx, yy, m_pView->m_zsee);
		}
	}
	m_pView->m_valid = 1;
	m_pView->m_grays = true;
	m_pView->m_color = false;
	
	{
		m_zlevel.SetRange(0,zs-1,1); // from z
		m_cstartdir.SetRange(0,ys-1,1); // from y
		m_cstopdir.SetRange(0,xs-1,1); // from x
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
	//m_pView->m_imageOp.GetXTestBound(m_pView->m_xsee/* - 1*/, m_pView->mEstimator.mBoundaryContour);
	//if (m_pView->mEstimator.mAreaEikonal.m_distance[0].xs > 0) {
	//	int xs(m_pView->mEstimator.mAreaEikonal.m_distance[0].xs), ys(m_pView->mEstimator.mAreaEikonal.m_distance[0].ys), zs(m_pView->mEstimator.mAreaEikonal.m_distance[0].zs);
	//	SVoxImg<SWorkImg<realnum>> & passi = m_pView->m_imageOp.GetIniMap(xs, ys, zs, m_pView->m_xsee/* - 1*/);
	//	//m_pview->mEstimator.mTransportFunctionCalculator.Initialize(m_pView->mEstimator.mAreaEikonal.m_smoothdist[0],passi, m_pView->mEstimator.mCurrentDistance[0]);
	//	m_pview->mEstimator.mTransportFunctionCalculator.Initialize(m_pView->mEstimator.mAreaEikonal.m_distance[0], passi, m_pView->mEstimator.mAreaEikonal.mCurrentDistance);

		//--m_pView->m_xsee;
		/*
		m_pview->mEstimator.mTransportFunctionCalculator.GetDispSlice(Talox, m_pView->m_xsee, m_pView->m_dispd2); // TaloX - enum
		m_pview->mEstimator.mTransportFunctionCalculator.GetDispSlice(Taloy, m_pView->m_ysee, m_pView->m_dispd1);
		m_pview->mEstimator.mTransportFunctionCalculator.GetDispSlice(Taloz, m_pView->m_zsee, m_pView->m_disp);*/
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
	vector<uint32_t> size = m_pView->mInputImage.GetSize();
	if (!m_bini) goto ret;
		m_pView->m_zsee = m_zlevel.GetPos();
		char txt[22];
		sprintf_s(txt,20,"%d",m_pView->m_zsee);
		m_ezlevel.SetWindowTextA(txt);

		int xs = size[0], ys = size[1], zs = size[2];
		double* im_buffer = m_pView->mInputImage.GetBufferAsDouble();
		m_pView->m_work.Set(xs, ys);
		for (int yy = 0; yy < ys; yy++) {
			for (int xx = 0; xx < xs; xx++) {
				m_pView->m_work[yy][xx] = BUF_IDX3D(im_buffer, xs, ys, zs, xx, yy, m_pView->m_zsee);
			}
		}

	if (!m_pView->m_btransportview) {
		m_pView->m_work.GetDispImg(m_pView->m_disp);
	}
	else {
		GetDispSliceFromTransportFunction(m_pView->mEstimator.GetTransportFunctionCalculator(), Taloz, m_pView->m_zsee, m_pView->m_disp);
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
	sitk::Image &testimg = m_pView->mInputImage;
	double* testimg_buffer = testimg.GetBufferAsDouble();

	if (!m_pView->m_btransportview) {
		vector<uint32_t> size = m_pView->mInputImage.GetSize();
		int xs = size[0], ys = size[1], zs = size[2];
		SWorkImg<realnum>& yslice = m_pView->m_intey;
		yslice.Set(xs, zs);
		int ylev = m_cstartdir.GetPos();

		for (int zz = 0; zz < zs; ++zz) {
			for (int xx = 0; xx < xs; ++xx) {
				yslice[zz][xx] = BUF_IDX3D(testimg_buffer, xs, ys, zs, xx, ylev, zz);
			}
		}

		yslice.GetDispImg(m_pView->m_dispd1);
		m_pView->m_ysee = m_cstartdir.GetPos();
	}
	else {
		m_pView->m_ysee = m_cstartdir.GetPos();
		GetDispSliceFromTransportFunction(m_pView->mEstimator.GetTransportFunctionCalculator(), Taloy, m_pView->m_ysee, m_pView->m_dispd1);
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


	sitk::Image &testimg = m_pView->mInputImage;
	double* testimg_buffer = testimg.GetBufferAsDouble();
	if (!m_pView->m_btransportview) {
		vector<uint32_t> size = m_pView->mInputImage.GetSize();
		int xs = size[0], ys = size[1], zs = size[2]; 
		SWorkImg<realnum>& xslice = m_pView->m_intex;
		xslice.Set(ys, zs);
		int xlev = m_cstopdir.GetPos();

		for (int zz = 0; zz < zs; ++zz) {
			for (int yy = 0; yy < ys; ++yy) {
				xslice[zz][yy] = BUF_IDX3D(testimg_buffer, xs, ys, zs, xlev, yy, zz);
			}
		}

		xslice.GetDispImg(m_pView->m_dispd2);
		m_pView->m_xsee = m_cstopdir.GetPos();
	}
	else {
		m_pView->m_xsee = m_cstopdir.GetPos();
		GetDispSliceFromTransportFunction(m_pView->mEstimator.GetTransportFunctionCalculator(), Talox, m_pView->m_xsee, m_pView->m_dispd2);
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
	m_pView->m_imageOp.CreateTestImage2();
	{
		vector<uint32_t> size = m_pView->mInputImage.GetSize();
		int xs = size[0], ys = size[1], zs = size[2];
		double* im_buffer = m_pView->mInputImage.GetBufferAsDouble();
		m_pView->m_work.Set(xs, ys);
		for (int yy = 0; yy < ys; yy++) {
			for (int xx = 0; xx < xs; xx++) {
				m_pView->m_work[yy][xx] = BUF_IDX3D(im_buffer, xs, ys, zs, xx, yy, m_pView->m_zsee);
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
	if (m_pView->mAlgorithmStage) return;
	static bool g_modeswitch = false;
	g_modeswitch = !g_modeswitch;
	m_pView->mEstimator.SetUsesCorrection(g_modeswitch);
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

		GetDispSliceFromTransportFunction(m_pView->mEstimator.GetTransportFunctionCalculator(), Talox, m_pView->m_xsee, m_pView->m_dispd2);
		GetDispSliceFromTransportFunction(m_pView->mEstimator.GetTransportFunctionCalculator(), Taloy, m_pView->m_ysee, m_pView->m_dispd1);
		GetDispSliceFromTransportFunction(m_pView->mEstimator.GetTransportFunctionCalculator(), Taloz, m_pView->m_zsee, m_pView->m_disp);

	}

	m_pView->Invalidate();
}


void CControlDlg::OnBnClickedTestBuild()
{
	// TODO: Add your control notification handler code here
	/*m_pView->mEstimator.mTransportFunctionCalculator.Calculate();
	if (m_pView->m_btransportview) {

		m_pView->mEstimator.mTransportFunctionCalculator.GetDispSlice(Talox, m_pView->m_xsee, m_pView->m_dispd2);
		m_pView->mEstimator.mTransportFunctionCalculator.GetDispSlice(Taloy, m_pView->m_ysee, m_pView->m_dispd1);
		m_pView->mEstimator.mTransportFunctionCalculator.GetDispSlice(Taloz, m_pView->m_zsee, m_pView->m_disp);

	}
	m_pView->Invalidate();*/
}

void GetDispSliceFromTransportFunction(const TransportFunction& transportFunction, int along, int at, SDisImg& r) // better colorization!
{
	if (!transportFunction.m_active) return;
	realnum zlim(1e-22);
	int bcol(0);
	const sitk::Image& trf = transportFunction.GetTransportFunction(0);
	const double* trf_buffer = trf.GetBufferAsDouble();
	const int* bound_buffer = transportFunction.GetReadOnlyMap().GetBufferAsInt32();
	std::vector<uint32_t> size = trf.GetSize();
	int xs(size[0]), ys(size[1]), zs(size[2]);
	double min_val(trf_buffer[0]); // boundary pixels were initialized to the minimum an won't be updated

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
				if (BUF_IDX3D(bound_buffer, xs, ys, zs, xx, yy, zz) == -1) // not initialized
					r[zz][yy] = (0xff << 16) + (0xff << 8) + 0xff;
				else { //if (mReadOnlyMap[zz][yy][xx] == 1) // boundary or internal
					if (BUF_IDX3D(trf_buffer, xs, ys, zs, xx, yy, zz) < -zlim) // -mMin: maximal value
						r[zz][yy] = (bcol << 16) + (bcol << 0) + (int)256 * (127 + 0xff * (int)((0.5 * BUF_IDX3D(trf_buffer, xs, ys, zs, xx, yy, zz) / min_val)));
					else if (BUF_IDX3D(trf_buffer, xs, ys, zs, xx, yy, zz) > zlim)
						r[zz][yy] = ((int)(127 - 0xff * (int)(0.5 * BUF_IDX3D(trf_buffer, xs, ys, zs, xx, yy, zz) / min_val)) << 16) + (bcol << 8) + bcol;
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
				if (BUF_IDX3D(bound_buffer, xs, ys, zs, xx, yy, zz) == -1) // not initialized
					r[zz][xx] = (0xff << 16) + (0xff << 8) + 0xff;
				else { //if (mReadOnlyMap[zz][yy][xx] == 1) // boundary or internal
					if (BUF_IDX3D(trf_buffer, xs, ys, zs, xx, yy, zz) < -zlim) // -mMin: maximal value
						r[zz][xx] = (bcol << 16) + (bcol << 0) + 256 * (int)(127 + 0xff * ((0.5 * BUF_IDX3D(trf_buffer, xs, ys, zs, xx, yy, zz) / min_val)));
					else if (BUF_IDX3D(trf_buffer, xs, ys, zs, xx, yy, zz) > zlim)
						r[zz][xx] = ((int)(127 - 0xff * (int)(0.5 * BUF_IDX3D(trf_buffer, xs, ys, zs, xx, yy, zz) / min_val)) << 16) + (bcol << 8) + bcol;
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
				if (BUF_IDX3D(bound_buffer, xs, ys, zs, xx, yy, zz) == -1) // not initialized
					r[yy][xx] = (0xff << 16) + (0xff << 8) + 0xff;
				else { //if (mReadOnlyMap[zz][yy][xx] == 1) // boundary or internal
					if (BUF_IDX3D(trf_buffer, xs, ys, zs, xx, yy, zz) < -zlim) // -mMin: maximal value
						r[yy][xx] = (bcol << 16) + (bcol << 0) + 256 * (int)(127 + 0xff * ((0.5 * BUF_IDX3D(trf_buffer, xs, ys, zs, xx, yy, zz) / min_val)));
					else if (BUF_IDX3D(trf_buffer, xs, ys, zs, xx, yy, zz) > zlim)
						r[yy][xx] = ((int)(127 - 0xff * (int)(0.5 * BUF_IDX3D(trf_buffer, xs, ys, zs, xx, yy, zz) / min_val)) << 16) + (bcol << 8) + bcol;
					else
						r[yy][xx] = (0xff << 16) + (0xff << 8) + 0xff;
				}
			}
		}

	}


}