
// ChildView.cpp : implementation of the CChildView class
//

#include "stdafx.h"
#include "MinArea.h"
#include "ChildView.h"
#include "FitPlane.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif
#include <iostream>
using namespace std;
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
}

CChildView::~CChildView()
{
}


BEGIN_MESSAGE_MAP(CChildView, CWnd)
	ON_WM_PAINT()
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
			std::vector<CVec3> &bound = m_liftedEikonal.m_boundcontour;
			int si = bound.size();
			if (si) {
				for (int ii = 0; ii < si; ++ii) {
					dc.SetPixelV(offx+(int)bound[ii].y,(int)bound[ii].z,0xffff);
				}
				//SWorkImg<realnum>& ar = m_liftedEikonal.m_imageOp.GetXTestSigDis(); // just returns m_loc
				//ar.GetDispImg(m_dispd3);
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

	if (m_threadactivated)  {
		//char txt[100];
		//sprintf_s(txt, 99, "%d %d         ", g_ign, g_ian);
		//dc.TextOutA(300, 10, txt);
		for (int ii = 0; ii < 2; ++ii) {
			SVoxImg<SWorkImg<realnum>> &field = m_liftedEikonal.m_phasefield.m_field[ii];
			SVoxImg<SWorkImg<int>>& meeting_plane = m_liftedEikonal.m_phasefield.meeting_plane_positions;

			SVoxImg<SWorkImg<realnum>> &dismap = (!m_bsee && !m_curvsee) ? m_liftedEikonal.m_phasefield.m_distance[ii]:
				(!m_bsee && m_curvsee == 1) ? m_liftedEikonal.m_phasefield.m_Sumcurvature:
				(!m_bsee && m_curvsee == 2) ?m_liftedEikonal.m_phasefield.m_Sumcurvature:
				m_liftedEikonal.m_phasefield.m_smoothdist
			;
			/*m_Sumcurvature[ii],m_expdist[ii],m_smoothdist[ii],m_distance[ii],m_thicknes[ii]*/

			int zsee = 0;
			if (m_pControl) {
				//char distex[22]; sprintf_s(distex,21,"%f",m_liftedEikonal.m_currentdistance[0]); m_pControl->m_cdist.SetWindowTextA(distex);
				zsee = m_zsee;
			}


			if (m_threadactivated != 3) {
				if (m_disp.xs && !m_btransportview) {
					for (int yy = 0; yy < field.ys; ++yy) {
						for (int xx = 0; xx < field.xs; ++xx) {
							if (m_threadactivated == 1 || m_threadactivated == 5) {
								if (abs(plane_normal.second.x*xx+plane_normal.second.y*yy+plane_normal.second.z*zsee+plane_offset) < 1.0) {
									dc.SetPixelV(xx, yy, 0xff00ff);
								}
								else if (meeting_plane[zsee][yy][xx]) {
									dc.SetPixelV(xx, yy, 0xff0000);
								}
								else if (field[zsee][yy][xx] > -0.6f && field[zsee][yy][xx] < 0.6f)
									dc.SetPixelV(xx, yy, !ii ? 0xff : 0xff00);
							}
						}
					}
				}
			}

			if (m_dispd1.xs) for (int zz = 0; zz < field.zs; ++zz) {
				for (int xx = 0; xx < field.xs; ++xx) {
					if (m_threadactivated  == 1 || m_threadactivated == 5) {
						if (abs(plane_normal.second.x * xx + plane_normal.second.y * m_ysee + plane_normal.second.z * zz + plane_offset) < 1.0) {
							dc.SetPixelV(m_disp.xs + 1 + xx, zz, 0xff00ff);
						}
						else if (meeting_plane[zz][m_ysee][xx]) {
							dc.SetPixelV(m_disp.xs + 1 + xx, zz, 0xff0000);
						}
						else if (field[zz][m_ysee][xx] > -0.6f && field[zz][m_ysee][xx] < 0.6f)
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
				for (int zz = 0; zz < field.zs; ++zz) {
					for (int yy = 0; yy < field.ys; ++yy) {
						if (m_threadactivated == 1 || m_threadactivated == 5) {
							if (abs(plane_normal.second.x * m_xsee + plane_normal.second.y * yy + plane_normal.second.z * zz + plane_offset) < 1.0) {
								dc.SetPixelV(2 * (m_disp.xs + 1) + yy, zz, 0xff00ff);
							}
							else if (meeting_plane[zz][yy][m_xsee]) {
								dc.SetPixelV(2 * (m_disp.xs + 1) + yy, zz, 0xff0000);
							}
							else if (field[zz][yy][m_xsee] > -0.6f && field[zz][yy][m_xsee] < 0.6f)
								dc.SetPixelV(2 * (m_disp.xs + 1) + yy, zz, !ii ? 0xff : 0xff00);
						}
						if (m_threadactivated == 3) {
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
		}


		{
			SVoxImg<SWorkImg<realnum>> &field = m_liftedEikonal.m_phasefield.m_field[0];

			int xx(0), yy(0), zz(0);
			dc.MoveTo(xx,300-(int)(10*field[m_zsee][m_ysee][xx]));
			for (xx = 1; xx < field.xs; ++xx) dc.LineTo(xx*5,300-(int)(10*field[m_zsee][m_ysee][xx]));
			dc.MoveTo(yy,350-(int)(10*field[m_zsee][yy][m_xsee]));
			for (yy = 1; yy < field.ys; ++yy) dc.LineTo(yy*5,350-(int)(10*field[m_zsee][yy][m_xsee]));
			dc.MoveTo(zz,400-(int)(10*field[zz][m_ysee][m_xsee]));
			for (zz = 1; zz < field.zs; ++zz) dc.LineTo(zz*5,400-(int)(10*field[zz][m_ysee][m_xsee]));

		}

	}
			

	if (m_disp.xs) {
		if (m_bispoint > 0) {
			if (abs(m_zsee-m_zdeparture) < 3) dc.FillSolidRect(m_fieldinit.x-1,m_fieldinit.y-1,3,3,0xff);
			if (abs(m_ysee-m_fieldinit.y) < 3) dc.FillSolidRect(m_fieldinit.x+m_disp.xs,m_zdeparture-1,3,3,0xff);
			if (abs(m_xsee-m_fieldinit.x) < 3) dc.FillSolidRect(m_fieldinit.y+2*m_disp.xs+1,m_zdeparture-1,3,3,0xff);
		}
		if (m_bispoint > 1) {
			if (abs(m_zsee-m_zarrival) < 3) dc.FillSolidRect(m_arrival.x-1,m_arrival.y-1,3,3,0xffff00);
			if (abs(m_ysee-m_arrival.y) < 3) dc.FillSolidRect(m_arrival.x+m_disp.xs,m_zarrival-1,3,3,0xffff00);
			if (abs(m_xsee-m_arrival.x) < 3) dc.FillSolidRect(m_arrival.y+2*m_disp.xs+1,m_zarrival-1,3,3,0xffff00);
		}
	}


	// Do not call CWnd::OnPaint() for painting messages
}

void CChildView::InitTransport()
{
	//m_liftedEikonal.m_imageOp.GetXTestBound(m_xsee, m_liftedEikonal.m_boundcontour);
	m_liftedEikonal.m_imageOp.GetPlaneDistMap(m_liftedEikonal.m_inicountourCalculator.RetrieveBound());//
	if (m_liftedEikonal.m_phasefield.m_distance[0].xs > 0) {
		SVoxImg<SWorkImg<realnum>>& passi = m_liftedEikonal.m_imageOp.GetIniMap(plane_normal.first.x);
		//realnum maxdist = min(m_liftedEikonal.m_currentdistance[0], m_liftedEikonal.m_currentdistance[1]);
		realnum maxdist = m_liftedEikonal.m_phasefield.m_currentdistance;
		m_transport.TrInit(m_liftedEikonal.m_phasefield.m_combined_distance, passi, maxdist);

		m_transport.GetDispSlice(Talox, m_xsee, m_dispd2); // TaloX - enum
		m_transport.GetDispSlice(Taloy, m_ysee, m_dispd1);
		m_transport.GetDispSlice(Taloz, m_zsee,m_disp);
	}

	Invalidate();
	m_transport.TrControl(333);
	Invalidate();

}

UINT BackgroundThread(LPVOID params)
{
	CChildView* view = (CChildView*)params;
	int cyc = 0;
	extern bool g_modeswitch;

	while (view->m_threadactivated) {

		if (view->m_threadactivated == 1) { // 3D

			view->m_liftedEikonal.m_phasefield.Iterate(g_modeswitch);
			if (++cyc > 7) { view->Invalidate(FALSE); cyc = 0; }

			if (view->m_liftedEikonal.m_phasefield.m_bdone) {
				for (int zz = 1; zz < view->m_liftedEikonal.m_phasefield.meeting_plane_positions.zs - 1; ++zz) {
					for (int yy = 1; yy < view->m_liftedEikonal.m_phasefield.meeting_plane_positions.ys - 1; ++yy) {
						for (int xx = 1; xx < view->m_liftedEikonal.m_phasefield.meeting_plane_positions.xs - 1; ++xx) {
							if (view->m_liftedEikonal.m_phasefield.meeting_plane_positions[zz][yy][xx])
								view->m_liftedEikonal.meeting_plane.insert(IPoi3<double>(xx, yy, zz));
							view->m_liftedEikonal.m_phasefield.m_combined_distance[zz][yy][xx] =
								view->m_liftedEikonal.m_phasefield.m_distance[0][zz][yy][xx] >= 0 ? view->m_liftedEikonal.m_phasefield.m_distance[0][zz][yy][xx] : view->m_liftedEikonal.m_phasefield.m_distance[1][zz][yy][xx];
						}
					}
				}
				view->plane_normal = best_plane_from_points(view->m_liftedEikonal.meeting_plane);
				view->plane_offset = -(view->plane_normal.first.x * view->plane_normal.second.x + view->plane_normal.first.y * view->plane_normal.second.y + view->plane_normal.first.z * view->plane_normal.second.z);
				view->m_threadactivated = 2;
				view->Invalidate(FALSE);
				
			}

		}
		else if (view->m_threadactivated == 2) {
			view->m_liftedEikonal.m_inicountourCalculator.GetDistancemean(view->m_liftedEikonal.m_phasefield.m_combined_distance, view->plane_normal.first.x);
			Sleep(300);
			if (view->m_liftedEikonal.m_inicountourCalculator.m_bInited)
				view->m_threadactivated = 3;
		}
		else if (view->m_threadactivated == 3) { // plane
			view->m_liftedEikonal.m_inicountourCalculator.Iterate();
			if (!view->m_liftedEikonal.m_inicountourCalculator.m_nect)
				view->m_threadactivated = 4;
			if (++cyc > 7) { view->Invalidate(FALSE); cyc = 0; }
		}
		else if (view->m_threadactivated == 4) {
			view->InitTransport();
			//view->m_prevthreadactivated = 2;
			Sleep(300);
			view->m_threadactivated = 5;
		}
		else if (view->m_threadactivated == 5) {
			Sleep(300);
		}

	}

	return 0;
}

void CChildView::InitThread()
{
	if (m_threadactivated) return;
	if (!m_bispoint) return;//m_fieldinit = CPoint(m_liftedEikonal.m_work.xs/2,m_liftedEikonal.m_work.ys/2);
	
	if (m_bispoint <= 1) { m_arrival = CPoint(-100,-100); }
	if (m_bispoint == 1) {
		int xto = m_xsee; //if (xto) ++xto;
		m_liftedEikonal.PhaseInit(IPoi(m_fieldinit.x, m_fieldinit.y), IPoi(m_arrival.x, m_arrival.y), m_zdeparture, xto, xto);
	}
	else
		//m_liftedEikonal.PhaseInit(IPoi(m_fieldinit.x,m_fieldinit.y),IPoi(m_arrival.x,m_arrival.y),m_zdeparture, m_xsee);
		m_liftedEikonal.PhaseInit(IPoi(m_fieldinit.x,m_fieldinit.y),IPoi(m_arrival.x,m_arrival.y),m_zdeparture,m_zarrival, m_xsee);

	CWinThread* thread = AfxBeginThread(BackgroundThread,this,THREAD_PRIORITY_NORMAL,0,CREATE_SUSPENDED,0);
	if(thread == 0) return;
	m_threadactivated = 1;
	thread->m_bAutoDelete = true;
	thread->ResumeThread();
	Invalidate();

}

void CChildView::PauseThread(int threadstat) // 1-(re-)start 2-suspended
{
	if (!m_threadactivated) return;
	if (m_threadactivated < 5) {
		m_prevthreadactivated = m_threadactivated;
		if (m_prevthreadactivated == 3) m_prevthreadactivated = 2;
		m_threadactivated = 5;
	}
	else {
		m_threadactivated = m_prevthreadactivated;
	}
	Invalidate();
}

void CChildView::StopThread()
{
	m_threadactivated = 0;
	Invalidate();
}


void CChildView::OnLButtonUp(UINT nFlags, CPoint point)
{
	// TODO: Add your message handler code here and/or call default

	if (!m_threadactivated) {
		if (m_disp.xs) {
			if (!m_bispoint) {
				m_fieldinit = point; m_zdeparture = m_zsee; ++m_bispoint;
			}
			else {
				m_arrival = point; m_zarrival = m_zsee; ++m_bispoint;
			}
		}
		Invalidate();
	}

	CWnd::OnLButtonUp(nFlags, point);
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
	ON_BN_CLICKED(IDC_BUTTON1, &CControlDlg::OnBnClickedButton1)
	ON_BN_CLICKED(IDC_BUTTON2, &CControlDlg::OnBnClickedButton2)
	ON_BN_CLICKED(IDC_BUTTON3, &CControlDlg::OnBnClickedButton3)
	ON_BN_CLICKED(IDC_BUTTON4, &CControlDlg::OnBnClickedButton4)
	ON_BN_CLICKED(IDC_BUTTON5, &CControlDlg::OnBnClickedButton5)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER1, &CControlDlg::OnNMCustomdrawSlider1)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER2, &CControlDlg::OnNMCustomdrawSlider2)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER3, &CControlDlg::OnNMCustomdrawSlider3)
	ON_EN_CHANGE(IDC_EDIT3, &CControlDlg::OnEnChangeEdit3)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_SLIDER4, &CControlDlg::OnNMCustomdrawSlider4)
	ON_EN_CHANGE(IDC_EDIT2, &CControlDlg::OnEnChangeEdit2)
	ON_BN_CLICKED(IDC_CHECK1, &CControlDlg::OnBnClickedCheck1)
	ON_EN_CHANGE(IDC_EDIT1, &CControlDlg::OnEnChangeEdit1)
	ON_BN_CLICKED(IDC_BUTTON6, &CControlDlg::OnBnClickedButton6)
	ON_BN_CLICKED(IDC_CHECK2, &CControlDlg::OnBnClickedCheck2)
	ON_BN_CLICKED(IDC_BUTTON7, &CControlDlg::OnBnClickedButton7)
	ON_BN_CLICKED(IDC_BUTTON8, &CControlDlg::OnBnClickedButton8)
	ON_BN_CLICKED(IDC_CHECK3, &CControlDlg::OnBnClickedCheck3)
	ON_BN_CLICKED(IDC_BUTTON10, &CControlDlg::OnBnClickedButton10)
END_MESSAGE_MAP()


// CControlDlg message handlers

void CControlDlg::OnBnClickedButton1()
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

			CImage ci; ci.Load(LPCTSTR(cfn));
			int xs = ci.GetWidth(), ys = ci.GetHeight();
			if (!xs || !ys) return;
			int pitch = ci.GetPitch();
			byte *buf = (byte*)ci.GetBits();
			int mod = abs(pitch)%xs;

			if (!pos && b1st) {
				bool bcolor = false;
				if (ci.GetBPP() == 8) {
					m_pView->m_liftedEikonal.m_img.Set(xs,ys,buf,mod,pitch);
				}
				else if (ci.GetBPP() == 24) {
					m_pView->m_liftedEikonal.m_img.SetColor(xs,ys,buf,mod,pitch);
					bool bok = m_pView->m_liftedEikonal.m_workr.GetAligned(m_pView->m_liftedEikonal.m_img,8,0);
					bok = m_pView->m_liftedEikonal.m_workg.GetAligned(m_pView->m_liftedEikonal.m_img,8,1);
					bok = m_pView->m_liftedEikonal.m_workb.GetAligned(m_pView->m_liftedEikonal.m_img,8,2);
					if (!bok) m_pView->m_liftedEikonal.m_valid = 0;
					else {
						m_pView->m_liftedEikonal.m_workb.GetDispChannel(m_pView->m_disp,0);
						m_pView->m_liftedEikonal.m_workg.GetDispChannel(m_pView->m_disp,1);
						m_pView->m_liftedEikonal.m_workr.GetDispChannel(m_pView->m_disp,2);
						m_pView->m_liftedEikonal.m_valid = 1;
						m_pView->m_liftedEikonal.m_color = true;
						m_pView->m_liftedEikonal.m_grays = false;
					}
					bcolor = true;
				}
				else {
					m_pView->m_liftedEikonal.m_valid = 0;	
					return;
				}
				if (!bcolor) {
					bool bok = m_pView->m_liftedEikonal.m_work.GetAligned(m_pView->m_liftedEikonal.m_img,8); // normed 0-1
					if (bok) {
						m_pView->m_liftedEikonal.m_work.GetDispImg(m_pView->m_disp);
						m_pView->m_liftedEikonal.m_valid = 1;
						m_pView->m_liftedEikonal.m_grays = true;
						m_pView->m_liftedEikonal.m_color = false;
					}
					else m_pView->m_liftedEikonal.m_valid = 0;
				}
			}
		
			b1st = false;

		} while(pos);
		// some init ok here
		m_zlevel.SetRange(0,ZS_-1,1);
		m_cstartdir.SetRange(0,ZS_-1,1);
		m_cstopdir.SetRange(0,ZS_-1,1);
		//m_cdatafac.SetRange(1+10,10+10,1);
		m_cdatafac.SetRange(1,10,1);
		m_cdatafac.SetPos(m_pView->m_liftedEikonal.m_expfac);
		char txt[22];
		sprintf_s(txt,20,"%d",m_flowray);
		m_cflowray.SetWindowTextA(txt);
		m_bini = true;

		m_pView->Invalidate();


	}


}

void CControlDlg::OnBnClickedButton2() // synth. image 1
{
	// TODO: Add your control notification handler code here
	if (m_bini) return;
	m_pView->m_liftedEikonal.m_imageOp.CreateTestImage(XS_,YS_,ZS_);
	{
		m_pView->m_liftedEikonal.m_work = m_pView->m_liftedEikonal.m_imageOp.m_testimage[m_pView->m_zsee];
		m_pView->m_liftedEikonal.m_valid = 1;
		m_pView->m_liftedEikonal.m_grays = true;
		m_pView->m_liftedEikonal.m_color = false;
	}
	{
		m_zlevel.SetRange(0,ZS_-1,1); // from z
		m_cstartdir.SetRange(0,YS_-1,1); // from y
		m_cstopdir.SetRange(0,XS_-1,1); // from x
		m_bini = true;

		m_cdatafac.SetRange(1,10,1); m_cdatafac.SetPos(m_pView->m_liftedEikonal.m_expfac);
		char txt[22]; sprintf_s(txt,20,"%d",m_flowray);	m_cflowray.SetWindowTextA(txt);
	}
	m_pView->m_liftedEikonal.m_work.GetDispImg(m_pView->m_disp);
	OnNMCustomdrawSlider2(0,0);
	OnNMCustomdrawSlider3(0,0);
	m_pView->Invalidate();

}

void CControlDlg::OnBnClickedButton3() // contour on plane
{
	// TODO: Add your control notification handler code here
	m_pView->m_liftedEikonal.m_imageOp.GetXTestBound(m_pView->m_xsee/* - 1*/, m_pView->m_liftedEikonal.m_boundcontour);
	if (m_pView->m_liftedEikonal.m_phasefield.m_distance[0].xs > 0) {
		SVoxImg<SWorkImg<realnum>> & passi = m_pView->m_liftedEikonal.m_imageOp.GetIniMap(m_pView->m_xsee/* - 1*/);
		//m_pView->m_transport.TrInit(m_pView->m_liftedEikonal.m_phasefield.m_smoothdist[0],passi, m_pView->m_liftedEikonal.m_currentdistance[0]);
		m_pView->m_transport.TrInit(m_pView->m_liftedEikonal.m_phasefield.m_distance[0], passi, m_pView->m_liftedEikonal.m_phasefield.m_currentdistance);

		//--m_pView->m_xsee;
		m_pView->m_transport.GetDispSlice(Talox, m_pView->m_xsee, m_pView->m_dispd2); // TaloX - enum
		m_pView->m_transport.GetDispSlice(Taloy, m_pView->m_ysee, m_pView->m_dispd1);
		m_pView->m_transport.GetDispSlice(Taloz, m_pView->m_zsee, m_pView->m_disp);
	}

	m_pView->Invalidate();
}

void CControlDlg::OnBnClickedButton4()
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


void CControlDlg::OnBnClickedButton5() // pause/continue
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
	if (!m_bini) goto ret;
		m_pView->m_zsee = m_zlevel.GetPos();
		char txt[22];
		sprintf_s(txt,20,"%d",m_pView->m_zsee);
		m_ezlevel.SetWindowTextA(txt);

	m_pView->m_liftedEikonal.m_work = m_pView->m_liftedEikonal.m_imageOp.m_testimage[m_pView->m_zsee];
	/*if (!m_pView->m_bsee)
		m_pView->m_liftedEikonal.m_work = m_pView->m_liftedEikonal.m_imageOp.m_testimage[m_pView->m_zsee];
	else
		m_pView->m_liftedEikonal.m_work = m_pView->m_liftedEikonal.m_imageOp.m_testinput[m_pView->m_zsee];
	*/
	if (!m_pView->m_btransportview) {
		m_pView->m_liftedEikonal.m_work.GetDispImg(m_pView->m_disp);
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
	SVoxImg<SWorkImg<realnum>> &testimg = m_pView->m_liftedEikonal.m_imageOp.m_testimage;
	/*
	SVoxImg<SWorkImg<realnum>> &testimg = m_pView->m_bsee?m_pView->m_liftedEikonal.m_imageOp.m_testinput:m_pView->m_liftedEikonal.m_imageOp.m_testimage;
	*/
	if (!m_pView->m_btransportview) {
		int xs = testimg.xs, zs = testimg.zs;
		SWorkImg<realnum>& yslice = m_pView->m_liftedEikonal.m_intey;
		yslice.Set(xs, zs);
		int ylev = m_cstartdir.GetPos();

		for (int zz = 0; zz < zs; ++zz) {
			for (int xx = 0; xx < xs; ++xx) {
				yslice[zz][xx] = testimg[zz][ylev][xx];
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

	//m_pView->m_transport.GetDispSlice(Talox, m_pView->m_xsee - 1, m_pView->m_dispd3); // new (Transport)


	SVoxImg<SWorkImg<realnum>> &testimg = m_pView->m_liftedEikonal.m_imageOp.m_testimage;
	/*
	SVoxImg<SWorkImg<realnum>> &testimg = m_pView->m_bsee?m_pView->m_liftedEikonal.m_imageOp.m_testinput:m_pView->m_liftedEikonal.m_imageOp.m_testimage;
	*/
	if (!m_pView->m_btransportview) {
		int ys = testimg.ys, zs = testimg.zs;
		SWorkImg<realnum>& xslice = m_pView->m_liftedEikonal.m_intex;
		xslice.Set(ys, zs);
		int xlev = m_cstopdir.GetPos();

		for (int zz = 0; zz < zs; ++zz) {
			for (int yy = 0; yy < ys; ++yy) {
				xslice[zz][yy] = testimg[zz][yy][xlev];
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
	//if (m_bini && !m_startstate) OnBnClickedButton2();
}

void CControlDlg::OnNMCustomdrawSlider4(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMCUSTOMDRAW pNMCD = reinterpret_cast<LPNMCUSTOMDRAW>(pNMHDR);
	// TODO: Add your control notification handler code here
	if (m_bini) {
		m_pView->m_liftedEikonal.m_expfac = m_cdatafac.GetPos();
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

void CControlDlg::OnBnClickedCheck1()
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

void CControlDlg::OnBnClickedButton6() // get minimal paths
{
	// TODO: Add your control notification handler code here
	//m_pView->GetAllPathX();
}

void CControlDlg::OnBnClickedCheck2()
{
	// TODO: Add your control notification handler code here
	m_pView->m_curvsee = m_ccursee.GetCheck();
	m_pView->Invalidate();
}

void CControlDlg::OnBnClickedButton7() // synth. image 2
{
	// TODO: Add your control notification handler code here
	if (m_bini) return;
	m_pView->m_liftedEikonal.m_imageOp.CreateTestImage2(XS_,YS_,ZS_);
	{
		m_pView->m_liftedEikonal.m_work = m_pView->m_liftedEikonal.m_imageOp.m_testimage[m_pView->m_zsee];
		m_pView->m_liftedEikonal.m_valid = 1;
		m_pView->m_liftedEikonal.m_grays = true;
		m_pView->m_liftedEikonal.m_color = false;
	}
	{
		m_zlevel.SetRange(0,ZS_-1,1); // from z
		m_cstartdir.SetRange(0,YS_-1,1); // from y
		m_cstopdir.SetRange(0,XS_-1,1); // from x
		m_bini = true;

		m_cdatafac.SetRange(1,10,1); m_cdatafac.SetPos(m_pView->m_liftedEikonal.m_expfac);
		char txt[22]; sprintf_s(txt,20,"%d",m_flowray);	m_cflowray.SetWindowTextA(txt);
	}
	m_pView->m_liftedEikonal.m_work.GetDispImg(m_pView->m_disp);
	OnNMCustomdrawSlider2(0,0);
	OnNMCustomdrawSlider3(0,0);
	m_pView->Invalidate();
}

void CControlDlg::OnBnClickedButton8() // correction using mean surf metric
{
	// TODO: Add your control notification handler code here
	extern bool g_modeswitch;
	g_modeswitch = !g_modeswitch;
	if (g_modeswitch)
		m_cmsw.SetWindowTextA("With corr");
	else
		m_cmsw.SetWindowTextA("W/O corr");
}


void CControlDlg::OnBnClickedCheck3()
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


void CControlDlg::OnBnClickedButton10()
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
