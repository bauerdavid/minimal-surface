
// ChildView.h : interface of the CChildView class
//


#pragma once

#include "AreaEikonal.h"
#include "Transport.h"

#include "afxwin.h"
#include "afxcmn.h"

#define EXP_COEF_DEF 9

class CControlDlg;

// CChildView window

class CChildView : public CWnd
{
	bool pressed = false;
// Construction
public:
	CChildView();

// Attributes

public:
	CControlDlg *m_pControl;
	CCurvEikonal m_liftedEikonal;
	SDisImg m_disp;
	SDisImg m_dispd1;
	SDisImg m_dispd2;
	SDisImg m_dispd3;
	SDisImg m_dislic;

	CTransport m_transport;
	CImageOp m_imageOp;

	void DispImage(CPaintDC &dc, SDisImg &disp, int offx, int offy);

	void InitThread();
	void PauseThread(int threadstat);
	void StopThread();
	int m_prevthreadactivated;
	int m_threadactivated;

	CVec3 m_start_point;
	CVec3 m_end_point;
	
	int m_bispoint;
	/*int m_fieldliney;
	CPoint m_setcenter;
	int m_setray;
	bool m_bset;*/
	int m_zsee;
	int m_ysee;
	int m_xsee;
	void Z2Dir(int z, int &dx, int &dy);

	bool m_bsee;
	int m_curvsee;
	bool m_btransportview;
	

	// image general
	int m_expfac;
	int m_valid;
	bool m_grays;
	bool m_color;

	SDisImg m_img;
	SWorkImg<realnum> m_work;
	SWorkImg<realnum> m_workr;
	SWorkImg<realnum> m_workg;
	SWorkImg<realnum> m_workb;

	SWorkImg<realnum> m_intex;
	SWorkImg<realnum> m_intey;
// Operations
public:

// Overrides
	protected:
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);

// Implementation
public:
	virtual ~CChildView();

	// Generated message map functions
protected:
	afx_msg void OnPaint();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
};

#pragma once


// CControlDlg dialog

class CControlDlg : public CDialog
{
	DECLARE_DYNAMIC(CControlDlg)

public:
	CControlDlg(CWnd* pParent = NULL);   // standard constructor
	virtual ~CControlDlg();

	CChildView *m_pView;
	int m_startstate;
	int m_fromstate;
	int m_flowray;
	bool m_bini;
	bool m_bforced;

// Dialog Data
	enum { IDD = IDD_DIALOG1 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedOpen();
	afx_msg void OnBnClickedIntImage1();
	afx_msg void OnBnClickedGetBC();
	afx_msg void OnBnClickedStartStop();
	CButton m_startstop;
	afx_msg void OnBnClickedSuspend();
	CEdit m_cdist;
	CSliderCtrl m_zlevel;
	afx_msg void OnNMCustomdrawSlider1(NMHDR *pNMHDR, LRESULT *pResult);
	CEdit m_ezlevel;
	CSliderCtrl m_cstartdir;
	CSliderCtrl m_cstopdir;
	afx_msg void OnNMCustomdrawSlider2(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnNMCustomdrawSlider3(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnEnChangeEdit3();
	CEdit m_cflowray;
	CSliderCtrl m_cdatafac;
	afx_msg void OnNMCustomdrawSlider4(NMHDR *pNMHDR, LRESULT *pResult);
	afx_msg void OnEnChangeEdit2();
	afx_msg void OnBnClickedCheckSee();
	afx_msg void OnEnChangeEdit1();
	afx_msg void OnBnClickedButton6();
	afx_msg void OnBnClickedCSee();
	CButton m_ccursee;
	afx_msg void OnBnClickedIntImage2();
	afx_msg void OnBnClickedWWOCorr();
	CButton m_cmsw;
	afx_msg void OnBnClickedSwitchView();
	afx_msg void OnBnClickedTestBuild();
};
