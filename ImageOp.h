#pragma once
#include "commontype.h"
#include <vector>
#include <unordered_set>
#include "Transport.h"

class CImageOp
{
public:
	CImageOp(void);
	~CImageOp(void);
	void CImageOp::GaussFilterImage(int nGauss);
	void CalcOneComponent(realnum r, SWorkImg<realnum> &comp);
	void CalculateCircularGradientFlow(realnum r, SWorkImg<realnum> &img);
	int GetDirectionalFlow(realnum r, CVec2 dir, SWorkImg<realnum> &img, SWorkImg<realnum> &out, bool bForced = false, int expcoef = 5, bool bReverse = false, int nG = 3);

	SWorkImg<realnum> m_gxx;
	SWorkImg<realnum> m_gxy;
	SWorkImg<realnum> m_gyy;
	SWorkImg<realnum> m_aux;
	SWorkImg<realnum> m_loc;

	std::unordered_set<unsigned long> m_bound;
	void GetXTestDisTrans();
	void GetXTestIntern();
	void XTestfill(short x, short y);
	static const int m_bandthick = 25;
public:
	SVoxImg<SWorkImg<realnum>>& CImageOp::GetIniMap(int ix);
	void GetXTestBound(int ix, std::vector<CVec3>& out);
	void GetPlaneDistMap(std::unordered_set<unsigned long> &boundset);
	SWorkImg<realnum> &GetXTestSigDis()
	{
		return m_loc;
	}
	SVoxImg<SWorkImg<realnum>>& GetTestImage() {
		return m_testimage;
	}
	SVoxImg<SWorkImg<realnum>>& GetTestInput() {
		return m_testinput;
	}
	SVoxImg<SWorkImg<realnum>>& CreateTestImage(int xs, int ys, int zs, int expcoef = 9+5); // 9,19
	SVoxImg<SWorkImg<realnum>>& CreateTestImage2(int xs, int ys, int zs, int expcoef = 9); // 9,19
	SVoxImg<SWorkImg<realnum>> m_testimage;
	SVoxImg<SWorkImg<realnum>> m_testinput;
	SVoxImg<SWorkImg<realnum>> m_inimap;
	SWorkImg<realnum>& GetXImageSlice(int ix);
private:
	void GauTest(bool bodd);

	bool m_bOk;
};
