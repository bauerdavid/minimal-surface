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
	
	SWorkImg<realnum> m_gxx;
	SWorkImg<realnum> m_gxy;
	SWorkImg<realnum> m_gyy;
	SWorkImg<realnum> m_aux;
	SWorkImg<realnum> m_loc;

	std::unordered_set<unsigned long> m_bound;
	// Replace object mask with distances from the bounds
	void GetXTestDisTrans();
	// fill m_loc with ones inside the bound (including the bounding pixels), and set to zero everywhere else
	void XTestfill(short x, short y);
	static const int m_bandthick = 25;
public:
	/*
	Create data, where every values is _NO_BOU_, except for the selected x slice, where it will be the same as m_loc
	Should be initialized based on the meeting of the two flow
	*/
	SVoxImg<SWorkImg<realnum>>& CImageOp::GetIniMap(int ix);
	void GetXTestBound(int ix, std::vector<CVec3>& out);
	// Calculate the distance from the bounds. Distances will be positive inside the object, and negative outside of it.
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
