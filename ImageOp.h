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
	sitk::Image& CImageOp::GetIniMap(int xs, int ys, int zs, int ix);
	void GetXTestBound(int ix, std::vector<CVec3>& out);
	// Calculate the distance from the bounds. Distances will be positive inside the object, and negative outside of it.
	void GetPlaneDistMap(int distmap_ys, int distmap_zs, std::unordered_set<unsigned long> &boundset);
	SWorkImg<realnum> &GetXTestSigDis()
	{
		return m_loc;
	}
	sitk::Image& GetTestImage() {
		return m_testimage;
	}
	sitk::Image& GetTestInput() {
		return m_testinput;
	}
	sitk::Image& CreateTestImage(int xs, int ys, int zs, int expcoef = 9+5); // 9,19
	sitk::Image& CreateTestImage2(int xs, int ys, int zs, int expcoef = 9); // 9,19
	sitk::Image& CreateTestInput(double expcoef=14);
	sitk::Image m_testimage;
	sitk::Image m_testinput;
	sitk::Image m_inimap;
	SWorkImg<realnum>& GetXImageSlice(int ix);
	void GauTest(bool bodd);
private:

	bool m_bOk;
};
