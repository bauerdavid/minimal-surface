#pragma once
#include "commontype.h"
#define _NO_BOU_ 1e11
#include "SimpleITK.h"
namespace sitk = itk::simple;

enum {Talox, Taloy, Taloz};

class CTransport
{
public:
	CTransport():m_active(0) {}

	int m_active;
	realnum m_min;
	// Initialize gradients fro, distance map, transform function from inimap, and boundaries (considering inimap, and setting data edges as boundaries)
	void TrInit(sitk::Image& distmap, sitk::Image& inimap, realnum maxdistance);
	void TrControl(int nIter);
	void GetDispSlice(int along, int at, SDisImg& r);
	sitk::Image m_transportfunction[2];
	//SVoxImg<SWorkImg<double>> m_transportfunction[2];
	/*
	Boundary state
	 1: read-only
	 0: calculation is needed
	 -1: not initialized
	*/
	//SVoxImg<SWorkImg<int>> m_isboundary;
	sitk::Image m_isboundary;
	sitk::Image m_gx;
	sitk::Image m_gy;
	sitk::Image m_gz;
private:
	bool TrIterate(int bev);
	



};

