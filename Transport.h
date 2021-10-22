#pragma once
#include "commontype.h"
#define _NO_BOU_ 1e11

enum {Talox, Taloy, Taloz};

class CTransport
{
public:
	CTransport():m_active(0) {}

	int m_active;
	realnum m_min;
	// Initialize gradients fro, distance map, transform function from inimap, and boundaries (considering inimap, and setting data edges as boundaries)
	void TrInit(SVoxImg<SWorkImg<realnum>>& distmap, SVoxImg<SWorkImg<realnum>>& inimap, realnum maxdistance);
	void TrControl(int nIter);
	void GetDispSlice(int along, int at, SDisImg& r);
	SVoxImg<SWorkImg<realnum>> m_transportfunction[2];
	/*
	Boundary state
	 1: read-only
	 0: calculation is needed
	 -1: not initialized
	*/
	SVoxImg<SWorkImg<int>> m_isboundary;
private:
	void TrIterate(int bev);
	SVoxImg<SWorkImg<realnum>> m_gx;
	SVoxImg<SWorkImg<realnum>> m_gy;
	SVoxImg<SWorkImg<realnum>> m_gz;



};

