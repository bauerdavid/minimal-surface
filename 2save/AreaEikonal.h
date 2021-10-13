#pragma once

#include "commontype.h"
#include <vector>

struct IPoi {
	int x;
	int y;
	IPoi():x(0),y(0) {}
	IPoi(int xx, int yy):x(xx),y(yy) {}
};

class CPhaseContainer
{
public:
	// Resolve path
	SWorkImg<realnum> m_distance;
	SWorkImg<realnum> m_distgrad[2];
	// Working field
	SWorkImg<realnum> m_field;
	SWorkImg<realnum> m_gradlen;
	SWorkImg<realnum> m_n[2];
	//SDisImg m_fieldgrad[2];
	SWorkImg<realnum> m_aux;

	SWorkImg<realnum> m_velo;
};

class CAreaEikonal
{
public:
	CAreaEikonal(void);
	~CAreaEikonal(void);

	// phase field stuff
	CPhaseContainer m_phasefield;
	void PhaseInit(IPoi reginit, IPoi arrival);
	void RegularizePhaseField(SWorkImg<realnum> &field, SWorkImg<realnum> &velo);
	void CalculateFundQuant();
	void Iterate();
	realnum m_currentdistance;
	CVec2 m_reference;
	void ResolvePath(realnum x,realnum y, bool bClear = true);
	std::vector<CVec2> m_minpath;
	/**/
	int m_resolvready;
	int m_inittype;
	IPoi m_distancefrom;
	IPoi m_distanceto;
	/**/

	// image general
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
	SWorkImg<realnum> m_divcheck;

	void GetIntImage(SWorkImg<realnum> &work, SWorkImg<realnum> &vx, SWorkImg<realnum> &vy);
	void GetDivergence(SWorkImg<realnum> &vx, SWorkImg<realnum> &vy, SWorkImg<realnum> &div);
};
