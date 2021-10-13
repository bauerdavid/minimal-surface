#pragma once
#include "StdAfx.h"
#include "math.h"
typedef double realnum;
#define PI 3.1415926536


struct SWorkData
{
	SWorkData():xs(0),ys(0) {
		data = 0;
	}
	void WorkInit(int x, int y, int inres = 1) {
		if (data) {
			if (xs != x || ys != y) {
				WorkFree();
			}
		}
		xs = x;
		ys = y;
		res = inres;
		del = realnum(res);
		idel = 1.0f/del;
		if (!data) {
			data = new realnum[xs*ys];
		}
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				data[yy*xs+xx] = 0;
			}
		}

	}	
	SWorkData(int x, int y, int r = 1) {
		WorkInit(x,y,r);
	}

	void WorkFree() {
		delete[] data;
		data = 0; xs = ys = 0;
	}
	SWorkData &operator=(SWorkData &s) {
		if (&s == this) return *this;
		WorkInit(s.xs, s.ys, s.res);
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				data[yy/xs+xx] = s[yy][xx];
			}
		}
		return *this;
	}
	void GetLaplace(SWorkData &s) {
		for (int yy = 1; yy < ys-1; ++yy) {
			for (int xx = 1; xx < xs-1; ++xx) {
				data[yy*xs+xx] = -4*s[yy][xx]+s[yy][xx+1]+s[yy][xx-1]+s[yy+1][xx]+s[yy-1][xx];
			}
		}
		int y = ys-1;
		for (int xx = 1; xx < xs-1; ++xx) {
			data[xx] = -4*s[0][xx]+s[0][xx+1]+s[0][xx-1]+s[1][xx]+s[1][xx];
			data[y*xs+xx] = -4*s[y][xx]+s[y][xx+1]+s[y][xx-1]+s[y-1][xx]+s[y-1][xx];
		}
		int x = xs-1;
		for (int yy = 1; yy < ys-1; ++yy) {
			data[yy*xs+0] = -4*s[yy][0]+s[yy+1][0]+s[yy-1][0]+s[yy][1]+s[yy][1];
			data[yy*xs+x] = -4*s[yy][x]+s[yy+1][x]+s[yy-1][x]+s[yy][x-1]+s[yy][x-1];
		}
		data[0*xs+0] = 0.5f*(data[0*xs+1]+data[1*xs+0]);
		data[0*xs+x] = 0.5f*(data[0*xs+x-1]+data[1*xs+x]);
		data[y*xs+0] = 0.5f*(data[(y-1)*xs+0]+data[y*xs+1]);
		data[y*xs+x] = 0.5f*(data[(y-1)*xs+x]+data[y*xs+x-1]);

		*this *= idel*idel;
		//return *this;
	}
	void GetHesseXY(SWorkData &s) {
		for (int yy = 1; yy < ys-1; ++yy) {
			for (int xx = 1; xx < xs-1; ++xx) {
				data[yy*xs+xx] = 0.25f*(s[yy+1][xx+1]+s[yy-1][xx-1]-s[yy+1][xx-1]-s[yy-1][xx+1]);
			}
		}
		int y = ys-1;
		for (int xx = 1; xx < xs-1; ++xx) {
			data[xx] = 0;
			data[y*xs+xx] = 0;
		}
		int x = xs-1;
		for (int yy = 1; yy < ys-1; ++yy) {
			data[yy*xs+0] = 0;
			data[yy*xs+x] = 0;
		}
		data[0*xs+0] = 0.5f*(data[0*xs+1]+data[1*xs+0]);
		data[0*xs+x] = 0.5f*(data[0*xs+x-1]+data[1*xs+x]);
		data[y*xs+0] = 0.5f*(data[(y-1)*xs+0]+data[y*xs+1]);
		data[y*xs+x] = 0.5f*(data[(y-1)*xs+x]+data[y*xs+x-1]);

		*this *= idel*idel;
		//return *this;
	}
	void GetHesseXX(SWorkData &s) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 1; xx < xs-1; ++xx) {
				data[yy*xs+xx] = s[yy][xx+1]+s[yy][xx-1]-2*s[yy][xx];
			}
		}
		int x = xs-1;
		for (int yy = 0; yy < ys; ++yy) {
			data[yy*xs+0] = 0;
			data[yy*xs+x] = 0;
		}

		*this *= idel*idel;
		//return *this;
	}
	void GetHesseYY(SWorkData &s) {
		for (int yy = 1; yy < ys-1; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				data[yy*xs+xx] = s[yy+1][xx]+s[yy-1][xx]-2*s[yy][xx];
			}
		}
		int y = ys-1;
		for (int xx = 0; xx < xs; ++xx) {
			data[0*xs+xx] = 0;
			data[y*xs+xx] = 0;
		}

		*this *= idel*idel;
		//return *this;
	}
	void GetUpGradX(SWorkData &s) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 1; xx < xs-1; ++xx) {
				data[yy*xs+xx] = s[yy][xx+1]-s[yy][xx];
			}
		}
		int x = xs-1;
		for (int yy = 0; yy < ys; ++yy) {
			data[yy*xs+0] = s[yy][1]-s[yy][0];
			data[yy*xs+x] = s[yy][x]-s[yy][x-1];
		}
	}
	void GetUpGradY(SWorkData &s) {
		for (int yy = 1; yy < ys-1; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				data[yy*xs+xx] = s[yy+1][xx]-s[yy][xx];
			}
		}
		int y = ys-1;
		for (int xx = 0; xx < xs; ++xx) {
			data[0*xs+xx] = s[1][xx]-s[0][xx];
			data[y*xs+xx] = s[y][xx]-s[y-1][xx];
		}
	}
	void GetGradX(SWorkData &s) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 1; xx < xs-1; ++xx) {
				data[yy*xs+xx] = 0.5f*(s[yy][xx+1]-s[yy][xx-1]);
			}
		}
		int x = xs-1;
		for (int yy = 0; yy < ys; ++yy) {
			data[yy*xs+0] = s[yy][1]-s[yy][0];
			data[yy*xs+x] = s[yy][x]-s[yy][x-1];
		}
		//*this *= idel;
		//return *this;
	}
	void GetGradY(SWorkData &s) {
		for (int yy = 1; yy < ys-1; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				data[yy*xs+xx] = 0.5f*(s[yy+1][xx]-s[yy-1][xx]);
			}
		}
		int y = ys-1;
		for (int xx = 0; xx < xs; ++xx) {
			data[0*xs+xx] = s[1][xx]-s[0][xx];
			data[y*xs+xx] = s[y][xx]-s[y-1][xx];
		}
		//*this *= idel;
		//return *this;
	}
	SWorkData& operator*= (realnum f) {
		for (int yy = 0; yy < ys; ++yy) {
			for (int xx = 0; xx < xs; ++xx) {
				data[yy*xs+xx] *= f;
			}
		}
		return *this;
	}

	realnum *operator[](int y) {
		return &data[y*xs];
	}

	~SWorkData() {
		WorkFree();
	}

	int xs;
	int ys;
	int res;
	realnum del;
	realnum idel;
	realnum *data;
};

class CVec2  
{
public:
	CVec2() {};
	CVec2(double x, double y) {
		v[0] = x;
		v[1] = y;
	};
	~CVec2() {};
	
	union {
		double v[2];
		struct {
			double x,y;
		};
	};
	
	double &operator[](int i) {
		return v[i];
	}
	CVec2 &operator*=(double s) {
		v[0] *= s;
		v[1] *= s;
		return *this;
	}
	CVec2 &operator+=(CVec2 &w) {
		v[0] += w[0];
		v[1] += w[1];
		return *this;
	}
	CVec2 &operator-=(CVec2 &w) {
		v[0] -= w[0];
		v[1] -= w[1];
		return *this;
	}
	
	CVec2 &Norm() {
		double l = 1.0/sqrt(v[0]*v[0]+v[1]*v[1]);
		v[0] *= l;
		v[1] *= l;
		return *this;
	}

	double Len() {
		return sqrt(v[0]*v[0]+v[1]*v[1]);
	}

};
CVec2 operator *(double f, CVec2 &s);
CVec2 operator *(int f, CVec2 &s);
realnum operator *(CVec2 &s1, CVec2 &s2);

class CVec3  
{
public:
	CVec3() {
		v[0] = v[1] = v[2] = 0;
	};
	CVec3(double x, double y, double z) {
		v[0] = x;
		v[1] = y;
		v[2] = z;
	};
	~CVec3() {};
	
	union {
		double v[3];
		struct {
			double x,y,z;
		};
	};
	
	double &operator[](int i) {
		return v[i];
	}
	CVec3 &operator*=(double s) {
		v[0] *= s;
		v[1] *= s;
		v[2] *= s;
		return *this;
	}
	CVec3 &operator+=(CVec3 &w) {
		v[0] += w[0];
		v[1] += w[1];
		v[2] += w[2];
		return *this;
	}
	CVec3 &operator-=(CVec3 &w) {
		v[0] -= w[0];
		v[1] -= w[1];
		v[2] -= w[2];
		return *this;
	}
	
	CVec3 &Norm() {
		double l = 1.0/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[3]);
		v[0] *= l;
		v[1] *= l;
		v[2] *= l;
		return *this;
	}

	double Len() {
		return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[3]);
	}

};

CVec3 operator *(double f, CVec3 &v);

class CMat3  
{
public:
	CMat3() {
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				m[i][j] = 0;
	};

	~CMat3() {};
	
	union {
		double m[3][3];
		struct {
			CVec3 r[3];
		};
	};

	void SetId() {
		m[0][0] = 1; m[0][1] = 0; m[0][2] = 0;
		m[1][0] = 0; m[1][1] = 1; m[1][2] = 0;
		m[2][0] = 0; m[2][1] = 0; m[2][2] = 1;
	}
	
	CVec3 &operator[](int i) {
		return r[i];
	}
	CMat3 &operator=(const CMat3 &s) {
		m[0][0] = s.m[0][0]; m[0][1] = s.m[0][1]; m[0][2] = s.m[0][2];
		m[1][0] = s.m[1][0]; m[1][1] = s.m[1][1]; m[1][2] = s.m[1][2];
		m[2][0] = s.m[2][0]; m[2][1] = s.m[2][1]; m[2][2] = s.m[2][2];
		return *this;
	}
	CMat3 &operator*=(double s) {
		r[0] *= s;
		r[1] *= s;
		r[2] *= s;
		return *this;
	}
	CMat3 &operator+=(CMat3 &w) {
		r[0] += w[0];
		r[1] += w[1];
		r[2] += w[2];
		return *this;
	}
	CMat3 &operator-=(CMat3 &w) {
		r[0] -= w[0];
		r[1] -= w[1];
		r[2] -= w[2];
		return *this;
	}
	
	CMat3 &operator*=(CMat3 &w)  {
		CMat3 tm;
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				for (int k = 0; k < 3; ++k)
					tm[i][j] += m[i][k]*w[k][j];
		/*for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				m[i][j] = tm[i][j];*/
		*this = tm;
		return *this;
	}
	CMat3 &Transpose()  {
		double a;
		a = m[0][1]; m[0][1] = m[1][0]; m[1][0] = a;
		a = m[0][2]; m[0][2] = m[2][0]; m[2][0] = a;
		a = m[1][2]; m[1][2] = m[2][1]; m[2][1] = a;
		return *this;
	}

};

/*CMat3 operator *(double f, CMat3 &m);
CMat3 operator *(CMat3 &m1, CMat3 &m2);
CVec3 operator *(CMat3 &m, CVec3 &v);
CVec3 operator *(CVec3 &v, CMat3 &m);*/

///////////////////////////////////////////////////////////////////////////////
#define HFRAME 8

struct SDisImg {
	SDisImg() {
		dat = 0;
		xs = ys = 0;
	}
	SDisImg(int x, int y, byte *buf, int mod = 0, int pitch = 0) {
		dat = 0;
		xs = ys = 0;
		Set(x,y,buf,mod, pitch);
	}
	void Reset(int x, int y) {
		if (dat && (x != xs || y != ys))
			Clean();
		xs = x;
		ys = y;
		if (!dat && xs && ys) {
			dat = new unsigned long[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return;
			}
		}
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] = 0;
			}
		}
	}
	void Set(int x, int y, byte *buf, int mod = 0, int pitch = 0) {
		if (dat && (x != xs || y != ys))
			Clean();
		xs = x;
		ys = y;
		if (!dat && xs && ys) {
			dat = new unsigned long[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return;
			}
		}
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				byte r, g, b;
				r = g = b = *buf++;
				dat[q*xs+p] = (r<<16)+(g<<8)+b;
			}
			buf += mod;
			if (pitch < 0) buf += pitch*2;
		}
	}
	// kaggle
	void InitKaggleSaveImg(int x, int y) {
		xs = x;
		ys = y;
		dat = new unsigned long[ys*xs];
		if (!dat) {
			xs = ys = 0;
			return;
		}
		for (int q = 0; q < ys; ++q) 
			for (int p = 0; p < xs; ++p)
				dat[q*ys+p] = 0;
	}
	void Set32Kaggle(int x, int y, byte *buf, int mod = 0, int pitch = 0) {

		if (dat && (x != xs || y != ys))
			Clean();
		xs = x+HFRAME*2;
		ys = y+HFRAME*2;
		if (!dat && xs && ys) {
			dat = new unsigned long[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return;
			}
		}
		int avg = 0;
		for (int q = 0+HFRAME; q < ys-HFRAME; ++q) {
			for (int p = 0+HFRAME; p < xs-HFRAME; ++p) {
				realnum gray;
				byte r, g, b;
				r = *buf++; 
				g = *buf++;
				b = *buf++;
				gray = 0.3f*r+0.5f*g+0.2f*b; 
				if (gray > 255) gray = 255;
				r = g = b = byte(gray);
				dat[q*xs+p] = (r<<16)+(g<<8)+b;
				avg += b;
				buf++;
			}
			buf += mod;
			if (pitch < 0) buf += pitch*2;
		}
		avg /= (xs-HFRAME*2)*(ys-HFRAME*2);
		int hist[256];
		for (int i = 0; i < 256; ++i) hist[i] = 0;
		for (int q = 0+HFRAME; q < ys-HFRAME; ++q)
			for (int p = 0+HFRAME; p < xs-HFRAME; ++p)
				hist[dat[q*xs+p]&0xff]++;
		int max = 0; byte mi = 0;
		for (int i = 0; i < 256; ++i) {
			if (hist[i] > max) { max = hist[i]; mi = i; }
		}
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				if (p >= HFRAME && p < xs-HFRAME && q >= HFRAME && q < ys-HFRAME) continue;
				dat[q*xs+p] = mi;
			}
		}
		if (mi > avg) {
			for (int q = 0; q < ys; ++q)
				for (int p = 0; p < xs; ++p)
					dat[q*xs+p] = 255-dat[q*xs+p];
		}


	}
	// kaggle
	void SetColor(int x, int y, byte *buf, int mod = 0, int pitch = 0) {
		if (dat && (x != xs || y != ys))
			Clean();
		xs = x;
		ys = y;
		if (!dat && xs && ys) {
			dat = new unsigned long[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return;
			}
		}
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				byte r, g, b;
				r = *buf++;
				g = *buf++;
				b = *buf++;
				dat[q*xs+p] = (r<<16)+(g<<8)+b;
			}
			buf += mod;
			if (pitch < 0) buf += pitch*2;
		}
	}
	void Clean() {
		if (dat) {
			delete[] dat;
			dat = 0;
			xs = ys = 0;
		}
	}
	~SDisImg() {
		Clean();
	}
	unsigned long* operator[](int y) {
		if (y >= ys) y = ys-1;
		else if (y < 0) y = 0;
		//if (!dat) return &safe[1][1];
		return &dat[y*xs];
	}
	void operator=(SDisImg &tc) {
		if (&tc == this) return;
		if (xs != tc.xs || ys != tc.ys) {
			Clean();
			xs = tc.xs;
			ys = tc.ys;
			dat = new unsigned long[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return;
			}
		}

		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] = tc[q][p];
			}
		}

	}
	unsigned long *dat;
	int xs;
	int ys;
};

struct SSavImg {
	SSavImg() {
		dat = 0;
		xs = ys = 0;
	}
	~SSavImg() {
		if (dat) {
			delete[] dat;
			dat = 0;
		}
	}

	unsigned char *dat;
	int xs;
	int ys;
};

// working with real numbers

template <class T> struct SWorkImg {
	SWorkImg() {
		dat = 0;
		xs = ys = 0;
	}
	SWorkImg(int x, int y) {
		dat = 0;
		xs = ys = 0;
		Set(x,y);
	}
	SWorkImg(int x, int y, byte *buf) {
		dat = 0;
		xs = ys = 0;
		Set(x,y,buf);
	}
	void Set(int x, int y) {
		if (dat && (x != xs || y != ys))
			Clean();
		xs = x;
		ys = y;
		if (!dat && xs && ys) {
			dat = new T[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return;
			}
		}
		maxval = minval = avgval = 0;
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] = (T)0;
			}
		}
	}
	void Set(int x, int y, T fill) {
		if (dat && (x != xs || y != ys))
			Clean();
		xs = x;
		ys = y;
		if (!dat && xs && ys) {
			dat = new T[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return;
			}
		}
		maxval = minval = avgval = 0;
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] = fill;
			}
		}
	}
	void Set(int x, int y, byte *buf) {
		if (dat && (x != xs || y != ys))
			Clean();
		xs = x;
		ys = y;
		if (!dat && xs && ys) {
			dat = new T[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return;
			}
		}
		maxval = 0;
		minval = (T)10000;
		avgval = 0;
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				T t = (T)*buf++;
				dat[q*xs+p] = t;
				if (maxval < t) maxval = t;
				if (minval > t) minval = t;
				avgval += t;
			}
		}
		avgval /= xs*ys;
		Norm();
	}
	void SetBound() {
		if (!xs || !ys || !dat) return;
		maxval = (T)-1e11;
		minval = (T)1e11;
		avgval = 0;
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				T t = dat[q*xs+p];
				if (maxval < t) maxval = t;
				if (minval > t) minval = t;
				avgval += t;
			}
		}
		avgval /= xs*ys;
	}
	void Clean() {
		if (dat) {
			delete[] dat;
			dat = 0;
			xs = ys = 0;
		}
	}
	~SWorkImg() {
		Clean();
	}
	T* operator[] (int y) {
		if (y >= ys) y = ys-1;
		else if (y < 0) y = 0;
		//if (!dat) return &safe[1][1];
		return &dat[y*xs];
	}
	const T* operator[](int y) const {
		int iy = y;
		if (iy >= ys) iy = ys-1;
		else if (iy < 0) iy = 0;
		return &dat[iy*xs];
	}
	SWorkImg& operator= (SWorkImg &tc) {
		if (&tc == this) return *this;
		if (xs != tc.xs || ys != tc.ys) {
			Clean();
			xs = tc.xs;
			ys = tc.ys;
			dat = new T[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return *this;
			}
		}

		maxval = tc.maxval;
		minval = tc.minval;
		avgval = tc.avgval;
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] = tc[q][p];
			}
		}
		return *this;
	}
	void operator=(SDisImg &tc) {
		if (xs != tc.xs || ys != tc.ys) {
			Clean();
			xs = tc.xs;
			ys = tc.ys;
			dat = new T[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return;
			}
		}

		maxval = 0;
		minval = (T)10000;
		avgval = 0;
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				T t;
				t = (T)(tc[q][p]&0xff);

				dat[q*xs+p] = t;
				if (maxval < t) maxval = t;
				if (minval > t) minval = t;
				avgval += t;
			}
		}
		avgval /= xs*ys;
		Norm();
	}
	bool GetAligned(SDisImg &tc, int align, int channel = 0) {
		int tcxs = tc.xs, tcys = tc.ys;
		int modxs = tcxs%align, modys = tcys%align;
		tcxs -= modxs; tcys -= modys;
		if (tcxs < align || tcys < align) return false;
		if (xs != tcxs || ys != tcys) {
			Clean();
			xs = tcxs;
			ys = tcys;
			dat = new T[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return false;
			}
		}

		maxval = 0;
		minval = (T)10000;
		avgval = 0;
		modxs >>= 1; modys >>= 1;
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				T t;
				if (!channel)
					t = (T)(tc[q+modys][p+modxs]&0xff);
				else if (channel == 1) 
					t = (T)((tc[q+modys][p+modxs]&0xff00)>>8);
				else 
					t = (T)((tc[q+modys][p+modxs]&0xff0000)>>16);
				dat[q*xs+p] = t;
				if (maxval < t) maxval = t;
				if (minval > t) minval = t;
				avgval += t;
			}
		}
		avgval /= xs*ys;
		Norm();
		return true;
	}
	bool GetAligned8(SDisImg &tc, int align) {
		int tcxs = tc.xs, tcys = tc.ys;
		int modxs = tcxs%align, modys = tcys%align;
		tcxs -= modxs; tcys -= modys;
		if (tcxs < align || tcys < align) return false;
		if (xs != tcxs || ys != tcys) {
			Clean();
			xs = tcxs;
			ys = tcys;
			dat = new T[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return false;
			}
		}

		maxval = 0;
		minval = (T)10000;
		avgval = 0;
		modxs >>= 1; modys >>= 1;
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] = (T)(tc[q+modys][p+modxs]);
				if (maxval < t) maxval = t;
				if (minval > t) minval = t;
				avgval += t;
			}
		}
		avgval /= xs*ys;
		Norm();
		return true;
	}
	// kaggle
	bool GetRenormed(SDisImg &tc) {
		int tcxs = tc.xs, tcys = tc.ys;
		if (tcxs < HFRAME || tcys < HFRAME) return false;
		if (xs != tcxs || ys != tcys) {
			Clean();
			xs = tcxs;
			ys = tcys;
			dat = new T[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return false;
			}
		}

		maxval = 0;
		minval = (T)10000;
		avgval = 0;
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				T t;
				t = (T)(tc[q][p]&0xff);
				dat[q*xs+p] = t;
				if (maxval < t) maxval = t;
				if (minval > t) minval = t;
				avgval += t;
			}
		}
		avgval /= xs*ys;
		Renorm();
		return true;
	}
	// kaggle
	void Norm()
	{
		if (maxval <= 0) return;
		T recip = ((T)1.0)/((T)maxval);
		//T recip = ((T)1.0)/((T)255);
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] *= recip;
			}
		}
	}
	void Renorm()
	{
		T min = (T)1e11f, max = (T)-1e11f;
		for (int q = 0; q < ys; ++q) 
			for (int p = 0; p < xs; ++p) {
				T v = dat[q*xs+p];
				if (v < min) min = v;
				if (v > max) max = v;
			}

		if (abs((T)(max-min)) < ((T)1e-11f)) return;

		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] -= (T)min;
				dat[q*xs+p] /= (T)(max-min);		
			}
		}
		SetBound();

	}
	void Invert()
	{
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] = ((T)1.0)-dat[q*xs+p];
			}
		}
	}
	void GetReduced2(SWorkImg &r)
	{
		int xx = xs/2, yy = ys/2;

		if (!xx || !yy) return;
		if (xx != r.xs || yy != r.ys) {
			r.Clean();
			r.dat = new T[yy*xx];
			if (!r.dat) return;
			r.xs = xx;
			r.ys = yy;
		}

		r.maxval = 0;
		r.minval = (T)10000;
		r.avgval = 0;
		for (int q = 0; q < yy; ++q) {
			for (int p = 0; p < xx; ++p) {
				int qq = 2*q, pp = 2*p;
				T t = dat[qq*xs+pp];
				t += dat[qq*xs+pp+1];
				t += dat[(qq+1)*xs+pp];
				t += dat[(qq+1)*xs+pp+1];
				t *= ((T)0.25);
				r.dat[q*r.xs+p] = t;
				if (r.maxval < t) r.maxval = t;
				if (r.minval > t) r.minval = t;
				r.avgval += t;
			}
		}
		r.avgval /= r.xs*r.ys;
		r.Norm();
		return;
	}
	void GetAugmented2(SWorkImg &r)
	{
		int xx = xs*2, yy = ys*2;

		if (!xx || !yy) return;
		if (xx != r.xs || yy != r.ys) {
			r.Clean();
			r.dat = new T[yy*xx];
			if (!r.dat) return;
			r.xs = xx;
			r.ys = yy;
		}

		r.maxval = 0;
		r.minval = (T)10000;
		r.avgval = 0;
		for (int q = 0; q < ys; ++q) {
			int qq = q*2;			
			int q1 = q+1;
			if (q1 >= ys) q1 = ys-1;
			for (int p = 0; p < xs; ++p) {
				int pp = p*2;
				int p1 = p+1;
				if (p1 >= xs) p1 = xs-1;
				T t00 = dat[q*xs+p];
				T t01 = dat[q*xs+p1];
				T t10 = dat[q1*xs+p];
				T t11 = dat[q1*xs+p1];
				r.dat[qq*r.xs+pp] = t00;
				r.dat[qq*r.xs+pp+1] = ((T)0.5)*(t00+t01);
				r.dat[(qq+1)*r.xs+pp] = ((T)0.5)*(t00+t10);
				r.dat[(qq+1)*r.xs+pp+1] = ((T)0.25)*(t00+t01+t10+t11);
			}
		}
		r.avgval = avgval;
		r.maxval = maxval;
		r.minval = minval;
	}
	void GetDispChannel(SDisImg &r, int channel = 0)
	{
		if (!xs || !ys) return;
		if (!channel) {
			if (xs != r.xs || ys != r.ys) {
				r.Clean();
				r.dat = new unsigned long[ys*xs];
				if (!r.dat) return;
				r.xs = xs;
				r.ys = ys;
			}
		}

		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				T t = dat[q*xs+p];
				t *= 0xff;
				if (t > (T)0xff) t = (T)0xff;
				byte b;
				b = (byte)t;
				if (!channel) 
					r.dat[q*xs+p] = b;
				else if (channel == 1)
					r.dat[q*xs+p] += b<<8;
				else
					r.dat[q*xs+p] += b<<16;
			}
		}

	}
	void GetDispImg(SDisImg &r)
	{
		if (!xs || !ys) return;
		if (xs != r.xs || ys != r.ys) {
			r.Clean();
			r.dat = new unsigned long[ys*xs];
			if (!r.dat) return;
			r.xs = xs;
			r.ys = ys;
		}

		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				T t = dat[q*xs+p];
				t *= 0xff;
				if (t > (T)0xff) t = (T)0xff;
				byte R, g, b;
				R = g = b = (byte)t;
				r.dat[q*xs+p] = (R<<16)+(g<<8)+b;
			}
		}

	}
	void GetNormalizedDispImg(SDisImg &r)
	{
		if (!xs || !ys) return;
		if (xs != r.xs || ys != r.ys) {
			r.Clean();
			r.dat = new unsigned long[ys*xs];
			if (!r.dat) return;
			r.xs = xs;
			r.ys = ys;
		}

		T min = (T)1e11f, max = (T)-1e11f;
		for (int q = 0; q < ys; ++q) 
			for (int p = 0; p < xs; ++p) {
				T v = dat[q*xs+p];
				if (v < min) min = v;
				if (v > max) max = v;
			}

		if (abs((T)(max-min)) < ((T)1e-11f)) return;

		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				T t = dat[q*xs+p]-min; t /= max-min;
				t *= 0xff;
				if (t > (T)0xff) t = (T)0xff;
				byte R, g, b;
				R = g = b = (byte)t;
				r.dat[q*xs+p] = (R<<16)+(g<<8)+b;
			}
		}

	}
	void GetSignDispImg(SDisImg &r)
	{
		if (!xs || !ys) return;
		if (xs != r.xs || ys != r.ys) {
			r.Clean();
			r.dat = new unsigned long[ys*xs];
			if (!r.dat) return;
			r.xs = xs;
			r.ys = ys;
		}

		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				T t = dat[q*xs+p];
				t *= 0xff;
				if (t > (T)0xff) t = (T)0xff;
				byte R, G, B;
				if (t < 0) { B = -(byte)t; G = -(byte)(2*t/3); R = 0; }
				else { R = (byte)t; G = (byte)t/3; B = 0; }
				r.dat[q*xs+p] = (R<<16)+(G<<8)+B;
			}
		}

	}
	void GetColorDispImg(SDisImg &r)
	{
		if (!xs || !ys) return;
		if (xs != r.xs || ys != r.ys) {
			r.Clean();
			r.dat = new unsigned long[ys*xs];
			if (!r.dat) return;
			r.xs = xs;
			r.ys = ys;
		}

		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				T t = dat[q*xs+p];
				t *= 0xff;
				if (t > (T)0xff) t = (T)0xff;
				byte R(0), G(0), B(0);
				if (t > 5) {
					R = (byte)t;
					B = 0xff-(byte)t;
					G = 0x7f;
				}
				else if (t > 0) {
					B = 0x7f;
				}
				r.dat[q*xs+p] = (R<<16)+(G<<8)+B;
			}
		}

	}
	void GetSaveImg(SSavImg &r)
	{
		if (!xs || !ys) return;
		{
			if (r.dat) delete[] r.dat;
			r.dat = new unsigned char[3*ys*xs];
			if (!r.dat) return;
			r.xs = xs;
			r.ys = ys;
		}

		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				T t = dat[q*xs+p];
				t *= 0xff;
				if (t > (T)0xff) t = (T)0xff;
				byte R, g, b;
				R = g = b = (byte)t;
				r.dat[(q*xs+p)*3] = R;
				r.dat[(q*xs+p)*3+1] = R;
				r.dat[(q*xs+p)*3+2] = R;
			}
		}

	}
	SWorkImg& GetDiv(SWorkImg &u, SWorkImg &v, bool bZeroBound = true) {
		if (u.xs != v.xs || u.ys != v.ys) 
			return *this;
		if (xs != u.xs || ys != u.ys) {
			Clean();
			xs = u.xs;
			ys = u.ys;
			dat = new T[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return *this;
			}
		}
#if 1
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				T div;

				if (!p) 
					div = u[q][p];
				else if (p == xs-1)
					div = -u[q][p-1];
				else 
					div = u[q][p]-u[q][p-1];

				if (!q) 
					div += v[q][p];
				else if (q == ys-1)
					div += -v[q-1][p];
				else 
					div += v[q][p]-v[q-1][p];

				dat[q*xs+p] = div;	
			}
		}
#else
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				T div;

				if (!p) 
					div = u[q][p];
				else if (p == xs-1)
					div = -u[q][p-1];
				else 
					div = u[q][p+1]-u[q][p];

				if (!q) 
					div += v[q][p];
				else if (q == ys-1)
					div += -v[q-1][p];
				else 
					div += v[q][p+1]-v[q][p];

				dat[q*xs+p] = div;	
			}
		}
#endif
		if (bZeroBound) {
			for (int q = 0; q < ys; ++q) {
				dat[(q+1)*xs-1] = dat[q*xs] =  0;
			}
			for (int p = 0; p < xs; ++p) {
				dat[(ys-1)*xs+p] = dat[p] =  0;
			}
		}

		return *this;
	}
	SWorkImg& GetCDiv(SWorkImg &u, SWorkImg &v, bool bZeroBound = true) {
		if (u.xs != v.xs || u.ys != v.ys) 
			return *this;
		if (xs != u.xs || ys != u.ys) {
			Clean();
			xs = u.xs;
			ys = u.ys;
			dat = new T[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return *this;
			}
		}
		for (int q = 1; q < ys-1; ++q) {
			for (int p = 1; p < xs-1; ++p) {
				T div;
				div = u[q][p+1]-u[q][p-1];
				div += v[q+1][p]-v[q-1][p];				
				dat[q*xs+p] = T(0.5)*div;	
			}
		}

		if (bZeroBound) {
			for (int q = 0; q < ys; ++q) {
				dat[(q+1)*xs-1] = dat[q*xs] =  0;
			}
			for (int p = 0; p < xs; ++p) {
				dat[(ys-1)*xs+p] = dat[p] =  0;
			}
		}
		else {
			for (int q = 0; q < ys; ++q) {
				dat[(q+1)*xs-1] = dat[(q+1)*xs-2];
				dat[q*xs] =  dat[q*xs+1];
			}
			for (int p = 0; p < xs; ++p) {
				dat[(ys-1)*xs+p] = dat[(ys-2)*xs+p];
				dat[p] = dat[xs+p];
			}
		}

		return *this;
	}
	void GetImgGrad(SWorkImg &gx, SWorkImg &gy, bool bZeroBound = true)
	{
		if (!xs || !ys) return;
		if (xs != gx.xs || ys != gx.ys) {
			gx.Clean();
			gx.dat = new T[ys*xs];
			if (!gx.dat) return;
			gx.xs = xs;
			gx.ys = ys;
		}
		if (xs != gy.xs || ys != gy.ys) {
			gy.Clean();
			gy.dat = new T[ys*xs];
			if (!gy.dat) return;
			gy.xs = xs;
			gy.ys = ys;
		}

		for (int q = 0; q < ys; ++q)
			for (int p = 0; p < xs-1; ++p)
				gx[q][p] = dat[q*xs+p+1]-dat[q*xs+p];
		if (bZeroBound) {
			for (int q = 0; q < ys; ++q) {
				gx[q][xs-1] = 0;
				gx[q][0] = 0;
			}
		}
		else {
			for (int q = 0; q < ys; ++q)
				gx[q][xs-1] = gx[q][xs-2];
		}

		for (int q = 0; q < ys-1; ++q)
			for (int p = 0; p < xs; ++p)
				gy[q][p] = dat[(q+1)*xs+p]-dat[q*xs+p];
		if (bZeroBound) {
			for (int p = 0; p < xs; ++p) {
				gy[ys-1][p] = 0;
				gy[0][p] = 0;
			}
		}
		else {
			for (int p = 0; p < xs; ++p)
				gy[ys-1][p] = gy[ys-2][p];
		}

	}
	void GetImgCGrad(SWorkImg &gx, SWorkImg &gy, bool bZeroBound = true)
	{
		if (!xs || !ys) return;
		if (xs != gx.xs || ys != gx.ys) {
			gx.Clean();
			gx.dat = new T[ys*xs];
			if (!gx.dat) return;
			gx.xs = xs;
			gx.ys = ys;
		}
		if (xs != gy.xs || ys != gy.ys) {
			gy.Clean();
			gy.dat = new T[ys*xs];
			if (!gy.dat) return;
			gy.xs = xs;
			gy.ys = ys;
		}

		for (int q = 1; q < ys-1; ++q)
			for (int p = 1; p < xs-1; ++p)
				gx[q][p] = T(0.5)*(dat[q*xs+p+1]-dat[q*xs+p-1]);
		for (int q = 1; q < ys-1; ++q)
			for (int p = 1; p < xs-1; ++p)
				gy[q][p] = T(0.5)*(dat[(q+1)*xs+p]-dat[(q-1)*xs+p]);

		if (bZeroBound) {
			for (int q = 0; q < ys; ++q) {
				gx[q][xs-1] = gx[q][0] =  0;
				gy[q][xs-1] = gy[q][0] =  0;
			}
			for (int p = 0; p < xs; ++p) {
				gx[ys-1][p] = gx[0][p] =  0;
				gy[ys-1][p] = gy[0][p] =  0;
			}
		}
		else {
			for (int q = 0; q < ys; ++q) {
				gx[q][xs-1] = gx[q][xs-2];
				gx[q][0] =  gx[q][1];
				gy[q][xs-1] = gy[q][xs-2];
				gy[q][0] =  gy[q][1];
			}
			for (int p = 0; p < xs; ++p) {
				gx[ys-1][p] = gx[ys-2][p];
				gx[0][p] =  gx[1][p];
				gy[ys-1][p] = gy[ys-2][p];
				gy[0][p] =  gy[1][p];
			}
		}
	}
	inline T GetInterpolated(T x, T y)
	{
		int ix = (int)x;
		int iy = (int)y;
		T mx = x-ix;
		T my = y-iy;
		if (ix < 0) ix = 0;
		if (ix > xs-2) ix = xs-2;
		if (iy < 0) iy = 0;
		if (iy > ys-2) iy = ys-2;
		T &r00 = dat[iy*xs+ix];
		T &r10 = dat[iy*xs+ix+1];
		T &r01 = dat[(iy+1)*xs+ix];
		T &r11 = dat[(iy+1)*xs+ix+1];
		T wx = ((T)1)-mx;
		T wy = ((T)1)-my;
		return r00*wx*wy+r10*mx*wy+r01*my*wx+r11*mx*my;
	}
	void Get1stDirectionalC(SWorkImg &gd, T dx, T dy)
	{
		if (!xs || !ys) return;
		if (xs != gd.xs || ys != gd.ys) {
			gd.Clean();
			gd.dat = new T[ys*xs];
			if (!gd.dat) return;
			gd.xs = xs;
			gd.ys = ys;
		}
		for (int q = 0; q < ys; ++q)
			for (int p = 0; p < xs; ++p) {
				T vp = GetInterpolated(p+dx,q+dy);
				T vm = GetInterpolated(p-dx,q-dy);
				gd[q][p] = T(0.5)*(vp-vm);
			}

		for (int q = 0; q < ys; ++q) {
			gd[q][xs-1] = gd[q][0] =  0;
		}
		for (int p = 0; p < xs; ++p) {
			gd[ys-1][p] = gd[0][p] =  0;
		}
	}
	void Get2ndDirectionalC(SWorkImg &gd, T dx, T dy)
	{
		if (!xs || !ys) return;
		if (xs != gd.xs || ys != gd.ys) {
			gd.Clean();
			gd.dat = new T[ys*xs];
			if (!gd.dat) return;
			gd.xs = xs;
			gd.ys = ys;
		}
		for (int q = 0; q < ys; ++q)
			for (int p = 0; p < xs; ++p) {
				T vp = GetInterpolated(p+dx,q+dy);
				T vm = GetInterpolated(p-dx,q-dy);
				gd[q][p] = vp+vm-2*dat[q*xs+p];
			}

		for (int q = 0; q < ys; ++q) {
			gd[q][xs-1] = gd[q][0] =  0;
		}
		for (int p = 0; p < xs; ++p) {
			gd[ys-1][p] = gd[0][p] =  0;
		}
	}
	void GetImgHesse(SWorkImg &gxx, SWorkImg &gxy, SWorkImg &gyy)
	{
		if (!xs || !ys) return;
		if (xs != gxx.xs || ys != gxx.ys) {
			gxx.Clean();
			gxx.dat = new T[ys*xs];
			if (!gxx.dat) return;
			gxx.xs = xs;
			gxx.ys = ys;
		}
		if (xs != gxy.xs || ys != gxy.ys) {
			gxy.Clean();
			gxy.dat = new T[ys*xs];
			if (!gxy.dat) return;
			gxy.xs = xs;
			gxy.ys = ys;
		}
		if (xs != gyy.xs || ys != gyy.ys) {
			gyy.Clean();
			gyy.dat = new T[ys*xs];
			if (!gyy.dat) return;
			gyy.xs = xs;
			gyy.ys = ys;
		}
		for (int q = 1; q < ys-1; ++q)
			for (int p = 1; p < xs-1; ++p)
				gxx[q][p] = (dat[q*xs+p+1]+dat[q*xs+p-1]-2*dat[q*xs+p]);
		for (int q = 1; q < ys-1; ++q)
			for (int p = 1; p < xs-1; ++p)
				gxy[q][p] = T(0.25)*(dat[(q+1)*xs+p+1]+dat[(q-1)*xs+p-1]-dat[(q+1)*xs+p-1]-dat[(q-1)*xs+p+1]);
		for (int q = 1; q < ys-1; ++q)
			for (int p = 1; p < xs-1; ++p)
				gyy[q][p] = (dat[(q+1)*xs+p]+dat[(q-1)*xs+p]-2*dat[q*xs+p]);

		for (int q = 0; q < ys; ++q) {
			gxx[q][xs-1] = gxx[q][0] =  0;
			gxy[q][xs-1] = gxy[q][0] =  0;
			gyy[q][xs-1] = gyy[q][0] =  0;
		}
		for (int p = 0; p < xs; ++p) {
			gxx[ys-1][p] = gxx[0][p] =  0;
			gxy[ys-1][p] = gxy[0][p] =  0;
			gyy[ys-1][p] = gyy[0][p] =  0;
		}
	}
	SWorkImg& GetMagnitude(SWorkImg &u, SWorkImg &v, T beta = T(0)) {
		if (u.xs != v.xs || u.ys != v.ys) 
			return *this;
		if (xs != u.xs || ys != u.ys) {
			Clean();
			xs = u.xs;
			ys = u.ys;
			dat = new T[ys*xs];
			if (!dat) {
				xs = ys = 0;
				return *this;
			}
		}
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				T mag = u[q][p]*u[q][p]+v[q][p]*v[q][p];

				dat[q*xs+p] = sqrt(mag+beta);	
			}
		}

		return *this;
	}
	void GetLaplace(SWorkImg &s) {
		if (xs != s.xs || ys != s.ys) 
			return;
		for (int yy = 1; yy < ys-1; ++yy) {
			for (int xx = 1; xx < xs-1; ++xx) {
				dat[yy*xs+xx] = -4*s[yy][xx]+s[yy][xx+1]+s[yy][xx-1]+s[yy+1][xx]+s[yy-1][xx];
			}
		}
		int y = ys-1;
		for (int xx = 1; xx < xs-1; ++xx) {
			dat[xx] = -4*s[0][xx]+s[0][xx+1]+s[0][xx-1]+s[1][xx]+s[1][xx];
			dat[y*xs+xx] = -4*s[y][xx]+s[y][xx+1]+s[y][xx-1]+s[y-1][xx]+s[y-1][xx];
		}
		int x = xs-1;
		for (int yy = 1; yy < ys-1; ++yy) {
			dat[yy*xs+0] = -4*s[yy][0]+s[yy+1][0]+s[yy-1][0]+s[yy][1]+s[yy][1];
			dat[yy*xs+x] = -4*s[yy][x]+s[yy+1][x]+s[yy-1][x]+s[yy][x-1]+s[yy][x-1];
		}
		dat[0*xs+0] = 0.5f*(dat[0*xs+1]+dat[1*xs+0]);
		dat[0*xs+x] = 0.5f*(dat[0*xs+x-1]+dat[1*xs+x]);
		dat[y*xs+0] = 0.5f*(dat[(y-1)*xs+0]+dat[y*xs+1]);
		dat[y*xs+x] = 0.5f*(dat[(y-1)*xs+x]+dat[y*xs+x-1]);

	}


	SWorkImg& operator= (T r) {
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] = r;	
			}
		}
		return *this;
	}
	SWorkImg& operator*= (T r) {
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] *= r;	
			}
		}
		return *this;
	}
	SWorkImg& operator+= (T r) {
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] += r;	
			}
		}
		return *this;
	}
	SWorkImg& operator-= (T r) {
		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] -= r;	
			}
		}
		return *this;
	}
	SWorkImg& operator+= (SWorkImg &tc) {
		if (xs != tc.xs || ys != tc.ys) {
			return *this;
		}

		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] += tc[q][p];
			}
		}
		return *this;
	}
	SWorkImg& operator-= (SWorkImg &tc) {
		if (xs != tc.xs || ys != tc.ys) {
			return *this;
		}

		for (int q = 0; q < ys; ++q) {
			for (int p = 0; p < xs; ++p) {
				dat[q*xs+p] -= tc[q][p];
			}
		}
		return *this;
	}

	T *dat;
	T maxval;
	T minval;
	T avgval;
	int xs;
	int ys;
	//bool bnormed;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


template <class T> struct SVoxImg {
	SVoxImg() {
		dat = 0;
		xs = ys = zs = 0;
	}
	SVoxImg(int x, int y, int z) {
		dat = 0;
		xs = ys = zs = 0;
		Set(x,y,z);
	}
	void Set(int x, int y, int z) {
		if (dat && (x != xs || y != ys || z != zs))
			Clean();
		xs = x;
		ys = y;
		zs = z;

		if (!dat && xs && ys && zs) {
			dat = new T[zs];
			if (!dat) {
				xs = ys = zs = 0;
				return;
			}
		}
		for (int q = 0; q < zs; ++q) {
			dat[q].Set(xs,ys,-1);
		}
	}
	void Set0(int x, int y, int z) {
		if (dat && (x != xs || y != ys || z != zs))
			Clean();
		xs = x;
		ys = y;
		zs = z;

		if (!dat && xs && ys && zs) {
			dat = new T[zs];
			if (!dat) {
				xs = ys = zs = 0;
				return;
			}
		}
		for (int q = 0; q < zs; ++q) {
			dat[q].Set(xs,ys,0.0);
		}
	}
	SVoxImg& operator= (SVoxImg &tc) {
		if (&tc == this) return *this;
		if (xs != tc.xs || ys != tc.ys || zs != tc.zs) {
			Clean();
			xs = tc.xs;
			ys = tc.ys;
			zs = tc.zs;
			dat = new T[zs];
			if (!dat) {
				xs = ys = zs = 0;
				return *this;
			}
		}

		for (int q = 0; q < zs; ++q) {
			dat[q] = tc.dat[q];
		}
		return *this;
	}

	void GetZLoopLaplace(SVoxImg &s) {
		if (xs != s.xs || ys != s.ys || zs != s.zs) 
			return;
#pragma omp parallel for
		for (int zz = 0; zz < zs; ++zz) {
			int zp(zz+1); if (zp == zs) zp = 0;
			int zm(zz-1); if (zm == -1) zm = zs-1;
			for (int yy = 0; yy < ys; ++yy) {
				int yp(yy+1); if (yp == ys) yp = ys-2;
				int ym(yy-1); if (ym == -1) ym = 1;
				for (int xx = 0; xx < xs; ++xx) {
					int xp(xx+1); if (xp == xs) xp = xs-2;
					int xm(xx-1); if (xm == -1) xm = 1;
					dat[zz][yy][xx] = -6*s[zz][yy][xx]
					+s[zz][yy][xp]+s[zz][yy][xm]+s[zz][yp][xx]+s[zz][ym][xx]+s[zp][yy][xx]+s[zm][yy][xx];
				}
			}
		}

	}

	void GetZCLoopGrad(SVoxImg &gx, SVoxImg &gy, SVoxImg &gz)
	{
		if (xs != gx.xs || ys != gx.ys || zs != gx.zs) 
			return;
		if (xs != gy.xs || ys != gy.ys || zs != gy.zs) 
			return;
		if (xs != gz.xs || ys != gz.ys || zs != gz.zs) 
			return;

#pragma omp parallel for
		for (int r = 0; r < zs; ++r) {
			int rm(r-1); if (rm == -1) rm = zs-1;
			int rp(r+1); if (rp == zs) rp = 0;

			for (int q = 0; q < ys; ++q) {
				int qp(q+1); if (qp == ys) qp = ys-2;
				int qm(q-1); if (qm == -1) qm = 1;
				for (int p = 0; p < xs; ++p) {
					int pp(p+1); if (pp == xs) pp = xs-2;
					int pm(p-1); if (pm == -1) pm = 1;
					gx[r][q][p] = 0.5*(dat[r][q][pp]-dat[r][q][pm]);
					gy[r][q][p] = 0.5*(dat[r][qp][p]-dat[r][qm][p]);
					gz[r][q][p] = 0.5*(dat[rp][q][p]-dat[rm][q][p]);
				}
			}
		}

		
	}
	void GetZFLoopGrad(SVoxImg &gx, SVoxImg &gy, SVoxImg &gz)
	{
		if (xs != gx.xs || ys != gx.ys || zs != gx.zs) 
			return;
		if (xs != gy.xs || ys != gy.ys || zs != gy.zs) 
			return;
		if (xs != gz.xs || ys != gz.ys || zs != gz.zs) 
			return;

		for (int r = 0; r < zs; ++r) {
			//int rm(r-1); if (rm == -1) rm = zs-1;
			int rp(r+1); if (rp == zs) rp = 0;

			for (int q = 0; q < ys; ++q) {
				int qp(q+1); if (qp == ys) qp = ys-2;
				//int qm(q-1); if (qm == -1) qm = 1;
				for (int p = 0; p < xs; ++p) {
					int pp(p+1); if (pp == xs) pp = xs-2;
					//int pm(p-1); if (pm == -1) pm = 1;
					gx[r][q][p] = (dat[r][q][pp]-dat[r][q][p]);
					gy[r][q][p] = (dat[r][qp][p]-dat[r][q][p]);
					gz[r][q][p] = (dat[rp][q][p]-dat[r][q][p]);
				}
			}
		}

		
	}

	void GetZBLoopGrad(SVoxImg &gx, SVoxImg &gy, SVoxImg &gz)
	{
		if (xs != gx.xs || ys != gx.ys || zs != gx.zs) 
			return;
		if (xs != gy.xs || ys != gy.ys || zs != gy.zs) 
			return;
		if (xs != gz.xs || ys != gz.ys || zs != gz.zs) 
			return;

		for (int r = 0; r < zs; ++r) {
			int rm(r-1); if (rm == -1) rm = zs-1;

			for (int q = 0; q < ys; ++q) {
				int qm(q-1); if (qm == -1) qm = 1;
				for (int p = 0; p < xs; ++p) {
					int pm(p-1); if (pm == -1) pm = 1;
					gx[r][q][p] = (dat[r][q][p]-dat[r][q][pm]);
					gy[r][q][p] = (dat[r][q][p]-dat[r][qm][p]);
					gz[r][q][p] = (dat[r][q][p]-dat[rm][q][p]);
				}
			}
		}

		
	}

	//////////////////////////////
	//// standard derivatives ////
	//////////////////////////////
	void GetLaplace(SVoxImg &s) {
		if (xs != s.xs || ys != s.ys || zs != s.zs) 
			return;
#pragma omp parallel for
		for (int zz = 0; zz < zs; ++zz) {
			int zp(zz+1); if (zp == zs) zp = zs-2;
			int zm(zz-1); if (zm == -1) zm = 1;
			for (int yy = 0; yy < ys; ++yy) {
				int yp(yy+1); if (yp == ys) yp = ys-2;
				int ym(yy-1); if (ym == -1) ym = 1;
				for (int xx = 0; xx < xs; ++xx) {
					int xp(xx+1); if (xp == xs) xp = xs-2;
					int xm(xx-1); if (xm == -1) xm = 1;
					dat[zz][yy][xx] = -6*s[zz][yy][xx]
					+s[zz][yy][xp]+s[zz][yy][xm]+s[zz][yp][xx]+s[zz][ym][xx]+s[zp][yy][xx]+s[zm][yy][xx];
				}
			}
		}

	}
	void GetGrad(SVoxImg &gx, SVoxImg &gy, SVoxImg &gz)
	{
		if (xs != gx.xs || ys != gx.ys || zs != gx.zs) 
			return;
		if (xs != gy.xs || ys != gy.ys || zs != gy.zs) 
			return;
		if (xs != gz.xs || ys != gz.ys || zs != gz.zs) 
			return;

#pragma omp parallel for
		for (int z = 0; z < zs; ++z) {
			int zm(z-1); if (zm == -1) zm = 1;
			int zp(z+1); if (zp == zs) zp = zs-2;
			for (int y = 0; y < ys; ++y) {
				int yp(y+1); if (yp == ys) yp = ys-2;
				int ym(y-1); if (ym == -1) ym = 1;
				for (int x = 0; x < xs; ++x) {
					int xp(x+1); if (xp == xs) xp = xs-2;
					int xm(x-1); if (xm == -1) xm = 1;
					gx[z][y][x] = 0.5*(dat[z][y][xp]-dat[z][y][xm]);
					gy[z][y][x] = 0.5*(dat[z][yp][x]-dat[z][ym][x]);
					gz[z][y][x] = 0.5*(dat[zp][y][x]-dat[zm][y][x]);
				}
			}
		}
		
	}
	void GetHess(SVoxImg &hxx, SVoxImg &hyy, SVoxImg &hzz, SVoxImg &hxy, SVoxImg &hxz, SVoxImg &hyz)
	{
		if (xs != hxx.xs || ys != hxx.ys || zs != hxx.zs) 
			return;
		if (xs != hyy.xs || ys != hyy.ys || zs != hyy.zs) 
			return;
		if (xs != hzz.xs || ys != hzz.ys || zs != hzz.zs) 
			return;
		if (xs != hxy.xs || ys != hxy.ys || zs != hxy.zs) 
			return;
		if (xs != hxz.xs || ys != hxz.ys || zs != hxz.zs) 
			return;
		if (xs != hyz.xs || ys != hyz.ys || zs != hyz.zs) 
			return;

#pragma omp parallel for
		for (int z = 0; z < zs; ++z) {
			int zm(z-1); if (zm == -1) zm = 1;
			int zp(z+1); if (zp == zs) zp = zs-2;
			for (int y = 0; y < ys; ++y) {
				int yp(y+1); if (yp == ys) yp = ys-2;
				int ym(y-1); if (ym == -1) ym = 1;
				for (int x = 0; x < xs; ++x) {
					int xp(x+1); if (xp == xs) xp = xs-2;
					int xm(x-1); if (xm == -1) xm = 1;


					hxx[z][y][x] = dat[z][y][xp]+dat[z][y][xm]-2*dat[z][y][x];
					hyy[z][y][x] = dat[z][yp][x]+dat[z][ym][x]-2*dat[z][y][x];
					hzz[z][y][x] = dat[zp][y][x]+dat[zm][y][x]-2*dat[z][y][x];
					
					hxy[z][y][x] = 0.25*(dat[z][yp][xp]+dat[z][ym][xm]-dat[z][yp][xm]-dat[z][ym][xp]);
					hxz[z][y][x] = 0.25*(dat[zp][y][xp]+dat[zm][y][xm]-dat[zp][y][xm]-dat[zm][y][xp]);
					hyz[z][y][x] = 0.25*(dat[zp][yp][x]+dat[zm][ym][x]-dat[zp][ym][x]-dat[zm][yp][x]);
				}
			}
		}
		
	}


	void GetGradLen(SVoxImg &glen)
	{
		if (xs != glen.xs || ys != glen.ys || zs != glen.zs) 
			return;

#pragma omp parallel for
		for (int r = 0; r < zs; ++r) {
			int rm(r-1); if (rm == -1) rm = 1;
			int rp(r+1); if (rp == zs) rp = zs-2;
			for (int q = 0; q < ys; ++q) {
				int qp(q+1); if (qp == ys) qp = ys-2;
				int qm(q-1); if (qm == -1) qm = 1;
				for (int p = 0; p < xs; ++p) {
					int pp(p+1); if (pp == xs) pp = xs-2;
					int pm(p-1); if (pm == -1) pm = 1;

					realnum gx = 0.5*(dat[r][q][pp]-dat[r][q][pm]);
					realnum gy = 0.5*(dat[r][qp][p]-dat[r][qm][p]);
					realnum gz = 0.5*(dat[rp][q][p]-dat[rm][q][p]);
					glen[r][q][p] = sqrt(gx*gx+gy*gy+gz*gz);
				}
			}
		}

		
	}


	SVoxImg& operator*= (realnum r) {
		for (int q = 0; q < zs; ++q) {
			dat[q] *= r;	
		}
		return *this;
	}


	void Clean() {
		if (dat) {
			delete[] dat;
			dat = 0;
			xs = ys = zs = 0;
		}
	}
	~SVoxImg() {
		Clean();
	}




	T& operator[] (int z) {
		if (z >= zs) z = zs-1;
		else if (z < 0) z = 0;
		return dat[z];
	}
	const T& operator[](int z) const {
		int iz = z;
		if (iz >= zs) iz = zs-1;
		else if (iz < 0) iz = 0;
		return dat[iz];
	}
	/**/



	T *dat;
	int xs;
	int ys;
	int zs;
};