// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef LLS_h
#define LLS_h

#include "Array.h"
#include "Stack.h"
#include "Pool.h"

typedef Array<double> FVector;

// class FVector : public Array<double> {
//   public:
// 	FVector(int pn=0) : Array<double>(pn) { }
// 	~FVector() { }
// };

extern double mag2(const FVector& v);
extern double dot(const FVector& va, const FVector& vb);

class FMatrix : public SArray<FVector> {
  public:
	FMatrix(int pm, int pn) : SArray<FVector>(pm) {
		for (int i=0;i<pm;i++) a[i].init(pn);
	}
	~FMatrix() { }
};

class LLS {
  public:
	// virtual constructor
	static LLS* create(int m, int n, int nd, double nonzerofrac);
	virtual ~LLS() { }
	// All entries will be zero unless entered as below.
	virtual void enter(int r, int c, double val)=0; // r<m, c<n
	void enterRh(int r, const double* dv);	       // r<m, dv[nd]
	void enterEstimate(int vn, const double* dv);   // vn<n, dv[nd]
	virtual void solve(double* rssb=0, double* rssa=0)=0;
	void getValue(int vn, double* dv); // vn<n, dv[nd]
  protected:
	LLS(int pm, int pn, int pnd);
	int m, n, nd;
	FMatrix rh;		// [nd][m]
	FMatrix sol;		// [nd][n];
  private:
	DISABLECOPY(LLS);
};

class SparseLLS : public LLS {
  public:
	SparseLLS(int pm, int pn, int pnd);
	~SparseLLS();
	void enter(int r, int c, double val);
	void solve(double* rssb=0, double* rssa=0);
	void tolerance(double ptolerance);
	void maxiter(int pmaxiter);
  private:
	struct Ival {
		POOLALLOCATION(Ival);
		Ival(int pi, double pv) : i(pi), v(pv) { }
		int i;
		double v;
	};
	SArray<Stack<Ival*> > rows;
	SArray<Stack<Ival*> > cols;
	double itolerance;
	int imaxiter;
	int nentries;
	void multmv(const FVector& vi, FVector& vo) const;
	void multmtv(const FVector& vi, FVector& vo) const;
	void docg(FVector& x, FVector& h,
		  double* prssb, double* prssa) const;
	DISABLECOPY(SparseLLS);
};

class FullLLS : public LLS {
  public:
	FullLLS(int pm, int pn, int pnd);
	~FullLLS();
	void enter(int r, int c, double val);
	void solve(double* rssb=0, double* rssa=0);
  protected:
	FMatrix u;		// [m][n]
	virtual void solveaux()=0; // abstract class
  private:
	double getrss();
};

// Slow, for validation only
class LudLLS : public FullLLS {
  public:
	LudLLS(int pm, int pn, int pnd);
	~LudLLS();
  protected:
	void solveaux();
};

class GivensLLS : public FullLLS {
  public:
	GivensLLS(int pm, int pn, int pnd);
	~GivensLLS();
  protected:
	void solveaux();
};

class SvdLLS : public FullLLS {
  public:
	SvdLLS(int pm, int pn, int pnd);
	~SvdLLS();
  protected:
	void solveaux();
};

class QrdLLS : public FullLLS {
  public:
	QrdLLS(int pm, int pn, int pnd);
	~QrdLLS();
  protected:
	void solveaux();
};

#endif
