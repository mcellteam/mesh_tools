// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "LLS.h"		// ludcmp(), lubksp(), svdcmp(), svbksb()
#include "Stack.h"
#include "Recipes.h"
#include "Linpack.h"		// sqrdc_(), sqrsl_()
#include "Timer.h"

static int debug=GetenvValue("DEBUG");

//*** linear algebra

double mag2(const FVector& v)
{
	register double sum=0;
	ForIndex(i,v.num()) { double a=v[i]; sum+=a*a; } EndFor;
	return sum;
}

double dot(const FVector& va, const FVector& vb)
{
	assertx(va.num()==vb.num());
	register double sum=0;
	// ForIndex(i,va.num()) { sum+=va[i]*vb[i]; } EndFor;
	// yuck, but that is what is needed.
	int s=4, n=va.num(), n1=n-s+1, i=0;
	const double* fa=va; const double* fb=vb;
	for (;i<n1;i+=s)
		sum+=(double(fa[i])*fb[i]+double(fa[i+1])*fb[i+1]+
		      double(fa[i+2])*fb[i+2]+double(fa[i+3])*fb[i+3]);
	for (;i<n;i++)
		sum+=double(fa[i])*fb[i];
	return sum;
}

//*** virtual constructor

LLS* LLS::create(int m, int n, int nd, double nonzerofrac)
{
	if (GetenvValue("LUDLLS"))	return new LudLLS(m,n,3);
	if (GetenvValue("GIVENSLLS"))	return new GivensLLS(m,n,3);
	if (GetenvValue("SVDLLS"))	return new SvdLLS(m,n,3);
	if (GetenvValue("QRDLLS"))	return new QrdLLS(m,n,3);
	if (GetenvValue("SPARSELLS"))	return new SparseLLS(m,n,3);
	double size=double(m)*n;
	if (size<1000*40) {	// small system
		return new QrdLLS(m,n,3);
	} else {		// large system
		if (nonzerofrac<.3) return new SparseLLS(m,n,3);
		return new QrdLLS(m,n,3);
	}
}

//*** LLS

LLS::LLS(int pm, int pn, int pnd)
: m(pm), n(pn), nd(pnd), rh(nd, m), sol(nd, n)
{
	ForIndex(di,nd) {
		ForIndex(i,m) { rh[di][i]=0; } EndFor;
		ForIndex(j,n) { sol[di][j]=0; } EndFor;
	} EndFor;
}

void LLS::enterRh(int r, const double* dv)
{
	ForIndex(di,nd) { rh[di][r]=dv[di]; } EndFor;
}

void LLS::enterEstimate(int vn, const double* dv)
{
	ForIndex(di,nd) { sol[di][vn]=dv[di]; } EndFor;
}

void LLS::getValue(int vn, double* dv)
{
	ForIndex(di,nd) { dv[di]=sol[di][vn]; } EndFor;
}

//*** SparseLLS

ALLOCATEPOOL(SparseLLS::Ival);

SparseLLS::SparseLLS(int pm, int pn, int pnd)
: LLS(pm,pn,pnd), rows(m), cols(n),
  itolerance(1e-10), imaxiter(0), nentries(0)
{ }

SparseLLS::~SparseLLS()
{
	ForIndex(i,m) {
		while (!rows[i].empty()) delete rows[i].pop();
	} EndFor;
	ForIndex(j,n) {
		while (!cols[j].empty()) delete cols[j].pop();
	} EndFor;
}

void SparseLLS::enter(int r, int c, double val)
{
	rows[r].push(new Ival(c,val));
	cols[c].push(new Ival(r,val));
	nentries++;
}

void SparseLLS::multmv(const FVector& vi, FVector& vo) const
{
	// vo[m]=u[m][n]*vi[n];
	const double* via=vi; double* voa=vo; // unchecked
	ForIndex(i,m) {
		register double sum=0;
		// ForStack(rows[i],Ival*,n) sum+=n->v*via[n->i]; EndFor;
		// GNU v2.4.5 doesnt expand inline Stack templates, so:
		const BStack& stack=rows[i];
		For(BStackIter si(stack);si;si.next()) {
			register Ival* n = Conv<Ival *>::d(si());
			sum+=n->v*via[n->i];
		} EndFor;
		voa[i]=sum;
	} EndFor;
}

void SparseLLS::multmtv(const FVector& vi, FVector& vo) const
{
	// vo[n]=uT[n][m]*vi[m];
	const double* via=vi; double* voa=vo; // unchecked
	ForIndex(j,n) {
		register double sum=0;
		// ForStack(cols[j],Ival*,n) sum+=n->v*via[n->i]; EndFor;
		// GNU v2.4.5 doesnt expand inline Stack templates, so:
		const BStack& stack=cols[j];
		For(BStackIter si(stack);si;si.next()) {
			register Ival* n = Conv<Ival*>::d(si());
			sum+=n->v*via[n->i];
		} EndFor;
		voa[j]=sum;
	} EndFor;
}

void SparseLLS::solve(double* prssb, double* prssa)
{
	ATIMER(_____sparseLLS);
	if (debug) SHOWF("SparseLLS: solving %dx%d system, nonzerofrac=%f\n",
			 m,n,double(nentries)/m/n);
	if (prssb) *prssb=0;
	if (prssa) *prssa=0;
	FVector x(n), rhv(m);
	ForIndex(di,nd) {
		ForIndex(i,m) { rhv[i]=rh[di][i]; } EndFor;
		ForIndex(j,n) { x[j]=sol[di][j]; } EndFor;
		docg(x,rhv,prssb,prssa);
		ForIndex(j,n) { sol[di][j]=x[j]; } EndFor;
	} EndFor;
}

void SparseLLS::docg(FVector& x, FVector& h,
		     double* prssb, double* prssa) const
{
	// x(n), h(m)
	int i,k;
	FVector rc(m), gc(n), gp(n), dc(n), tc(m);
	multmv(x,rc);
	for (i=0;i<m;i++) rc[i]-=h[i];
	double rssb=mag2(rc);
	multmtv(rc,gc);
	for (i=0;i<n;i++) gc[i]=-gc[i];
	for (k=0;k<n;k++) {
		double gm2=mag2(gc);
		if (gm2<itolerance || imaxiter && k==imaxiter) break;
		if (k>0) {
			double bi=gm2/mag2(gp);
			for (i=0;i<n;i++) dc[i]=gc[i]+bi*dc[i];
		} else {
			for (i=0;i<n;i++) dc[i]=gc[i];
		}
		multmv(dc,tc);
		double ai=gm2/mag2(tc);
		for (i=0;i<n;i++) x[i]+=ai*dc[i];
		for (i=0;i<m;i++) rc[i]+=ai*tc[i];
		for (i=0;i<n;i++) gp[i]=gc[i];
		multmtv(rc,gc);
		for (i=0;i<n;i++) gc[i]=-gc[i];
	}
	double rssa=mag2(rc);
	// Print final gradient norm squared and final residual norm squared.
	if (debug) SHOWF("CG: %d iter (gm2=%g, rssb=%g, rssa=%g)\n",
			 k,mag2(gc),rssb,rssa);
	if (prssb) *prssb+=rssb;
	if (prssa) *prssa+=rssa;
}

void SparseLLS::tolerance(double ptolerance)
{
	itolerance=ptolerance;
}

void SparseLLS::maxiter(int pmaxiter)
{
	imaxiter=pmaxiter;
}

//*** FullLLS

FullLLS::FullLLS(int pm, int pn, int pnd)
: LLS(pm,pn,pnd), u(m,n)
{
	ForIndex(i,m) { ForIndex(j,n) { u[i][j]=0; } EndFor; } EndFor;
}

FullLLS::~FullLLS() { }

void FullLLS::enter(int r, int c, double val)
{
	u[r][c]=val;
}

void FullLLS::solve(double* prssb, double* prssa)
{
	assertx(m>=n);
	if (prssb) *prssb=getrss();
	solveaux();
	if (prssa) *prssa=getrss();
}

double FullLLS::getrss()
{
	double rss=0;
	ForIndex(di,nd) {
		ForIndex(i,m) {
			double a=dot(u[i],sol[di])-double(rh[di][i]);
			rss+=a*a;
		} EndFor;
	} EndFor;
	return rss;
}

//*** LudLLS

LudLLS::LudLLS(int pm, int pn, int pnd)
: FullLLS(pm,pn,pnd) { }

LudLLS::~LudLLS() { }

void LudLLS::solveaux()
{
        int j;
	ATIMER(_____ludLLS);
	double** a=matrix(1,n,1,n);
	int* indx=ivector(1,n);
	double* b=vector(1,n);
	for (int i=1;i<=n;i++)
		for (int j=1;j<=n;j++) {
			double s=0;
			for (int k=1;k<=m;k++)
				s+=u[k-1][i-1]*u[k-1][j-1];
			a[i][j]=s;
		}
	/* a[1..n][1..n] -> a[1..n][1..n], indx[1..n], d={-1.,1.} */
	double d;
	ludcmp(a,n,indx,&d);
	ForIndex(di,nd) {
//		for (int j=1;j<=n;j++) {
		for (j=1;j<=n;j++) {
			double s=0;
			for (int i=1;i<=m;i++)
				s+=u[i-1][j-1]*rh[di][i-1];
			b[j]=s;
		}
		/* a,n,indx, b[1..n] -> b[1..n] */
		lubksb(a,n,indx,b);
		for (j=1;j<=n;j++)
			sol[di][j-1]=b[j];
	} EndFor;
	free_matrix(a,1,n,1,n);
	free_ivector(indx,1,n);
	free_vector(b,1,n);
}

//*** GivensLLS

GivensLLS::GivensLLS(int pm, int pn, int pnd)
: FullLLS(pm,pn,pnd) { }

GivensLLS::~GivensLLS() { }

void GivensLLS::solveaux()
{
	ATIMER(_____givensLLS);
	int nposs=0,ngivens=0;
	// Destroys both tu and trh
	FMatrix tu(m,n);
	FMatrix trh(m,nd);
	ForIndex(i,m) { ForIndex(j,n) { tu[i][j]=u[i][j]; } EndFor; } EndFor;
	ForIndex(di,nd) {
		ForIndex(i,m) { trh[i][di]=rh[di][i]; } EndFor;
	} EndFor;
	int i,j,k;
	for (i=0;i<n;i++) {
		for (k=i+1;k<m;k++) {
			nposs++;
			if (!tu[k][i]) continue;
			ngivens++;
			double xi=tu[i][i];
			double xk=tu[k][i];
			register double c,s;
			if (fabs(xk)>fabs(xi)) {
				double t=xi/xk; s=1/mysqrt(1+t*t); c=s*t;
			} else {
				double t=xk/xi; c=1/mysqrt(1+t*t); s=c*t;
			}
			for (j=i;j<n;j++) {
				double xij=tu[i][j], xkj=tu[k][j];
				tu[i][j]=c*xij+s*xkj;
				tu[k][j]= -s*xij+c*xkj;
			}
			ForIndex(di,nd) {
				double xij=trh[i][di], xkj=trh[k][di];
				trh[i][di]=c*xij+s*xkj;
				trh[k][di]= -s*xij+c*xkj;
			} EndFor;
		}
	}
	if (debug) SHOWF("Givens: %d/%d rotations done\n",ngivens,nposs);
	// Backsubstitutions
	ForIndex(di,nd) {
		for (i=n-1;i>=0;i--) {
			double sum=trh[i][di];
			for (j=i+1;j<n;j++)
				sum-=tu[i][j]*trh[j][di];
			trh[i][di]=sum/tu[i][i];
			sol[di][i]=trh[i][di];
		}
	} EndFor;
}

//*** SvdLLS

SvdLLS::SvdLLS(int pm, int pn, int pnd)
: FullLLS(pm,pn,pnd) { }

SvdLLS::~SvdLLS() { }

void SvdLLS::solveaux()
{
	ATIMER(_____svdLLS);
	int i,j;
	double** a=matrix(1,m,1,n);
	double* w=vector(1,n);
	double** v=matrix(1,n,1,n);
	for (i=1;i<=m;i++) for (j=1;j<=n;j++) a[i][j]=u[i-1][j-1];
	// Given a[1..m][1..n], compute singular value decomp.: a=u*w*v^t
	// store u in a; w[1..n]; v[1..n][1..n]
	svdcmp(a,m,n,w,v);
	// Backsubstitutions
	double* b=vector(1,m);
	double* x=vector(1,n);
	ForIndex(di,nd) {
		for (i=1;i<=m;i++) b[i]=rh[di][i-1];
		svbksb(a,w,v,m,n,b,x);
		for (j=1;j<=n;j++) sol[di][j-1]=x[j];
	} EndFor;
	free_vector(x,1,n);
	free_vector(b,1,m);
	free_matrix(v,1,n,1,n);
	free_vector(w,1,n);
	free_matrix(a,1,m,1,n);
}

//*** QrdLLS

QrdLLS::QrdLLS(int pm, int pn, int pnd)
: FullLLS(pm,pn,pnd) { }

QrdLLS::~QrdLLS() { }

void QrdLLS::solveaux()
{
	ATIMER(_____qrdLLS);
	Array<double> a(m*n);
	double* ap=a;
	ForIndex(j,n) { ForIndex(i,m) { *ap++=u[i][j]; } EndFor; } EndFor;
	Array<double> lrh(m), qraux(n), work(n), qty(m), b(n);
	Array<int> jpvt(n); ForIndex(i,n) { jpvt[i]=0; } EndFor;
	int pivot=1, modeb=100, info;
	sqrdc_(a,&m,&m,&n,qraux,jpvt,work,&pivot);
	ForIndex(j,n) { jpvt[j]--; } EndFor; // fortran indices start at 1
	ForIndex(di,nd) {
		ForIndex(i,m) { lrh[i]=rh[di][i]; } EndFor;
		sqrsl_(a,&m,&m,&n,qraux,lrh,0,qty,b,0,0,&modeb,&info);
		if (assertw1(!info)) continue; // do not change sol
		ForIndex(j,n) { sol[di][jpvt[j]]=b[j]; } EndFor;
	} EndFor;
}
