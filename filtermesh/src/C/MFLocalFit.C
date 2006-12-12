// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Meshfit.h"
#include "Facedistance.h"
#include "GeomOp.h"
#include "Set.h"
#include "Array.h"
#include "Stat.h"

//*** UPointLLS

/*
 UPointLLS:
 Solve a linear least squares problem involving a single point,
 in which the problem decomposes into 3 independent univariate LLS problems.
 The system is Ux=b,
 the design matrix U is a column vector with m rows
 (the m constraints coming from either the springs or the projected points).
 The unknown x is simply a scalar (actually, one scalar for each dimension).
 The vector b is also a column vector with m rows.
 To solve it using (UtU)x=Utb, note that UtU is simply the norm of U,
 and Utb is simply the dot product of U and b.
 To do this efficiently, UtU and Utb can be accumulated for all 3
 coordinates simultaneously, while traversing U and b row-by-row.
 To compute the rss ( ||Ux-b||^2 ), Werner observed that
 rss= ||b||^2-||Ux||^2 = ||b||^2 - x^2*||U||^2.
 */
class UPointLLS {
  public:
	UPointLLS(Point& pp);
	~UPointLLS() { }
	void enterspring(const Point& pother, double sqrttension);
	// Constraint between point pdata and the point on triangle (pp,p1,p2)
	// with barycentric coordinates (1-param1-param2,param1,param2)
	void enterprojection(const Point& pdata,
			     const Point& p1, const Point& p2,
			     double param1, double param2);
	void solve(double* rss0, double* rss1); // updates point pp!
  private:
	Point& p;
	double UtU[3], Utb[3], btb[3], rss0;
};

UPointLLS::UPointLLS(Point& pp) : p(pp), rss0(0)
{
	for (int c=0;c<3;c++) UtU[c]=Utb[c]=btb[c]=0;
}

inline void UPointLLS::enterspring(const Point& pother, double sqrttension)
{
	for (int c=0;c<3;c++) {
		double u=sqrttension;
		double b=double(pother[c])*sqrttension;
		UtU[c]+=u*u; Utb[c]+=u*b; btb[c]+=b*b;
		rss0+=square(u*p[c]-b);
	}
}

inline void UPointLLS::enterprojection(const Point& pdata,
				       const Point& p1, const Point& p2,
				       double param1, double param2)
{
	double u=1-param1-param2;
	double pa1=param1, pa2=param2;
	for (int c=0;c<3;c++) {
		double b=pdata[c]-pa1*p1[c]-pa2*p2[c];
		UtU[c]+=u*u; Utb[c]+=u*b; btb[c]+=b*b;
		rss0+=square(u*p[c]-b);
	}
}

void UPointLLS::solve(double* prss0, double* prss1)
{
	double rss1=0;
	for (int c=0;c<3;c++) {
		double newv=assertw1(UtU[c]) ? p[c] : Utb[c]/UtU[c];
		p[c]=newv;
		double a=btb[c]-UtU[c]*square(newv);
		assertw1(a>-1e-8);
		if (a>0) rss1+=a;
	}
	assertw1(rss1-rss0<1e-13);
	if (prss0) *prss0=rss0;
	if (prss1) *prss1=rss1;
}

//*** other

void GatherVertexRing(const GMesh& mesh, Vertex v, Array<const Point*>& wa)
{
	wa.init(0);
	Vertex w=mesh.mostClwVertex(v), wf=w;
	for (;;) {
		wa+=&mesh.point(w);
		w=mesh.ccwVertex(v,w);
		if (!w || w==wf) break;
	}
	if (w) wa+=wa[0];
}

void GatherEdgeRing(const GMesh& mesh, Edge e, Array<const Point*>& wa)
{
	Vertex v1=mesh.vertex1(e), v2=mesh.vertex2(e);
	// current Mesh implementation boundary Edge direction
	if (mesh.isBoundary(e)) assertx(mesh.mostClwVertex(v2)!=v1);
	Vertex cv=v2;		// current vertex for rotation
	if (!mesh.isBoundary(e) && mesh.isBoundary(v1)) cv=v1;
	Vertex ov=mesh.oppVertex(cv,e); // other vertex of rotation (v1 or v2)
	Vertex w=mesh.mostClwVertex(cv);
	if (w==ov) w=mesh.clwVertex(cv,w),assertx(!mesh.isBoundary(cv));
	Vertex wf=assertv(w);
	wa.init(0);
	for (;;) {
		wa+=&mesh.point(w);
		Vertex w2=mesh.ccwVertex(cv,w);
		if (w2==ov) {
			ov=cv; cv=w2;
			w2=mesh.ccwVertex(cv,w);
		}
		w=w2;
		if (!w || w==wf) break;
	}
	if (w==wf) wa+=wa[0];
	assertx(wa.num()>=(wa[0]==wa[wa.num()-1]?4:2));
}

// Fit a set of points to a ring of vertices while optimizing the center
// vertex.  Does niter iterations of projection+refit.
// Afterwards, should call ReprojectLocally() or equivalent to do final
// reprojection and update pt projections.
//  * Given:
// stpts: the set of point indices that project locally
// wa: array of vertex positions (wa[0]==wa[nw-1] if closed loop)
// niter: number of iterations to do
// newp: the initial position to use for center vertex
//  * Return:
// newp: the final fitted position
// rss0: energy after first projection (true)
// rss1: energy after final refit, before final reprojection (over-estimate)
void LocalFit(const DataPts& pt, double spring, double spbf,
	      const Stack<int>& stpts, const Array<const Point*>& wa,
	      int niter, Point& newp, double& prss0, double& prss1)
{
	int nw=wa.num();
	assertx(nw>1 && niter>0); // at least one face
	int closed=wa[0]==wa[nw-1];
	double sqrtit=mysqrt(spring), sqrtbt=mysqrt(spring*spbf);
	double rss1;
	for (int ni=0;ni<niter;ni++) {
		UPointLLS ulls(newp);
		ForStack(stpts,int,pi) {
			const Point& p=pt.co[pi];
			static Pqueue<int> pq;
			pq.clear();
			For (int i=0;i<nw-1;i++) {
				double d=LBDistPTri(p,newp,*wa[i],*wa[i+1]);
				pq.enterUnsorted(i,d*d);
			} EndFor;
			pq.sort();
			SSTAT(SLFconsid,nw-1);
			int nproj=0; double mind2=1e30; int mini=-1; Bary minb;
			for (;!pq.empty();) {
				if (pq.minpriority()>=mind2) break;
				int i=pq.removemin();
				double d2; Bary b;
				ProjectPTri(p,newp,*wa[i],*wa[i+1],&d2,0,&b,0);
				nproj++;
				if (d2<mind2) { mind2=d2; mini=i; minb=b; }
			}
			SSTAT(SLFproj,nproj);
			// Found closest face mini and corresponding minb
			assertx(mini>=0);
			ulls.enterprojection(p,*wa[mini],*wa[mini+1],
					     minb[1],minb[2]);
		} EndFor;
		for (int i=0;i<nw-closed;i++) {
			int isbe=!closed && (i==0 || i==nw-1);
			ulls.enterspring(*wa[i],isbe?sqrtbt:sqrtit);
		}
		double rss0;
		ulls.solve(&rss0,&rss1);
		if (!ni) prss0=rss0;
	}
	prss1=rss1;
}

double MinLocalDihedral(const Array<const Point*>& wa, const Point& newp)
{
	int nw=wa.num();
	assertx(nw>1);
	int open=wa[0]!=wa[nw-1];
	double mindic=2;
	for (int i=1;i<nw-open;i++) {
		int i1=i+1; if (i1==nw) i1=1;
		double dic=DihedralAngleCos(newp,*wa[i],*wa[i-1],*wa[i1]);
		mindic=min(mindic,dic);
	}
	return mindic;
}

double MinDihedralAboutVertex(const GMesh& mesh, Vertex v)
{
	Array<const Point*> wa;
	GatherVertexRing(mesh,v,wa);
	return MinLocalDihedral(wa,mesh.point(v));
}

void FitRing(GMesh& mesh, DataPts& pt, double spring, double spbf,
	     double mincos, Vertex v, int niter)
{
	assertx(niter>0);
	Stack<int> stpts;
	Stack<Face> stfaces;
	Array<const Point*> wa;
	GatherVertexRing(mesh,v,wa);
	ForVertexFace(mesh,v,f) {
		stfaces.push(f);
		ForSet(*pt.mfpts.get(f),int,pi) {
			stpts.push(pi);
		} EndFor;
	} EndFor;
	Point newp=mesh.point(v);
	double minb=MinLocalDihedral(wa,newp);
	double rss0,rss1;
	LocalFit(pt,spring,spbf,stpts,wa,niter,newp,rss0,rss1);
	double mina=MinLocalDihedral(wa,newp);
	if (mina<mincos && mina<minb) return; // change disallowed
	mesh.setPoint(v,newp);
	ReprojectLocally(mesh,pt,stpts,stfaces);
}
