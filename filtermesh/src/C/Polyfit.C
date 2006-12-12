// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include <csignal>
#include "Options.h"
#include "A3dStream.h"
#include "FileIO.h"
#include "Bbox.h"
#include "FrameIO.h"
#include "LLS.h"
#include "Principal.h"
#include "Polygon.h"
#include "Set.h"
#include "Array.h"
#include "Map.h"
#include "Stack.h"
#include "Random.h"
#include "Timer.h"

#include "fdstream.hpp"
#include <iomanip>
#include <iostream>
using std::cout;
using std::ostream;
using std::flush;

// Signals:
//  SIGUSR1 : dump current state to file and continue

struct mvertex; typedef mvertex *vertex;
struct mvertex {
	Point p;		// position of vertex
	vertex v[2];		// v[0]: previous vertex
				// v[1]: next vertex
				// either may be 0
	// vertex has an associated edge if v[1]!=0
	Set<int> pts;		// points projecting onto associated edge
};
static Set<vertex> verts;
static int closed;		// polyline closed?

static struct _pt {
	int n;
	Array<Point> co;	// position of point
	Array<vertex> cle;	// closest edge
	Array<double> dis2;	// distance squared to closest edge
	_pt() : n(0), co(0), cle(0), dis2(0) {} DISABLECOPY(_pt); // GNUG
} pt;

static Set<vertex> ecand;	// set of candidate edges in stoc

static double spring=0;
static double fliter=1;
static double crep=1e-5;
static int nooutput=0;
static int verb=1;

static WSA3dStream oa3d(cout);
static int sdebug=GetenvValue("DEBUG");
static Frame xform;		// original verts+pts -> verts+pts in unit cube
static Frame xformi;		// inverse
enum Operation { OP_ecol, OP_espl };
const char* opname[2]={"ecol", "espl"};
enum OResult { OR_success, OR_energy, OR_illegal };
const char* orname[3]={"success", "positive_energy", "illegal"};
static int gotsignalUSR1=0;
static struct {
	int na[2], ns[2], nor[3];
} opstat;

const double sprsche[]={1e-2,1e-3,1e-4,1e-8,-1};
const int MAXGFITITER=30;
static FILE* filespawn=0;
static ostream* ospawn=0;
static WSA3dStream* a3dspawn;

// Declarations
static int dofilename(int argc, char const** argv);
static int dopfilename(int argc, char const** argv);
static int dosample(int argc, char const** argv);
static int doopencurve(int argc, char const** argv);
static int doclosedcurve(int argc, char const** argv);
static void perhapsinitialize();
static void initialprojection();
static void doreconstruct();
static void dosimplify();
static int dogfit(int argc, char const** argv);
static void dostoc();
static int dolfit(int argc, char const** argv);
static int dooutpoly(int argc, char const** argv);
static int dospawn(int argc, char const** argv);
static void outputpoly(WSA3dStream& oa3d, int clearobject=0);
// stoc
static OResult tryop(vertex v, Operation op, double& edrss);
static OResult tryecol(vertex v, int ni, int nri, double& edrss);
static OResult tryespl(vertex v, int ni, int nri, double& edrss);
// other
static void analyzepoly(int indent, const char* s);
static double getedis();
static double getespr();
static double geterep();
static void checksignal();
// unordered

static void sigUSR1()
{
	signal(SIGUSR1,HHSIG_PF(sigUSR1)); // for ATT unix
	gotsignalUSR1=1;
}

int main(int argc, char const** argv)
{
	SHOWDF("Created using\n%s",CreateHeader(argc,argv));
	Options opts;
	OPTSDC(opts,pfilename,	"name : initial polygon (can be -)");
	OPTSDC(opts,sample,	"n : sample n points from polygon");
	OPTSDC(opts,filename,	"name : point data (can be -)");
	OPTSDC(opts,opencurve,	"n : create initial polyline from points");
	OPTSDC(opts,closedcurve,"n : create initial polygon from points");
	OPTSPC(opts,crep,	"v: set constant for representation energy");
	OPTSDC(opts,reconstruct,": apply reconstruction schedule");
	OPTSDC(opts,simplify,	": apply simplification schedule");
	opts.c("",":");
	OPTSPC(opts,spring,	"tension : set sprint constant");
	opts.c("",":");
	OPTSDC(opts,gfit,	"niter : do global fit (0=until convergence)");
	OPTSDC(opts,stoc,	"do stochastic operations");
	OPTSDC(opts,lfit,	"ni nli : do ni iters, each nli local fits");
	OPTSDC(opts,outpoly,	"filename : output current poly to file");
	OPTSDC(opts,spawn,	"'command': send record to popen");
	OPTSFC(opts,nooutput,	": don't print final poly on stdout");
	OPTSPC(opts,fliter,	"factor : modify # local iters done in stoc");
	OPTSPC(opts,verb,	"i : verbosity level (1=avg,2=more,3=lots)");
	signal(SIGUSR1,HHSIG_PF(sigUSR1));
	NEST {
		TIMER(Polyfit);
		opts.allmustparse();
		if (!opts.parse(argc,argv)) {
			opts.problem(argc,argv);
			return 1;
		}
		perhapsinitialize();
		SHOWDF("\n"); analyzepoly(0,"FINAL"); SHOWDF("\n");
	}
	SHOWDF("\n"),SHOWDF("EndPolyfit\n");
	if (filespawn) {
		delete a3dspawn;
		delete ospawn;
		pclose(filespawn);
	}
	CleanUp();
	if (!nooutput) outputpoly(oa3d);
	ForSet(verts,vertex,v) { delete v; } EndFor;
	verts.clear();
	return 0;
}

static void enterpoint(const Point& p)
{
	pt.co += p;
	pt.cle += static_cast<vertex>(0);
	pt.dis2 += static_cast<double>(0);
	pt.n++;
}

static int dofilename(int argc, char const** argv)
{
	TIMER(_filename);
	assertx(argc>1);
	assertx(!pt.n);
	RFile is(argv[1]);
	RSA3dStream ia3d(is());
	A3dElem el;
	for (;;) {
		ia3d.read(el);
		if (el.type()==A3dElem::TEndFile) break;
		if (el.type()==A3dElem::TComment) continue;
		if (assertw1(el.type()==A3dElem::TPoint)) continue;
		enterpoint(el[0].p);
	}
	SHOWDF("%d points read\n",pt.n);
	return 2;
}

static void initializepoly(const Polygon& poly)
{
	assertx(verts.empty());
	int num=poly.num();
	closed=!compare(poly[0],poly[num-1]);
	if (closed) num--;
	vertex vf=0, vl=0;
	for (int i=0;i<num;i++) {
		vertex v=new mvertex;
		verts.enter(v);
		v->p=poly[i];
		v->v[0]=vl;
		if (vl) vl->v[1]=v;
		vl=v;
		if (!vf) vf=v;
	}
	vf->v[0]=closed?vl:0;
	vl->v[1]=closed?vf:0;
	SHOWDF("created %s poly structure with %d vertices\n",
	       closed?"closed":"open",verts.num());
}

static int dopfilename(int argc, char const** argv)
{
	TIMER(_pfilename);
	assertx(argc>1);
	RFile is(argv[1]);
	RSA3dStream ia3d(is());
	A3dElem el;
	for (;;) {
		ia3d.read(el);
		if (el.type()==A3dElem::TEndFile) break;
		if (el.type()==A3dElem::TComment) continue;
		if (assertw1(el.type()==A3dElem::TPolyline)) continue;
		el.update(A3dElem::TPolygon);
		initializepoly(el);
	}
	assertx(verts.num());
	return 2;
}

static int dosample(int argc, char const** argv)
{
        vertex v;
	assertx(argc>1);
	assertx(!pt.n && verts.num());
	int n=atoi(argv[1]);
	// *Non-uniformly* sample points from polygon ??
	for (int i=0;i<n;i++) {
		for (;;) {
			SetIter<vertex> seti(verts,Random::G);
			v=seti();
			if (v->v[1]) break;
		}
		enterpoint(interp(v->p,v->v[1]->p,Random::G.unif()));
	}
	// Enter vertices as points
	ForSet(verts,vertex,v) {
		enterpoint(v->p);
	} EndFor;
	SHOWDF("%d points read\n",pt.n);
	return 2;
}

static void createpoly(int pclosed, int n)
{
	// initialize the manifold using principal components of data
	assertx(pt.n);
	Frame frame; double eimag[3];
	PrincipalComponents(pt.co,pt.n,frame,eimag);
	Polygon poly;
	if (!pclosed) {
		assertx(n>1);
		double stdv=1.1;
		// generate n points along principal line (-1,1)*stdv
		for (int i=0;i<n;i++)
			poly+=Point(-stdv+2*stdv*i/(n-1.),0,0)*frame;
	} else {
		assertx(n>2);
		double stdv=1.5;
		// generate n points along ellipse at radius stdv
		for (int i=0;i<=n;i++) {
			double a=double(i)/n*2*PI;
			poly+=Point(stdv*cos(a),stdv*sin(a),0)*frame;
		}
	}
	initializepoly(poly);
}

static int doopencurve(int argc, char const** argv)
{
	assertx(argc>1);
	createpoly(0,atoi(argv[1]));
	return 2;
}

static int doclosedcurve(int argc, char const** argv)
{
	assertx(argc>1);
	createpoly(1,atoi(argv[1]));
	return 2;
}

static void polytransform(const Frame& f)
{
	ForSet(verts,vertex,v) { v->p*=f; } EndFor;
}

static void computexform()
{
	Bbox bb;
	bb.clear();
	NEST { for (int i=0;i<pt.n;i++) bb.takeunion(pt.co[i]); }
	ForSet(verts,vertex,v) { bb.takeunion(v->p); } EndFor;
	xform=bb.getFrameToSmallCube();
	SHOWDF("Applying xform: %s",FrameIO::string(xform,1,0));
	xformi=~xform;
	NEST { for (int i=0;i<pt.n;i++) pt.co[i]*=xform; }
	polytransform(xform);
}

static void perhapsinitialize()
{
	static int done=0; if (done++) return;
	assertx(pt.n && verts.num());
	assertw1(spring>0);	// just warn user
	computexform();
	initialprojection();
}

/* Return squared distance from p to segment (p1,p2).
Optionally return barycentric coordinate with respect to p1
p==interp(p1,p2,bary)   (bary may lie outside [0,1]) */
static double distpe2(const Point& p, const Point& p1,
		     const Point& p2, double* bary)
{
	Vector v12=p2-p1, v1p=p-p1;
	double d122=mag2(v12);
	if (d122<1e-12) { if (bary) *bary=.5; return dist2(p,p1); }
	double d12=sqrt(d122);
	double don12=dot(v12,v1p)/d12;
	if (bary) *bary=1-don12/d12;
	if (don12<0) return dist2(p,p1);
	if (don12>d12) return dist2(p,p2);
	double d1p2=mag2(v1p);
	if (d1p2<1e-12) { if (bary) *bary=1; return dist2(p,p1); }
	return max(d1p2-don12*don12,0.);
}

static void initialprojection()
{
	TIMER(_initialproj);
	// ?? do it inefficiently for now
	for (int i=0;i<pt.n;i++) {
		double mind2=1e30;
		vertex mine=0;
		ForSet(verts,vertex,v) {
			if (!v->v[1]) continue;
			double d2=distpe2(pt.co[i],v->p,v->v[1]->p,0);
			if (d2<mind2) mind2=d2,mine=v;
		} EndFor;
		pt.cle[i]=assertv(mine);
		pt.dis2[i]=mind2;
		mine->pts.enter(i);
	}
	analyzepoly(0,"INITIAL");
}

static void reprojectlocally(int pi)
{
	vertex v=assertv(pt.cle[pi]);
	assertx(v->pts.contains(pi)); // optional
	assertx(v->v[1]);	      // optional
	double a, mind2=distpe2(pt.co[pi],v->p,v->v[1]->p,0);
	vertex mine=v;
	if (v->v[0] &&
	    (a=distpe2(pt.co[pi],v->v[0]->p,v->p,0))<mind2)
		mind2=a,mine=v->v[0];
	if (v->v[1]->v[1] &&
	    (a=distpe2(pt.co[pi],v->v[1]->p,v->v[1]->v[1]->p,0))<mind2)
		mind2=a,mine=v->v[1];
	pt.dis2[pi]=mind2;
	if (mine==v) return;
	assertx(v->pts.remove(pi));
	mine->pts.enter(pi);
	pt.cle[pi]=mine;
}

static void globalproject()
{
	// local projection
	for (int i=0;i<pt.n;i++) reprojectlocally(i);
}

static void globalfit()
{
	Map<vertex,int> mvi;
	Array<vertex> va;
	ForSet(verts,vertex,v) {
		mvi.enter(v,va.num());
		va+=v;
	} EndFor;
	int m=pt.n, n=verts.num();
	if (spring) m+=verts.num()-!closed;
	SHOWF("GlobalFit: about to solve a %dx%d LLS system\n",m,n);
	LLS* lls=LLS::create(m,n,3,2./n);
	// Add point constraints
	For (int i=0;i<pt.n;i++) {
		vertex cle=assertv(pt.cle[i]);
		vertex v[2]; v[0]=cle; v[1]=assertv(cle->v[1]);
		double bary;
		(void)distpe2(pt.co[i],v[0]->p,v[1]->p,&bary);
		bary=min(max(bary,0.),1.);
		for (int j=0;j<2;j++)
			lls->enter(i,mvi.get(v[j]),!j?bary:1-bary);
		lls->enterRh(i,&pt.co[i][0]);
	} EndFor;
	// Add spring constraints
	if (spring) {
		double sqrtit=mysqrt(spring);
		Vector vzero(0,0,0);
		int ri=pt.n;
		ForSet(verts,vertex,v) {
			if (!v->v[1]) continue;
			lls->enter(ri,mvi.get(v),sqrtit);
			lls->enter(ri,mvi.get(v->v[1]),-sqrtit);
			lls->enterRh(ri,&vzero[0]);
			ri++;
		} EndFor;
		assertx(ri==m);
	}
	// Suggest current solution
	int i;
	for (i=0;i<n;i++) lls->enterEstimate(i,&va[i]->p[0]);
	// Solve
	lls->solve();
	// Update solution
	for (i=0;i<n;i++) {
		Point p; lls->getValue(i,&p[0]);
		va[i]->p=p;
	}
	delete lls,lls=0;
}

// v[0] or v[1] may be 0
static void localfit(const Stack<int>& stpts, vertex v[2], int niter,
		     Point& newp, double& prss0, double& prss1)
{
        int c;
	assertx(v[0] || v[1]);
	double sqrtit=mysqrt(spring);
	for (int ni=0;ni<niter;ni++) {
		double UtU[3], Utb[3], btb[3], rss0=0, rss1=0;
//		for (int c=0;c<3;c++) UtU[c]=Utb[c]=btb[c]=0;
		for (c=0;c<3;c++) UtU[c]=Utb[c]=btb[c]=0;
		// Enter projections
		ForStack(stpts,int,pi) {
			const Point& p=pt.co[pi];
			double bary0, bary1;
			double d0=v[0]?distpe2(p,v[0]->p,newp,&bary0):1e30;
			double d1=v[1]?distpe2(p,v[1]->p,newp,&bary1):1e30;
			int mini=d0<d1?0:1; assertx(v[mini]);
			double bary=mini?bary1:bary0;
			bary=min(max(bary,0.),1.);
			double u=1-bary;
			for (int c=0;c<3;c++) {
				double b=p[c]-bary*v[mini]->p[c];
				UtU[c]+=u*u; Utb[c]+=u*b; btb[c]+=b*b;
				rss0+=square(u*newp[c]-b);
			}
		} EndFor;
		// Enter springs
		for (int i=0;i<2;i++) {
			if (!v[i]) continue;
			for (int c=0;c<3;c++) {
				double u=sqrtit, b=double(v[i]->p[c])*sqrtit;
				UtU[c]+=u*u; Utb[c]+=u*b; btb[c]+=b*b;
				rss0+=square(u*newp[c]-b);
			}
		}
		// Solve
		for (c=0;c<3;c++) {
			double newv=assertw1(UtU[c])?newp[c]:Utb[c]/UtU[c];
			newp[c]=newv;
			double a=btb[c]-UtU[c]*square(newv);
			assertw1(a>-1e-8);
			if (a>0) rss1+=a;
		}
		assertw1(rss1-rss0<1e-13);
		if (!ni) prss0=rss0;
		prss1=rss1;
	}
}

static void fitring(vertex v, int niter)
{
	assertx(v);
	double rss0, rss1;
	vertex va[2]; va[0]=v->v[0]; va[1]=v->v[1];
	Stack<int> stpts;
	if (va[0]) ForSet(va[0]->pts,int,pi) { stpts.push(pi); } EndFor;
	if (va[1]) ForSet(v->pts,int,pi) { stpts.push(pi); } EndFor;
	localfit(stpts,va,niter,v->p,rss0,rss1);
	ForStack(stpts,int,pi) { reprojectlocally(pi); } EndFor;
}

static void cleanupneighborhood(vertex v, int nri)
{
	assertx(v);
	if (!nri) return;
	vertex v0=v->v[0], v00=v0 && v0->v[0]?v0->v[0]:0;
	vertex v1=v->v[1], v11=v1 && v1->v[1]?v1->v[1]:0;
	fitring(v,nri);
	if (v0) fitring(v0,nri);
	if (v1) fitring(v1,nri);
	if (v00) fitring(v00,nri);
	if (v11) fitring(v11,nri);
	if (v0) fitring(v0,nri);
	if (v1) fitring(v1,nri);
	fitring(v,nri);
	if (v0) fitring(v0,nri);
	if (v1) fitring(v1,nri);
	fitring(v,nri);
}

static void applyschedule()
{
	const char* argv[3]; argv[1]="2"; argv[2]="3";
	for (;spring>sprsche[0];) {
		dolfit(3, argv);	// -lfit 2 3
		dostoc();		// -stoc
		dolfit(3, argv);	// -lfit 2 3
		spring*=.1;
	}
	for (int i=0;sprsche[i]>=0;i++) {
		spring=sprsche[i]; // -spring f
		dolfit(3, argv);	// -lfit 2 3
		dostoc();		// -stoc
		dolfit(3, argv);	// -lfit 2 3
	}
}

static void doreconstruct()
{
	TIMER(_reconstruct);
	if (!spring) spring=sprsche[0];
	perhapsinitialize();
	const char* argv[2]; argv[1]="0";
	dogfit(2, argv);
	applyschedule();
}

static void dosimplify()
{
	TIMER(_simplify);
	if (!spring) spring=sprsche[0];
	perhapsinitialize();
	applyschedule();
}

static int dogfit(int argc, char const** argv)
{
        int i;
	perhapsinitialize();
	TIMER(_gfit);
	assertx(argc>1);
	int niter=atoi(argv[1]);
	SHOWDF("\n");
	SHOWDF("Beginning gfit, %d iterations, spr=%g\n",niter,spring);
	double ecsc=getedis()+getespr();	// energy constant simplicial complex
//	for (int i=0;!niter || i<niter;) {
	for (i=0;!niter || i<niter;) {
		if (!niter && i>=MAXGFITITER) break;
		i++;
		SHOWDF("iter %d/%d\n",i,niter);
		cout << flush;
		NEST { ATIMER(__lls); globalfit(); }
		NEST { ATIMER(__project); globalproject(); }
		double necsc=getedis()+getespr();
		double echange=necsc-ecsc;
		assertw(echange<0);
		if (verb>=2) {
			analyzepoly(2,"gfit_iter");
			SHOWDF(" change in energy=%g\n",echange);
		}
		if (a3dspawn) outputpoly(*a3dspawn,1);
		if (!niter && echange>-1e-4) break;
		ecsc=necsc;
	}
	SHOWDF("Finished gfit, did %d iterations\n",i);
	analyzepoly(0,"after_gfit");
	return 2;
}

static void dostoc()
{
	perhapsinitialize();
	TIMER(_stoc);
	SHOWDF("\n");
	SHOWDF("Beginning stoc, spring=%g, fliter=%g\n",spring,fliter);
	NEST { for (int c=0;c<2;c++) opstat.na[c]=opstat.ns[c]=0; }
	NEST { for (int c=0;c<3;c++) opstat.nor[c]=0; }
	int i=0, nbad=0, lasti=-99999;
	ForSet(verts,vertex,v) {
		if (v->v[1]) ecand.enter(v);
	} EndFor;
	for (;!ecand.empty();) {
		i++;
		SetIter<vertex> si(ecand,Random::G);
		vertex v=si();
		assertx(ecand.remove(v));
		assertx(v->v[1]);
		Operation op; OResult oor=OR_illegal; double edrss;
		if (oor!=OR_success) op=OP_ecol,oor=tryop(v,op,edrss);
		if (oor!=OR_success) op=OP_espl,oor=tryop(v,op,edrss);
		if (oor!=OR_success) { nbad++; continue; }
		if (verb>=2 || verb>=1 && i>=lasti+100) {
			SHOWDF("it %5d, %s (after %3d) [%5d/%-5d] edrss=%e\n",
			       i,opname[op],nbad,
			       ecand.num(),verts.num(),edrss);
			lasti=i;
		}
		if (a3dspawn) outputpoly(*a3dspawn,1);
		nbad=0;
	}
	SHOWDF("it %d, last search: %d wasted attempts\n",i,nbad);
	int nat=opstat.na[0]+opstat.na[1];
	int nst=opstat.ns[0]+opstat.ns[1];
	SHOWDF("Endstoc:  (col=%d/%d, espl=%d/%d tot=%d/%d)\n",
	       opstat.ns[OP_ecol],opstat.na[OP_ecol],
	       opstat.ns[OP_espl],opstat.na[OP_espl],
	       nst,nat);
	SHOWDF("Result of %d attempted operations:\n",nat);
	ForIndex(i,3) {
		SHOWDF("  %5d %s\n",opstat.nor[i],orname[i]);
	} EndFor;
	analyzepoly(0,"after_stoc");
}

static int dolfit(int argc, char const** argv)
{
	perhapsinitialize();
	TIMER(_lfit);
	assertx(argc>2);
	int ni=atoi(argv[1]);
	int nli=atoi(argv[2]);
	SHOWDF("\n");
	SHOWDF("Beginning lfit, %d iters (nli=%d), spr=%g\n",ni,nli,spring);
	for (int i=0;i<ni;i++) {
		ForSet(verts,vertex,v) {
			checksignal();
			fitring(v,nli);
		} EndFor;
		if (a3dspawn) outputpoly(*a3dspawn,1);
	}
	SHOWDF("Finished lfit\n");
	analyzepoly(0,"after_lfit");
	return 3;
}

static void outputpoly(WSA3dStream& oa3d, int clearobject)
{
	perhapsinitialize();
	polytransform(xformi);
	if (clearobject) oa3d.writeClearObject();
	A3dElem el;
	el.init(A3dElem::TPolyline);
	vertex v=verts.getone();
	if (!closed) while (v->v[0]) v=v->v[0];
	for (vertex vf=v;v;) {
		el+=A3dVertex(v->p,Vector(0,0,0),A3dVertexColor::White);
		if (v==vf && el.num()>1) break;
		v=v->v[1];
	}
	oa3d.write(el);
	oa3d.os() << flush;
	polytransform(xform);
}

static int dooutpoly(int argc, char const** argv)
{
	TIMER(_outpoly);
	assertx(argc>1);
	WFile os(argv[1]);
	WSA3dStream oa3d(os());
	outputpoly(oa3d);
	return 2;
}

static int dospawn(int argc, char const** argv)
{
	filespawn=assertv(popen(argv[1],"w"));
	ospawn=new boost::fdostream(fileno(filespawn));
	a3dspawn=new WSA3dStream(*ospawn);
	outputpoly(*a3dspawn);
	return 2;
}

//*** stoc

static OResult tryop(vertex v, Operation op, double& edrss)
{
	checksignal();
	ATIMER(__tryop);
	OResult oor;
	oor=(op==OP_ecol?tryecol(v,int(4*fliter+.5),int(2*fliter+.5),edrss):
	     op==OP_espl?tryespl(v,int(3*fliter+.5),int(4*fliter+.5),edrss):
	     (assertnever(""),OR_success));
	opstat.na[op]++;
	if (oor==OR_success) opstat.ns[op]++;
	opstat.nor[oor]++;
	return oor;
}

//    v0      v      v1
//    va[0]                 va[1]
//    ev[0]   ev[1]  ev[2]
static OResult tryecol(vertex v, int ni, int nri, double& edrss)
{
	vertex v0=v->v[0], v1=assertv(v->v[1]), va[2], ev[3], ov=v0?v0:v1;
	va[0]=v0; va[1]=v1->v[1];
	ev[0]=v0; ev[1]=v; ev[2]=v1;
	// closed triangle or single open segment
	if (va[0]==va[1]) return OR_illegal;
	double rssf=0;
	Stack<int> stpts;
	For (int i=0;i<3;i++) {
		if (!ev[i] || !ev[i]->v[1]) continue;
		ForSet(ev[i]->pts,int,pi) { stpts.push(pi); } EndFor;
		rssf+=spring*dist2(ev[i]->p,ev[i]->v[1]->p);
	} EndFor;
	ForStack(stpts,int,pi) { rssf+=pt.dis2[pi]; } EndFor;
	// Find the best starting location by exploring one iteration.
	double minrss1=1e30; int minii=-1;
	for (int ii=0;ii<3;ii++) {
		Point newp=interp(v->p,v1->p,ii*.5);
		double rss0,rss1;
		localfit(stpts,va,1,newp,rss0,rss1);
		if (rss1<minrss1) minrss1=rss1,minii=ii;
	}
	if (minii<0) assertnever("");
	// Then, explore ni iterations from that chosen starting point
	Point newp=interp(v->p,v1->p,minii*.5);
	double rss0,rss1;
	localfit(stpts,va,ni,newp,rss0,rss1);
	double drss=rss1-rssf-double(crep);
	edrss=drss;
	if (verb>=3) SHOWF("# ecol: rssf=%lg rss1=%lg drss=%lg\n",
			   rssf,rss1,drss);
	if (drss>=0) return OR_energy; // energy function does not decrease
	// ALL SYSTEMS GO
	// move points off to other segment and reproject later
	ForSet(v->pts,int,pi) {
		ov->pts.enter(pi);
		pt.cle[pi]=ov;
	} EndFor;
	v->pts.clear();
	delete v;
	assertx(verts.remove(v));
	if (v0) v0->v[1]=v1;
	v1->v[0]=v0;
	if (v0) ecand.add(v0);
	if (v1->v[1]) ecand.add(v1);
	v1->p=newp;
	cleanupneighborhood(v1,nri);
	return OR_success;
}

static OResult tryespl(vertex v, int ni, int nri, double& edrss)
{
	// always legal
	vertex va[2]; va[0]=v; va[1]=assertv(v->v[1]);
	double rssf=spring*dist2(va[0]->p,va[1]->p);
	Stack<int> stpts;
	ForSet(v->pts,int,pi) { stpts.push(pi); } EndFor;
	ForStack(stpts,int,pi) { rssf+=pt.dis2[pi]; } EndFor;
	Point newp=interp(va[0]->p,va[1]->p);
	double rss0,rss1;
	localfit(stpts,va,ni,newp,rss0,rss1);
	double drss=rss1-rssf+double(crep);
	edrss=drss;
	if (verb>=3) SHOWF("# espl: rssf=%lg rss1=%lg drss=%lg\n",
			   rssf,rss1,drss);
	if (drss>=0) return OR_energy; // energy function does not decrease
	// ALL SYSTEMS GO
	vertex vn=new mvertex;
	verts.enter(vn);
	vn->p=newp;
	vn->v[0]=va[0];
	vn->v[1]=va[1];
	va[0]->v[1]=vn;
	va[1]->v[0]=vn;
	if (va[0]->v[0]) ecand.add(va[0]->v[0]);
	ecand.add(va[0]);
	ecand.add(vn);
	if (va[1]->v[1]) ecand.add(va[1]);
	// v->pts not cleared, let refit reproject points
	cleanupneighborhood(vn,nri);
	return OR_success;
}

//*** other

static void analyzepoly(int indent, const char* s)
{
	char g[100];
	g[0]=0; for (int i=0;i<indent;i++) strcat(g," ");
	int nv=verts.num();
	double edis=getedis(), espr=getespr(), erep=geterep();
	erep=crep*nv;
	double etot=edis+espr+erep;
	SHOWDF("%sPoly analysis: %s\n",g,s);
	SHOWDF("%s  poly: v=%d\n",g,nv);
	SHOWDF("%s  parameters: crep=%g spring=%g\n",g,crep,spring);
	SHOWDF("%s  energies: edis=%g espr=%g erep=%g etot=%g\n",
	       g,edis,espr,erep,etot);
}

static double getedis()
{
	double edis=0;
	for (int i=0;i<pt.n;i++) edis+=pt.dis2[i];
	return edis;
}

static double getespr()
{
	double espr=0;
	ForSet(verts,vertex,v) {
		if (v->v[1]) espr+=spring*dist2(v->p,v->v[1]->p);
	} EndFor;
	return espr;
}

static double geterep()
{
	return crep*verts.num();
}

static void checksignal()
{
	if (!gotsignalUSR1) return;
	gotsignalUSR1=0;
	SHOWS("***starting outpoly");
	NEST {
		const char* argv[2]; argv[1]="v.polyfit.m.Z";
		dooutpoly(2, argv);
	}
	SHOWS("***done with outpoly");
}
