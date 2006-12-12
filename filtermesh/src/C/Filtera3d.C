// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include <string>
#include "Hh.h"
#include "Options.h"
#include "FileIO.h"
#include "A3dStream.h"
#include "FrameIO.h"
#include "Stat.h"
#include "Bbox.h"
#include "Random.h"		// for -noise
#include "Polygon.h"		// for -intersect
#include "Kdtree.h"		// for -intersect
#include "HashPoint.h"		// for -joinlines
#include "Graph.h"		// for -joinlines
#include "Spatial.h"
#include "Stack.h"
#include "Array.h"

#include <iostream>
#include <sstream>

using std::string;
using std::cout;
using std::istringstream;

static WSA3dStream oa3d(cout);

static int nopolygons, nopolylines, nopoints, onlypoly;
static int nonormals, optnormals, nocolor;
static int fixdegen;
static int gnormalize, fixorient, flipnormals;
static double shownormals, stretch, noise, offset;
static int twosided, triangulate, tessellate, intersect;
static int statistics, box, boxframe, nooutput;
static int every, first, split;
static Point cusphc;
static double cusphr, mindis;
static double speedup, frdelay, eldelay;
static int ctoascii, ctobinary, ctoold;

static int ndegen=0;
static int isrestrictf;
static Frame crestrictf;
static int istransf;
static Frame ctransf;
static int isctransfinv;
static Frame ctransfinv;
static A3dColor cdiff, cspec, cphong;
static int nfixorient;		// # polygons flipped due to fixorient
static int joinlines;
static int ncullsphere;
static int nmindis;
static struct {
	int npolyb;
	int ntrib;
	int ntria;
} tri;
static struct {
	int ntrib;
	int ntria;
} tess;
static struct {
	Bbox bb;		// global bounding box of all polygons
	Stack<Polygon*> Spolyg;
	int nedges;
	Polygon* polyg;
} inter;
static struct _join {
	HashPoint* hp;
	Array<Point> pa;
	Graph<int>* graph;
	_join() : hp(0), pa(0), graph(0) {} DISABLECOPY(_join);	// GNUG
} join;
static STATNP(Slnvert);		// polyline # of vertices
static STATNP(Sledgel);		// polyline edge length
static STATNP(Spnvert);		// polygon # of vertices
static STATNP(Spedgel);		// polygon edge length
static STATNP(Sparea);		// polygon area
static STATNP(Splanar);		// polygon planarity (0=planar)
static STATNP(Sptnor);		// point, existence of normal
static Bbox bbox;		// box extent
static double fsplit;		// fsplit=split; { fsplit*=speedup; }
static A3dVertexColor zerovertexcolor;

// Declarations
static void doit(RSA3dStream& ia3d);
static int loop(A3dElem& el);
static void pass1(const A3dElem& el);
static void pass2(const A3dElem& el);
static void pass3(const A3dElem& el);
static void doout(const A3dElem& el);
static void delayframe();
static void delayelement();
static int isdegen(const A3dElem& el);
static int outofbounds(const A3dElem& el);
static int polygonneedsflip(const A3dElem& el);
static void flippolygon(A3dElem& el);
static void doshownormals(const A3dElem& el);
static void dostat(const A3dElem& el);
static void dointersect();
static int domindis(const Point& p);
static void dojoinlines();

int main(int argc, char const** argv)
{
	double diff[3],spec[3],phong[3]; diff[0]=spec[0]=phong[0]=-1; // undef
	double cullsphere[4]; cullsphere[3]=0;
	const char* restrictf=0;
	const char* transf=0;
	Options opts;
	
	OPTSFC(opts,onlypoly,	": remove special a3d commands");
	OPTSFC(opts,nopolygons,	": cull polygons");
	OPTSFC(opts,nopolylines,": cull polylines");
	OPTSFC(opts,nopoints,	": cull points");
	OPTSFC(opts,fixdegen,	": remove zero area triangles");
	OPTSPC(opts,restrictf,	"'frame' : cull elems outside unit frame");
	OPTSPC(opts,every,	"i : use only every ith element");
	OPTSPC(opts,first,	"i : use only first i elements");
	OPTSFC(opts,joinlines,	": join line segments into polyline");
	opts.c("",":");
	OPTSAC(opts,cullsphere,4,"x y z r : remove points within sphere");
	OPTSPC(opts,mindis,	"f : make no pair of points closer than f");
	opts.c("",":");
	OPTSFC(opts,nonormals,	": remove vertex normals");
	OPTSFC(opts,optnormals,	": remove unnecessary polygon normals");
	OPTSFC(opts,nocolor,	": remove color information");
	OPTSAC(opts,diff,3,	"r g b : set diffuse color");
	OPTSAC(opts,spec,3,	"r g b : set specular color");
	OPTSAC(opts,phong,3,	"r g b : set phong color");
	OPTSPC(opts,transf,	"'frame' : transform all elements");
	opts.f("-normalize",gnormalize,": normalize normals");
	OPTSFC(opts,fixorient,	": vertex normals -> orient polygon");
	OPTSFC(opts,flipnormals,": flip orientations of normals and polygons");
	OPTSPC(opts,stretch,	"factor : stretch polylines");
	OPTSPC(opts,shownormals,"fsize : print normals as small segments");
	OPTSPC(opts,offset,	"factor : move vertices along their normals");
	OPTSPC(opts,noise,	"sdv : add Gaussian noise to vertices");
	OPTSFC(opts,intersect,	": intersect polygons to produce lines");
	OPTSFC(opts,twosided,	": make polygons two-sided");
	OPTSFC(opts,triangulate,": triangulate all faces with >3 vertices");
	OPTSPC(opts,tessellate,	"n : subdivide each triangle into n*n faces");
	opts.c("",":");
	OPTSFC(opts,statistics,	": print statistics");
	OPTSFC(opts,box,	": show bounding box");
	OPTSFC(opts,boxframe,	": output frame that will box data");
	OPTSFC(opts,nooutput,	": turn off a3d output");
	opts.c("",":");
	OPTSPC(opts,split,	"i : output frame every ith element");
	OPTSPC(opts,speedup,	"factor : increase 'split' every frame");
	OPTSPC(opts,frdelay,	"fsec : pause after each frame");
	OPTSPC(opts,eldelay,	"fsec : pause after each element");
	OPTSFC(opts,ctoascii,	": make output be ascii text");
	OPTSFC(opts,ctobinary,	": make output be binary");
	OPTSFC(opts,ctoold,	": convert to old a3d format");
	
	const char* filename="-";
	if (argc>1 && argv[1][0]!='-') {
		filename=argv[1];
		argc--; argv++;
	}
	RFile is(filename);
	RSA3dStream ia3d(is());
	if (!opts.parse(argc,argv) || !opts.noargs(argc,argv)) {
		opts.problem(argc,argv);
		return 1;
	}
	if (restrictf) {
		int obn; double zoom; int bin;
		isrestrictf=1;
		string str(restrictf);
		istringstream istr(str);
		assertx(FrameIO::read(istr,crestrictf,obn,zoom,bin));
	}
	if ((cusphr=cullsphere[3]) != 0)
		cusphc=Point(cullsphere[0],cullsphere[1],cullsphere[2]);
	if (transf) {
		int obn; double zoom; int bin;
		istransf=1;
		string str(transf);
		istringstream istr(str);
		assertx(FrameIO::read(istr,ctransf,obn,zoom,bin));
		if (!(isctransfinv=invert(ctransf,ctransfinv)))
			SHOWDF("Warning: uninvertible frame, normals lost\n");
	}
	if (joinlines) {
		join.hp=new HashPoint;
		join.graph=new Graph<int>;
	}
	for (int c=0;c<3;c++)
		cdiff.c[c]=diff[c],cspec.c[c]=spec[c],cphong.c[c]=phong[c];
	if (statistics || box || boxframe) nooutput=1;
	bbox.clear();
	inter.bb.clear();
	fsplit=split;
	assertx(!(ctoascii && ctobinary));
	if (ctoascii) setenv("A3DBINARY","0",1);
	if (ctobinary) setenv("A3DBINARY","",1);
	if (ctoold) setenv("A3DOLD","",1);
	doit(ia3d);
	CleanUp();
	return 0;
}

static void doit(RSA3dStream& ia3d)
{
	A3dElem el;
	for (;;) {
		ia3d.read(el);
		if (loop(el)) break;
	}
	if (joinlines) dojoinlines();
	if (triangulate)
		SHOWDF("triangulation: %d polyg (%d triang) -> %d triang\n",
		       tri.npolyb,tri.ntrib,tri.ntria);
	if (tessellate)
		SHOWDF("tessellate: %d triang -> %d triang\n",
		      tess.ntrib,tess.ntria);
	if (intersect) {
		dointersect();
		SHOWDF("intersect: added %d edges\n",inter.nedges);
	}
	if (statistics) {
		SHOWDF("Polygons\n");
		SHOWDF(" %s",Spnvert.namestring());
		SHOWDF(" %s",Spedgel.namestring());
		SHOWDF(" %s",Sparea.namestring());
		SHOWDF(" tot_area: %g\n",Sparea.sum());
		SHOWDF(" %s",Splanar.namestring());
		if (Splanar.num()<Spnvert.num())
			SHOWDF("  some zero area polygons not counted!\n");
		SHOWDF("Polylines\n");
		SHOWDF(" %s",Slnvert.namestring());
		SHOWDF(" %s",Sledgel.namestring());
		SHOWDF("Points:\n");
		SHOWDF(" %s",Sptnor.namestring());
	}
	if (box) {
		SHOWF("%g %g %g\n",bbox[0][0],bbox[0][1],bbox[0][2]);
		SHOWF("%g %g %g\n",bbox[1][0],bbox[1][1],bbox[1][2]);
	}
	if (boxframe) FrameIO::write(cout,bbox.getFrameToCube(),0,0,0);
	if (ncullsphere) SHOWDF("ncullsphere=%d\n",ncullsphere);
	if (nmindis) SHOWDF("nmindis=%d\n",nmindis);
	if (nfixorient) SHOWDF("nfixorient=%d\n",nfixorient);
	if (ndegen) SHOWDF("ndegen=%d\n",ndegen);
}

static int loop(A3dElem& el)
{
        int i;
	if (el.type()==A3dElem::TEndFile) return 1;
	int polyg=el.type()==A3dElem::TPolygon;
	int polyl=el.type()==A3dElem::TPolyline;
	int point=el.type()==A3dElem::TPoint;
	if (!polyg && !polyl && !point) {
		if (!onlypoly) doout(el);
		if (el.type()==A3dElem::TEndFrame) delayframe();
		return 0;
	}
	if (joinlines && polyl) {
//		for (int i=0;i<el.num();i++) {
		for (i=0;i<el.num();i++) {
			int vi=join.hp->enter(el[i].p);
			assertx(vi<=join.pa.num());
			if (vi==join.pa.num()) {
				join.pa+=el[i].p;
				join.graph->enter(vi);
			}
		}
		for (i=0;i<el.num()-1;i++) {
			int v0=join.hp->enter(el[i].p);
			int v1=join.hp->enter(el[i+1].p);
			join.graph->enteru(v0,v1);
		}
		return 0;
	}
	if (polyg && nopolygons) return 0;
	if (polyl && nopolylines) return 0;
	if (point && nopoints) return 0;
	if (polyg && fixdegen && isdegen(el)) { ndegen++; return 0; }
	if (isrestrictf && outofbounds(el)) return 0;
	static int nevery=0;
	nevery++;
	if (every>1 && nevery%every!=1) return 0;
	if (first && nevery>first) return 1;
	if (cusphr && point && dist2(el[0].p,cusphc)<=square(cusphr)) {
		ncullsphere++;
		return 0;
	}
	if (mindis && point && domindis(el[0].p)) {
		nmindis++;
		return 0;
	}
	Vector pnor(3,0,0);
	if (optnormals && polyg) pnor=el.pnormal();
//	for (int i=0;i<el.num();i++) {
	for (i=0;i<el.num();i++) {
		if (nonormals || optnormals && !compare(el[i].n,pnor,1e-6))
			el[i].n=Vector(0,0,0);
		if (nocolor) el[i].c=zerovertexcolor;
		int validcol=0;
		if (cdiff.c[0]>=0) el[i].c.d=cdiff,validcol=1;
		if (cspec.c[0]>=0) el[i].c.s=cspec,validcol=1;
		if (validcol && !el[i].c.g.c[0]) el[i].c.g.c[0]=1;
		if (cphong.c[0]>=0) el[i].c.g=cphong;
		if (istransf) {
			el[i].p*=ctransf;
			if (!isctransfinv) {
				el[i].n=Vector(0,0,0);
			} else if (!el[i].n.iszero()) {
				el[i].n=ctransfinv*el[i].n;
				assertw1(el[i].n.normalize());
			}
		}
		if (gnormalize && !el[i].n.iszero())
			assertx(el[i].n.normalize());
	}
	if (fixorient && polyg && polygonneedsflip(el))
		flippolygon(el),nfixorient++;
	if (flipnormals) flippolygon(el);
	if (stretch && polyl && el.num()==2) {
		if (stretch>0) {
			el[1].p+=(el[1].p-el[0].p)*stretch;
		} else {
			el[0].p+=(el[0].p-el[1].p)*stretch;
		}
	}
	if (shownormals) doshownormals(el);
	for (i=0;i<el.num();i++) {
		if (offset) {
			Vector nor=el[i].n;
			if (nor.iszero() && polyg) {
				if (pnor[0]>2) pnor=el.pnormal();
				nor=pnor;
			}
			if (!nor.iszero()) el[i].p+=nor*offset;
		}
		if (noise)
			for (int c=0;c<3;c++)
				el[i].p[c]+=Random::G.gauss()*noise;
	}
	if (intersect && polyg) {
		Polygon* np=new Polygon(el);
		Bbox bb;
		np->getbbox(bb);
		inter.bb.takeunion(bb);
		inter.Spolyg.push(np);
	}
	if (twosided && polyg) { pass1(el); flippolygon(el); }
	pass1(el);
	return 0;
}

static A3dVertex A3dVertexAffComb(const A3dElem& el, double* w=0)
{
	A3dVertex vavg(Point(0,0,0),Vector(0,0,0),zerovertexcolor);
	Vector pnor=el.pnormal();
	for (int i=0;i<el.num();i++) {
		double a=w?w[i]:1/double(el.num());
		vavg.n+=(el[i].n.iszero()?pnor:el[i].n)*a;
		for (int c=0;c<3;c++) {
			vavg.p[c]+=el[i].p[c]*a;
			vavg.c.d.c[c]+=el[i].c.d.c[c]*a;
			vavg.c.s.c[c]+=el[i].c.s.c[c]*a;
			vavg.c.g.c[c]+=el[i].c.g.c[c]*a;
		}
	}
	return vavg;
}

// maybe do triangulate
static void pass1(const A3dElem& el)
{
        int i;
	if (!triangulate || el.type()!=A3dElem::TPolygon) { pass2(el);return; }
	tri.npolyb++;
	if (el.num()==3) { tri.ntrib++,tri.ntria++; pass2(el); return; }
	static A3dElem el2;
	el2.init(el.type(),el.binary(),3);
	if (el.num()==4) {
//		for (int i=1;i<el.num();i++)
		for (i=1;i<el.num();i++)
			if (compare(el[i].n,el[0].n,1e-6) ||
			    compare(el[i].c.d,el[0].c.d) ||
			    compare(el[i].c.s,el[0].c.s) ||
			    compare(el[i].c.g,el[0].c.g)) break;
		if (i==el.num()) { // they all match
			for (i=0;i<3;i++) el2[i]=el[i];
			pass2(el2);
			for (i=0;i<3;i++) el2[i]=el[(i+2)%4];
			pass2(el2);
			tri.ntria+=2;
			return;
		}
	}
	A3dVertex vavg=A3dVertexAffComb(el);
//	for (int i=0;i<el.num();i++) {
	for (i=0;i<el.num();i++) {
		el2[0]=el[i];
		el2[1]=el[(i+1)%el.num()];
		el2[2]=vavg;
		pass2(el2);
	}
	tri.ntria+=el.num();
}

static A3dVertex getvert(const A3dElem& el, int i, int j, int k, int nt)
{
	assertx(i>=0 && j>=0 && k>=0 && i<=nt && j<=nt && k<=nt);
	assertx(i+j+k==nt && el.num()==3);
	double w[3];
	w[0]=double(i)/nt;
	w[1]=double(j)/nt;
	w[2]=double(k)/nt;
	return A3dVertexAffComb(el,w);
}

// maybe do tessellate
static void pass2(const A3dElem& el)
{
	if (!tessellate || el.type()!=A3dElem::TPolygon) { pass3(el); return; }
	if (assertw1(el.num()==3)) { pass3(el); return; }
	int nt=tessellate;
	tess.ntrib++;
	static A3dElem el2;
	el2.init(el.type(),el.binary(),3);
	for (int i=0;i<nt;i++) {
		for (int j=0;j<nt-i;j++) {
			A3dVertex v0=getvert(el,nt-i-j,j,i,nt);
			A3dVertex v1=getvert(el,nt-i-j-1,j+1,i,nt);
			A3dVertex vn=getvert(el,nt-i-j-1,j,i+1,nt);
			el2[0]=v0; el2[1]=v1; el2[2]=vn;
			pass3(el2);
			if (i) {
				A3dVertex vp=getvert(el,nt-i-j,j+1,i-1,nt);
				el2[0]=v1; el2[1]=v0; el2[2]=vp;
				pass3(el2);
			}
		}
	}
	tess.ntria+=nt*nt;
}

// split and output statistics
static void pass3(const A3dElem& el)
{
	static int nelem=0;
	if (split && nelem++>=split) {
		oa3d.writeEndFrame(el.binary());
		delayframe();
		oa3d.writeEndObject(el.binary(),1,0);
		if (speedup) fsplit*=speedup;
		split=int(fsplit+.01);
		nelem=1;
	}
	if (statistics) dostat(el);
	if (box || boxframe)
		for (int i=0;i<el.num();i++) bbox.takeunion(el[i].p);
	doout(el);
	delayelement();
}

static void doout(const A3dElem& el)
{
	if (nooutput) return;
	oa3d.write(el);
}

static void delayframe()
{
	if (!frdelay) return;
	assertx(frdelay>0);
	Delay(frdelay);
}

static void delayelement()
{
	if (!eldelay) return;
	assertx(eldelay>0);
	Delay(eldelay);
}

static int isdegen(const A3dElem& el)
{
	if (assertw1(el.num()>=3)) return 1;
	Vector vt(0,0,0);
	for (int i=1;i<el.num()-1;i++)
		vt+=cross(el[0].p,el[i].p,el[i+1].p);
	double area=.5*mag(vt);
	return !area;
}

static int outofbounds(const A3dElem& el)
{
	for (int i=0;i<el.num();i++) {
		Point p=el[i].p*crestrictf;
		for (int c=0;c<3;c++) if (p[c]<0 || p[c]>1) return 1;
	}
	return 0;
}

static int polygonneedsflip(const A3dElem& el)
{
	Vector va(0,0,0);
	Vector pnor=el.pnormal();
	for (int i=0;i<el.num();i++)
		va+=el[i].n.iszero()?pnor:el[i].n;
	return dot(va,pnor)<0;
}

static void flippolygon(A3dElem& el)
{
        int i;
	A3dVertex v;
//	for (int i=1;i<=(el.num()-1)/2;i++) {
	for (i=1;i<=(el.num()-1)/2;i++) {
		v=el[i];
		el[i]=el[el.num()-i];
		el[el.num()-i]=v;
	}
	for (i=0;i<el.num();i++)
		if (!el[i].n.iszero()) el[i].n=-el[i].n;
}

static void doshownormals(const A3dElem& el)
{
        int i;
	double c=0;
	if (el.num()<2) {
		c=1;
	} else {
		Array<Point> pa;
//		for (int i=0;i<el.num();i++) pa+=el[i].p;
		for (i=0;i<el.num();i++) pa+=el[i].p;
		Point pavg=centroid(pa);
		for (i=0;i<el.num();i++) c+=dist(pa[i],pavg);
		c/=el.num();
	}
	static A3dElem el2;
	el2.init(A3dElem::TPolyline,el.binary(),2);
	Vector pnor(3,0,0);
//	for (int i=0;i<el.num();i++) {
	for (i=0;i<el.num();i++) {
		Vector nor=el[i].n;
		if (nor.iszero()) {
			if (el.type()!=A3dElem::TPolygon) continue;
			if (pnor[0]>2) pnor=el.pnormal();
			nor=pnor;
		}
		el2[0]=el2[1]=el[i];
		el2[0].n=el2[1].n=nor;
		el2[1].p+=nor*c*.2*shownormals;
		pass1(el2);
	}
}

static void dostat(const A3dElem& el)
{
        int i;
	if (el.type()==A3dElem::TPolyline) {
		Slnvert+=el.num();
		for (int i=0;i<el.num()-1;i++)
			Sledgel+=dist(el[i].p,el[i+1].p);
		return;
	}
	if (el.type()==A3dElem::TPoint) {
		Sptnor+=!el[0].n.iszero();
		return;
	}
	Spnvert+=el.num();
//	for (int i=0;i<el.num();i++)
	for (i=0;i<el.num();i++)
		Spedgel+=dist(el[i].p,el[(i+1)%el.num()].p);
	assertx(el.num()>=3);
	Vector vt(0,0,0);
	for (i=1;i<el.num()-1;i++)
		vt+=cross(el[0].p,el[i].p,el[i+1].p);
	double area=.5*mag(vt);
	Sparea+=area;
	if (area) {
		assertw1(vt.normalize());
		double sumd=0;
		for (i=0;i<el.num();i++) {
			const Point& p=el[i].p;
			sumd+=pvdot(p,vt);
		}
		double d=sumd/el.num();
		double tol=0;
		for (i=0;i<el.num();i++) {
			const Point& p=el[i].p;
			tol=max(tol,pvdot(p,vt)-d);
		}
		Splanar+=tol/mysqrt(area);
	}
}

static int considerpoly(Polygon* id, double [][2], KdtreeNode *)
{
	const Polygon& p1=*id;
	const Polygon& p2=*inter.polyg;
	static Array<Point> pa;
	IntersectPolyPoly(p1,p2,pa);
	if (!pa.num()) return 0;
	static A3dElem el2;
	el2.init(A3dElem::TPolyline,0,2);
	for (int i=0;i<pa.num();i+=2) {
		el2[0]=A3dVertex(pa[i],Vector(0,0,0),A3dVertexColor::White);
		el2[1]=A3dVertex(pa[i+1],Vector(0,0,0),A3dVertexColor::White);
		doout(el2);
		inter.nedges++;
	}
	return 0;
}

static void dointersect()
{
	if (inter.Spolyg.empty()) return;
	Frame f=inter.bb.getFrameToCube();
	Kdtree<Polygon*> kd(3,8);
	ForStack(inter.Spolyg,Polygon*,poly) {
		inter.polyg=poly;
		Bbox bb;
		inter.polyg->getbbox(bb);
		bb[0]*=f;
		bb[1]*=f;
		double bv[3][2];
		bv[0][0]=bb[0][0]; bv[1][0]=bb[0][1]; bv[2][0]=bb[0][2];
		bv[0][1]=bb[1][0]; bv[1][1]=bb[1][1]; bv[2][1]=bb[1][2];
		(void) kd.search(bv, considerpoly, 0);
		// cast for DECCXX
		kd.enter(inter.polyg, bv);
	} EndFor;
	while (!inter.Spolyg.empty()) delete inter.Spolyg.pop();
}

static int domindis(const Point& p)
{
	static int pn=1;
	static PointSpatial* SPp;
	if (!SPp) SPp=new PointSpatial(30); // no delete
	SpatialSearch ss(*SPp,p,mindis);    // look no farther than mindis
	if (!ss.done()) {
		double dis2;
		(void)Conv<int>::d(ss.next(&dis2));
		if (dis2<square(mindis)) return 1;
	}
	SPp->enter(Conv<int>::e(pn++),new Point(p)); // no delete
	return 0;
}

static void dojoinlines()
{
	joinlines=0;		// so that loop() does not capture new lines
	A3dElem el;
	Graph<int>& graph=*join.graph;
	Set<int> odd;		// vertices of odd degree
	For (GraphVertices<int> gv(graph);gv;gv.next()) {
		assertx(graph.outdegree(gv())>0);
		if (graph.outdegree(gv())%2) odd.enter(gv());
	} EndFor;
	for (;;) {
		GraphVertices<int> gv(graph);
		if (!gv) break;
		int vi=gv();
		if (!odd.empty()) {
			vi=odd.removeone();
		} else {
			odd.enter(vi);
		}
		el.init(A3dElem::TPolyline);
		for (;;) {
			el+=A3dVertex(join.pa[vi],Vector(0,0,0),
				      A3dVertexColor::White);
			GraphEdges<int> ge(graph,vi);
			if (!ge) break;
			int vn=ge();
			assertx(graph.removeu(vi,vn));
			if (!graph.outdegree(vi)) assertx(graph.remove(vi));
			vi=vn;
		}
		assertx(graph.remove(vi));
		assertx(odd.remove(vi));
		loop(el);
	}
	assertx(odd.empty());
	delete join.hp;
	delete join.graph;
}
