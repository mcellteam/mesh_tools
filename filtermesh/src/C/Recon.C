// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

// Recon2
// From a set of points, build an implicit function that represents the signed
// distance to the manifold.  Then extract its zero set.
// 03/16/93 ported from C version

// Notes:
//  The pcnor's should point outwards (I define distance function to be
//   negative inside the object).
//  The growing of pc's is not ideal.
//  The goals are:
//   - to only consider points on the current sheet
//   - to get a good tangent plane pcnor[], so need sufficient number of
//     points to form good estimator (e.g. 4 points near sharp corner can be
//     terrible) in the presence of noise, should consider a larger set
//   - to be adaptive to sampling density and noise

// 04/27/93: idea: use samplingd for everything, don't use npertp
//  to estimate tp, gather all points within samplingd
//  to connect tp, connect tp's whose original co's are within samplingd
//  signed distance is undefined if projp >samplingd from nearest point
//  no longer need EMST, assume this new graph represents components
// 09/09/93: all coordinates transformed internally into unit box.

#include "Hh.h"
#include "Options.h"
#include "A3dStream.h"
#include "Homogeneous.h"
#include "Mk3d.h"
#include "Mklib.h"
#include "Stat.h"
#include "Polygon.h"
#include "Spatial.h"
#include "Principal.h"
#include "Graph.h"
#include "GraphOp.h"
#include "Contour.h"
#include "GMesh.h"
#include "MeshOp.h"
#include "Timer.h"
#include "Array.h"
#include "Bbox.h"
#include "FrameIO.h"
#include "Set.h"
#include "Queue.h"
#include "Stack.h"

#include <iostream>

using std::cin;
using std::cout;
using std::ostream;
using std::ofstream;

// necessary
static const char* rootname=0;
static const char* what="m";
static double samplingd=0;
static int gridsize=0;

// optional
static double unsigneddis=0;
static int prop=2;
static int usenormals=0;
static int maxkintp=20;
static int minkintp=4;

static int num;			// # data points
static int is3D;		// is it a 3D problem (vs. 2D)
static int minora;		// minor axis: 1 in 2D, 2 in 3D
static int havenormals;		// data contains normal information

static char* header;		// header to place at top of output files
static Frame xform;		// original pts -> pts in unit cube
static Frame xformi;		// inverse
static PointSpatial* SPp;	// spatial partition on co
static PointSpatial* SPpc;	// spatial partition on pcorg
static STAT(Scorr);		// correlation in orientation propagation
static Graph<int>* gpcpseudo;	// Riemannian on pc centers (based on co)
static Graph<int>* gpcpath;	// path of orientation propagation
static GMesh mesh;

static Mk3d* iod;		// lines: co[i] to pcorg[i]
static Mk3d* iou;		// unoriented pc boxes (volumes)
static Mk3d* ioU;		// unoriented tangent planes (polygons)
static Mk3d* iof;		// flat tangent subspaces (polyline)
static Mk3d* iog;		// pseudograph over pcorg
static Mk3d* iop;		// orientation propagation path
static Mk3d* ioO;		// oriented tangent planes
static Mk3d* ioG;		// ioG with normals
static Mk3d* iol;		// projections co[i] to nearest tp
static Mk3d* ioC;		// visited marching cubes
static Mk3d* iom;		// final mesh (not a3d!), a3d curve if !is3D

static Array<Point> co;		// data points
static Array<Vector> nor;	// normal at point (!nor or iszero() if none)
// The following are allocated after points have been read in
static Array<Point> pcorg;	// origins of tangent planes
static Array<Vector> pcnor;	// normals of tangent planes
static Array<int> pciso;	// tangent plane has been oriented
static Array<Frame> pctrans;	// defined if ioO

// Declarations
static void doit();
static void doread();
static void computexform();
static void initoutput();
static void doarg(char ch, Mk3d*& mk);
static void closemk(Mk3d*& mk);
static void doprincipal();
static void computetp(int i, int& n, Frame& f);
static void printprincipal(const Frame& f);
static void drawpcextent(Mk3d* mk);
static void drawpclinear(Mk3d* mk);
static void doorient();
static void orientset(const Set<int>& nodes);
static void addexteriororientation(const Set<int>& nodes);
static void removeexteriororientation();
static double pc_corr(Univ e1, Univ e2);
static double pc_dot(int i, int j);
static void propagatealongpath(int i);
static void showprop(int i, int j, double dotp);
static void draworientedtps();
static void docontour();
static double evalpoint(const Point& p);
static void printdirectedseg(Mk3d* mk, const Point& p1, const Point& p2,
			     const A3dColor& col);
static double computesigned(const Point& p, Point& proj);
static void outputborder(const Polygon& poly);
static void outputcontour(const Polygon& poly);
static void prtgraph(Mk3d* mk, const Graph<int>* g, const Point* pa,
		     const Vector* pn);
static double computeunsigned(const Point& p, Point& proj);

int main(int argc, char const** argv)
{
	header=newString(hform("# Created using\n%s",CreateHeader(argc,argv)));
	Options opts;
	OPTSPC(opts,rootname,	"string : name for output files (optional)");
	OPTSPC(opts,what,	"string : output codes 'duUfgpoGlCm' ('m')");
	OPTSPC(opts,samplingd,	"f : sampling density + noise (delta+rho)");
	OPTSPC(opts,gridsize,	"n : contouring # grid cells (opt.)");
	OPTSPC(opts,maxkintp,	"k : max # points in tp (def. 20)");
	OPTSPC(opts,minkintp,	"k : min # points in tp (def. 4)");
	OPTSPC(opts,unsigneddis,"f : use unsigned distance, set value");
	OPTSPC(opts,prop,	"i : orient. prop. (0=naive,1=emst,2=mst)");
	OPTSPC(opts,usenormals,	"b : use normals given with data");
	if (!opts.parse(argc,argv) || !opts.noargs(argc,argv)) {
		opts.problem(argc,argv);
		return 1;
	}
	assertx(samplingd);
	doit();
	CleanUp();
	// iom closed here so mesh comes after everything else in file
	if (iom && is3D) {
		ForMeshVertex(mesh,v) {
			mesh.setPoint(v,mesh.point(v)*xformi);
		} EndFor;
		mesh.write(reinterpret_cast<WSA3dStream &>(iom->oa3d()).os());
	}
	closemk(iom);
	delete[] header;
	mesh.clear();
	return 0;
}

static void doit()
{
	TIMER(Recon);
	doread();
	computexform();
	if (!gridsize) gridsize=int(1./samplingd+.5);
	SHOWDF("gridsize=%d\n",gridsize);
	assertx(gridsize>1);
	initoutput();
	NEST {
		pcorg.init(num);
		pcnor.init(num);
		pciso.init(num);
		if (ioO) pctrans.init(num);
		for (int i=0;i<num;i++) pciso[i]=0;
	}
	NEST {
		TIMER(_SPp);
		int n=is3D?(num>5000?36:20):(num>1000?36:20);
		SPp=new PointSpatial(n);
		for (int i=0;i<num;i++)
			SPp->enter(Conv<int>::e(i),&co[i]);
	}
	if (!unsigneddis) {
		gpcpseudo=new Graph<int>;
		for (int i=0;i<num;i++) gpcpseudo->enter(i);
		doprincipal();
		NEST {
			TIMER(_SPpc);
			int n=is3D?(num>5000?36:20):(num>1000?36:20);
			SPpc=new PointSpatial(n);
			for (int i=0;i<num;i++)
				SPpc->enter(Conv<int>::e(i),&pcorg[i]);
		}
		doorient();
		delete gpcpseudo;
	}
	if (iol || ioC || iom) docontour();
	if (iom && is3D) SHOWDF(MeshGenusString(mesh));
	delete SPp;
	delete SPpc;
}

static void doread()
{
	TIMER(_read);
	RSA3dStream ia3d(cin);
	static A3dElem el;
	int nnor=0;
	for (;;) {
		ia3d.read(el);
		if (el.type()==A3dElem::TEndFile) break;
		if (el.type()==A3dElem::TComment) continue;
		assertx(el.type()==A3dElem::TPoint);
		co+=el[0].p;
		nor+=el[0].n;
		if (co[num][0]) is3D=1;
		if (!nor[num].iszero()) {
			nnor++;
			if (usenormals) assertw(nor[num].normalize());
		}
		num++;
	}
	SHOWDF("%d points (with %d normals), %dD analysis\n",num,nnor,2+is3D);
	assertx(num>1);
	minora=is3D?2:1;
	if (nnor>0 && !usenormals) SHOWDF("ignoring normals!\n");
	havenormals=usenormals && nnor>0;
	if (havenormals) assertx(prop==2);
	co.squeeze();
	if (!havenormals) nor.init(0);
	nor.squeeze();
}

static void computexform()
{
        int i;
	Bbox bb;
	bb.clear();
//	for (int i=0;i<num;i++) bb.takeunion(co[i]);
	for (i=0;i<num;i++) bb.takeunion(co[i]);
	xform=bb.getFrameToSmallCube();
	if (!is3D) xform.p[0]=0; // preserve x==0
	double xformscale=xform[0][0];
	SHOWDF("Applying xform: %s",FrameIO::string(xform,1,0));
	xformi=~xform;
	for (i=0;i<num;i++) co[i]*=xform;
	// nor[] is unchanged
	samplingd*=xformscale;
	unsigneddis*=xformscale;
}

static void initoutput()
{
	doarg('d',iod);
	doarg('u',iou);
	doarg('U',ioU);
	doarg('f',iof);
	doarg('g',iog);
	doarg('p',iop);
	doarg('O',ioO);
	doarg('G',ioG);
	doarg('l',iol);
	doarg('C',ioC);
	doarg('m',iom);
}

static void doarg(char ch, Mk3d*& mk)
{
	if (!strchr(what,ch)) return;
	ostream* os;
	if (rootname) {
		os=new ofstream(hform("%s.%c",rootname,ch));
	} else {
		assertx(strlen(what)==1);
		os=&cout;	// instead, could have used: ofstream(1)
	}
	*os << header;
	WSA3dStream* oa3d=new WSA3dStream(*os);
	mk=new Mk3d(*oa3d);
	mk->oa3d().writeComment(hform(" Output of option '%c'",ch));
	mk->oa3d().flush();
}

static void closemk(Mk3d*& mk)
{
	if (!mk) return;
	WSA3dStream* oa3d = reinterpret_cast<WSA3dStream *>(&mk->oa3d()); // sort of a hack
	ostream* os=&oa3d->os();
	delete mk; mk=0;
	delete oa3d;
	if (os!=&cout) delete os;
}

static void doprincipal()
{
	TIMER(_principal);
	// Statistics in reverse order of printout
	STAT(Sr21); STAT(Sr20); STAT(Sr10);
	STAT(Slen2); STAT(Slen1); STAT(Slen0);
	STAT(Snei);
	for (int i=0;i<num;i++) {
		int n;
		Frame f;
		computetp(i,n,f);
		if (ioO) pctrans[i]=f;
		Snei+=n;
		double len0=mag(f.v[0]), len1=mag(f.v[1]), len2=mag(f.v[2]);
		assertx(len2>0); // PrincipalComponents() should do this
		Slen0+=len0; Slen1+=len1; Slen2+=len2;
		Sr10+=len1/len0; Sr20+=len2/len0; Sr21+=len2/len1;
		pcorg[i]=f.p;
		pcnor[i]=f.v[minora];
		assertx(pcnor[i].normalize());
		printprincipal(f);
	}
	closemk(iou); closemk(iof); closemk(ioU);
	if (iod) {
		iod->diffuse(1,1,.3); iod->specular(0,0,0); iod->phong(1);
		for (int i=0;i<num;i++) {
			iod->point(co[i]*xformi);
			iod->point(pcorg[i]*xformi);
			iod->endpolyline();
		}
	}
	closemk(iod);
}

static void computetp(int i, int& n, Frame& f)
{
	static Array<Point> pa;
	pa.init(0);
	SpatialSearch ss(*SPp,co[i]);
	double eimag[3];
	for (;;) {
		assertx(!ss.done());
		double dis2; int pi=Conv<int>::d(ss.next(&dis2));
		if (pa.num()>=minkintp && dis2>square(samplingd) ||
		    pa.num()>=maxkintp) break;
		pa+=co[pi];
		if (pi!=i) gpcpseudo->enteru(i,pi);
	}
	PrincipalComponents(pa,pa.num(),f,eimag);
	n=pa.num();
}

static void printprincipal(const Frame& f)
{
	if (iou) {
		iou->diffuse(.8,.5,.4); iou->specular(.5,.5,.5); iou->phong(3);
		iou->push();
		iou->apply(f*xformi); drawpcextent(iou);
		iou->pop();
	}
	if (iof) {
		iof->diffuse(1,1,1); iof->specular(0,0,0); iof->phong(1);
		iof->beginForcePolyline(1);
		iof->push();
		iof->apply(f*xformi); drawpclinear(iof);
		iof->pop();
		iof->endForcePolyline();
	}
	if (ioU) {
		ioU->diffuse(.7,.3,.1); ioU->specular(.6,.6,.2); ioU->phong(3);
		ioU->push();
		ioU->apply(f*xformi); drawpclinear(ioU);
		ioU->pop();
	}
}

static void drawpcextent(Mk3d* mk)
{
	mk->push();
	mk->scale(2);
	Mklib mklib(*mk);
	// frame is defined in terms of principal frame!
	if (is3D) {
		mklib.cubeO();
	} else {
		mklib.squareU();
	}
	mk->pop();
}

static void drawpclinear(Mk3d* mk)
{
	if (is3D) {
		Mklib mklib(*mk);
		mk->push(); mk->scale(2); mklib.squareU(); mk->pop();
	} else {
		mk->point(+1,0,0);
		mk->point(-1,0,0);
		mk->endpolyline();
	}
}

static double pcdist(Univ e1, Univ e2)
{
	return dist(pcorg[Conv<int>::d(e1)],pcorg[Conv<int>::d(e2)]);
}

static void doorient()
{
	TIMER(_orient);
	NEST {
		Stat* st=newGraphEdgeStats(*gpcpseudo,pcdist);
		st->setname("gpcpseudo"); st->setprint(1);
		delete st;
	}
	NEST {
		TIMER(__graphnumcompon);
		int nc=GraphNumComponents(*gpcpseudo);
		SHOWDF("Number of components: %d\n",nc);
		if (nc>1) SHOWDF("*** #comp>1, may want larger -samp\n");
	}
	if (iog) {
		iog->diffuse(.7,.7,.7); iog->specular(0,0,0); iog->phong(1);
		prtgraph(iog,gpcpseudo,pcorg,0);
		closemk(iog);
	}
	// Now treat each connected component of gpcpseudo separately.
	Set<int> setnotvis; ForIndex(i,num) { setnotvis.enter(i); } EndFor;
	for (;!setnotvis.empty();) {
		Set<int> nodes;
		Queue<int> queue;
		int fi=setnotvis.getone();
		nodes.enter(fi); queue.enqueue(fi);
		for (;!queue.empty();) {
			int i=queue.dequeue();
			assertx(setnotvis.remove(i));
			For (GraphEdges<int> ge(*gpcpseudo,i);ge;ge.next()) {
				if (nodes.add(ge())) queue.enqueue(ge());
			} EndFor;
		}
		orientset(nodes);
	}
	Scorr.terminate();
	for (int i=0;i<num;i++) assertx(pciso[i]);
	closemk(iop);
	if (ioO) draworientedtps();
	closemk(ioO);
	if (ioG) {
		ioG->diffuse(.7,.7,.7); ioG->specular(0,0,0); ioG->phong(1);
		prtgraph(ioG,gpcpseudo,pcorg,pcnor);
		closemk(ioG);
	}
}

static void orientset(const Set<int>& nodes)
{
	SHOWDF("component with %d points\n",nodes.num());
	addexteriororientation(nodes);
	gpcpath=new Graph<int>;
	ForSet(nodes,int,i) { gpcpath->enter(i); } EndFor;
	gpcpath->enter(num);
	NEST {
		TIMER(__graphmst);
		// must be connected here!
		assertx(GraphMst(*gpcpseudo,pc_corr,*gpcpath));
	}
	int nextlink=gpcpath->outdegree(num);
	if (nextlink>1) SHOWDF(" num_exteriorlinks_used=%d\n",nextlink);
	propagatealongpath(num);
	delete gpcpath,gpcpath=0;
	removeexteriororientation();
}

static void addexteriororientation(const Set<int>& nodes)
{
	// vertex num is a pseudo-node used for outside orientation
	gpcpseudo->enter(num);
	if (havenormals) {
		// add pseudo-edges from "exterior" to points with normals
		ForSet(nodes,int,i) {
			if (nor[i].iszero()) continue;
			gpcpseudo->enteru(i,num);
		} EndFor;
	} else {
		// add 1 pseudo-edge to point with largest z value
		double maxz=-1e30; int maxi=-1;
		ForSet(nodes,int,i) {
			if (pcorg[i][2]>maxz) maxz=pcorg[i][2],maxi=i;
		} EndFor;
		gpcpseudo->enteru(maxi,num);
	}
}

static void removeexteriororientation()
{
	Stack<int> stack;
	For (GraphEdges<int> ge(*gpcpseudo,num);ge;ge.next()) {
		stack.push(ge());
	} EndFor;
	ForStack(stack,int,i) { gpcpseudo->removeu(num,i); } EndFor;
	gpcpseudo->remove(num);
}

static double pc_corr(Univ e1, Univ e2)
{
	int i=Conv<int>::d(e1), j=Conv<int>::d(e2);
	if (j==num && i<num) return pc_corr(e2,e1);
	assertx(i>=0 && j>=0 && i<=num && j<num);
	if (havenormals) assertx(prop==2);
	double vdot, corr;
	switch (prop) {
	  case 0:
		corr=1;		// any path is ok
	  bcase 1:
		if (i==num) {
			corr=0; // single exterior link
		} else {
			corr=dist2(pcorg[i],pcorg[j]);
		}
	  bcase 2:
		if (i==num) {
			if (havenormals) {
				assertx(!nor[j].iszero());
				vdot=dot(nor[j],pcnor[j]);
			} else {
				vdot=1; // single exterior link
			}
		} else {
			vdot=dot(pcnor[i],pcnor[j]);
		}
		corr=2-fabs(vdot);
	  bdefault:
		assertnever(""); corr=0;
	}
	return corr;
}

static double pc_dot(int i, int j)
{
	assertx(i>=0 && j>=0 && i<=num && j<num);
	if (i==num) {
		if (havenormals) {
			return dot(nor[j],pcnor[j]);
		} else {
			return pcnor[j][2]<0?-1:1;
		}
	} else {
		return dot(pcnor[i],pcnor[j]);
	}
}

// Propagate orientation along tree gpcpath from vertex i (orig. num) using
// recursive DFS.
static void propagatealongpath(int i)
{
	assertx(i>=0 && i<=num);
	if (i<num) assertx(pciso[i]);
	for (GraphEdges<int> ge(*gpcpath,i);ge;ge.next()) {
		int j=ge();
		assertx(j>=0 && j<=num);
		if (j==num || pciso[j]) continue; // immediate caller
		double corr=pc_dot(i,j);
		Scorr+=fabs(corr);
		if (corr<0) pcnor[j]=-pcnor[j];
		pciso[j]=1;
		showprop(i,j,fabs(corr));
		propagatealongpath(j);
	}
}

static void showprop(int i, int j, double dotp)
{
	assertx(i>=0 && i<=num && j>=0 && j<num && dotp>=0);
	if (i==num || !iop) return ;
	double f=min((1-dotp)*10.,1.); // low dotp, high f are signific.
	iop->diffuse(.2+.8*f,.8+.2*f,.5+.5*f);
	iop->specular(0,0,0); iop->phong(3);
	iop->point(pcorg[i]*xformi); iop->normal(pcnor[i]);
	iop->point(pcorg[j]*xformi); iop->normal(pcnor[j]);
	iop->endpolyline();
}

static void draworientedtps()
{
	ioO->diffuse(.7,.3,.1); ioO->specular(.6,.6,.2); ioO->phong(3);
	for (int i=0;i<num;i++) {
		Frame& f=pctrans[i];
		if (dot(f.v[minora],pcnor[i])<0) {
			// flip 2 of the axes to keep right hand rule
			f.v[minora]=-f.v[minora];
			f.v[0]=-f.v[0];
		}
		ioO->push();
		ioO->apply(f*xformi);
		drawpclinear(ioO);
		if (!is3D) {
			ioO->point(0,0,0);
			ioO->point(0,1,0);
			ioO->endpolyline();
		}
		ioO->pop();
	}
}

static void docontour()
{
        int i;
	TIMER(_contour);
	Contour contour(is3D,evalpoint,gridsize,&cout);
	if (iom &&  is3D) contour.setOutputMesh(mesh);
	if (iom && !is3D) contour.setOutputContour(outputcontour);
	if (ioC) contour.setOutputBorder(outputborder);
	if (!is3D) assertx(!pcorg[0][0]);
	if (unsigneddis) {
		Point p=co[0];
//		for (int i=1;i<num;i++) if (co[i][2]>p[2]) p=co[i];
		for (i=1;i<num;i++) if (co[i][2]>p[2]) p=co[i];
		for (i=0;i<20;i++) {
			p[2]+=double(i)/gridsize;
			if (contour.marchFrom(p)>1) break;
		}
	} else {
		for (int i=0;i<num;i++)
			(void)contour.marchFrom(pcorg[i]);
	}
	// iom closed back in main
	closemk(ioC);
	closemk(iol);
}

static double evalpoint(const Point& p)
{
	if (!is3D) assertx(!p[0]);
	Point proj;
	double dis=unsigneddis?computeunsigned(p,proj):computesigned(p,proj);
	if (dis==Contour::UNDEF) return dis;
	if (iol) {
		if (dis<0)
			printdirectedseg(iol,p,proj,A3dColor(1,1,.3));
		else
			printdirectedseg(iol,proj,p,A3dColor(1,1,1));
	}
	return dis;
}

static void printdirectedseg(Mk3d* mk, const Point& p1, const Point& p2,
			     const A3dColor& col)
{
	Vector v=p2-p1;
	assertx(v.normalize());
	mk->diffuse(col); mk->specular(0,0,0); mk->phong(3);
	mk->point(p1*xformi);
	mk->normal(v);
	mk->point(p2*xformi);
	mk->normal(v);
	mk->endpolyline();
}

// Find the unsigned distance to the centroid of the k nearest data points.
static double computeunsigned(const Point& p, Point& proj)
{
	SpatialSearch ss(*SPp,p);
	Homogeneous h;
	const int k=1;		// ?? parameter
	for (int i=0;i<k;i++) {
		assertx(!ss.done());
		int pi=Conv<int>::d(ss.next());
		h+=co[pi];
	}
	proj=toPoint(h/k);
	return dist(p,proj)-unsigneddis;
}

// Find the closest tangent plane origin and compute the signed distance to
// that tangent plane.
// Was: check to see if the projection onto the tangent plane lies farther
// than samplingd from any data point.
// Now: check to see if the sample point is farther than samplingd+cube_size
// from any data point.
static double computesigned(const Point& p, Point& proj)
{
	int tpi;
	NEST {
		SpatialSearch ss(*SPpc,p);
		tpi=Conv<int>::d(ss.next());
	}
	Vector vptopc=p-pcorg[tpi];
	double dis=dot(vptopc,pcnor[tpi]);
	proj=p-dis*pcnor[tpi];
	if (!is3D) assertx(!proj[0]);
	if (is3D && (proj[0]<=0 || proj[0]>=1) ||
	    proj[1]<=0 || proj[1]>=1 || proj[2]<=0 || proj[2]>=1)
		return Contour::UNDEF;
	if (1) {
		// check that projected point is close to a data point
		SpatialSearch ss(*SPp,proj);
		double dis2; (void)Conv<int>::d(ss.next(&dis2));
		if (dis2>square(samplingd)) return Contour::UNDEF;
	}
	if (prop) {
		// check that grid point is close to a data point
		SpatialSearch ss(*SPp,p);
		double dis2; (void)Conv<int>::d(ss.next(&dis2));
		double griddiagonal2=square(1./gridsize)*3.;
		const double fudge=1.2;
		if (dis2>griddiagonal2*square(fudge)) return Contour::UNDEF;
	}
	return dis;
}

static void outputborder(const Polygon& poly)
{
	assertx(ioC);
	static A3dElem el;
	if (is3D) {
		el.init(A3dElem::TPolygon);
	} else {
		el.init(A3dElem::TPolyline);
	}
	for (int i=0;i<poly.num();i++)
		el+=A3dVertex(poly[i]*xformi,Vector(0,0,0),
			      A3dVertexColor(A3dColor(1,.16,0),
					     A3dColor(1,.5,.3),
					     A3dColor(3,0,0)));
	ioC->oa3d().write(el);
}

static void outputcontour(const Polygon& poly)
{
	assertx(iom);
	static A3dElem el;
	el.init(is3D?A3dElem::TPolygon:A3dElem::TPolyline);
	for (int i=0;i<poly.num();i++)
		el+=A3dVertex(poly[i]*xformi,Vector(0,0,0),
			      A3dVertexColor(A3dColor(.7,.3,.1)));
	iom->oa3d().write(el);
}

static void prtgraph(Mk3d* mk, const Graph<int>* g, const Point* pa,
		     const Vector* pn)
{
	for (int i=0;i<num;i++) {
		for (GraphEdges<int> ge(*g,i);ge;ge.next()) {
			int j=ge();
			if (j>=i) continue;
			mk->point(pa[i]*xformi); if (pn) mk->normal(pn[i]);
			mk->point(pa[j]*xformi); if (pn) mk->normal(pn[j]);
			mk->endpolyline();
		}
	}
}
