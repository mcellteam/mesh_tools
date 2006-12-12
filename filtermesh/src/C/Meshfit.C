// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include <csignal>
#include "Meshfit.h"
#include "Options.h"
#include "GMesh.h"
#include "MeshOp.h"
#include "A3dStream.h"
#include "FileIO.h"
#include "GeomOp.h"		// for DihedralAngleCos()
#include "Bbox.h"
#include "Facedistance.h"
#include "FrameIO.h"
#include "Set.h"
#include "Array.h"
#include "Map.h"
#include "HashStruct.h"
#include "Stack.h"
#include "Random.h"
#include "Recipes.h"		// LineOptimization()
#include "Timer.h"

#include "fdstream.hpp"
#include <iostream>
#include <iomanip>
using std::cout;
using std::ostream;
using std::flush;

// Notes:
//  Could check dihedral(starbar) instead of dihedral(star)
//   using GetEdgeRing() ->vertices

// SIGUSR1: dump current mesh out to file and continue
// SIGUSR2: exit()

static GMesh mesh;
static DataPts pt;

static double spring=0;
static double spbf=1;
static double fliter=1;
static double crep=1e-5;
static double crbf=3;
static int nooutput=0;
static int verb=1;
static double feswaasym=.01;

static WSA3dStream oa3d(cout);
static int sdebug=GetenvValue("DEBUG");
static Frame xform;		// original mesh+pts -> mesh+pts in unit cube
static Frame xformi;		// inverse
enum Operation { OP_ecol, OP_espl, OP_eswa, OPNUM };
const char* opname[OPNUM]={"ecol", "espl", "eswa"};
enum OResult { OR_success, OR_energy, OR_dih, OR_illegal, ORNUM };
const char* orname[ORNUM]={"success", "positive_energy", "bad_dihedral",
			   "illegal_move"};
static struct {
	int na[OPNUM], ns[OPNUM], nor[ORNUM];
	int notswaps;
} opstat;

static double mincos=-1./3.;	// acos(109.471) == tetrahedron angle
const double sprsche[]={1e-2,1e-3,1e-4,1e-8,-1};
const int MAXGFITITER=30;
static FILE* filespawn=0;
static ostream* ospawn=0;

static Univ hashe(const MEdge* e) {
	return Conv<int>::e(mesh.vertexid(mesh.vertex1(const_cast<Edge>(e)))+
			    mesh.vertexid(mesh.vertex2(const_cast<Edge>(e)))*76541);
}
static int cmpe(const MEdge* e1, const MEdge* e2) { return e1!=e2; }
// Set of candidate edges in stoc
// HashStruct so that hashing does not use pointer values (portable random)
static HashStruct<MEdge> ecand(hashe,cmpe);

// Declarations
static int dofilename(int argc, char const** argv);
static int domfilename(int argc, char const** argv);
static int dooutmesh(int argc, char const** argv);
static int dospawn(int argc, char const** argv);
static int dogfit(int argc, char const** argv);
static int dofgfit(int argc, char const** argv);
static void dostoc();
static int dolfit(int argc, char const** argv);
static void perhapsinitialize();
static void computexform();
static void meshtransform(const Frame& f);
static void initialprojection();
static void doreconstruct();
static void dosimplify();
static void dofour1split();
static void dopclp();
static void dorecord();
// stoc
static OResult tryop(Edge e, Operation op, double& edrss);
static OResult tryecol(Edge e, int ni, int nri, double& edrss);
static OResult tryespl(Edge e, int ni, int nri, double& edrss);
static OResult tryeswa(Edge e, int ni, int nri, double& edrss);
static OResult checkhalfeswa(Edge e, Vertex vo1, Vertex v1,
			     Vertex vo2, Vertex v2, Face f2, Face f1,
			     int ni, int nri, double& edrss);
static void cleanupneighborhood(Vertex v, int nri);
//
static void sigUSR1();
static void sigUSR2();
// other
static void analyzemesh(const char* s);
static double getedis();
static double getespr();
static double geterep();
static int getnbv();
static double springenergy(Vertex v1, Vertex v2);
static double springenergy(Edge e);
static void pushfacepoints(Face f, Stack<int>& stpts);
static void removeface(Face f);
// unordered

int main(int argc, char const** argv)
{
	SHOWDF("Created using\n%s",CreateHeader(argc,argv));
	Options opts;
	OPTSDC(opts,filename,	"name : point data (can be -)");
	OPTSDC(opts,mfilename,	"name : initial mesh (can be -)");
	OPTSPC(opts,crep,	"v: set constant for representation energy");
	OPTSDC(opts,reconstruct,": apply surface reconstruction schedule");
	OPTSDC(opts,simplify,	": apply mesh simplification schedule");
	opts.c("",":");
	OPTSPC(opts,spring,	"tension : set spring constant");
	opts.c("",":");
	OPTSDC(opts,gfit,	"niter : do global fit (0=until convergence)");
	OPTSDC(opts,fgfit,	"niter : use conjugate gradients");
	OPTSDC(opts,stoc,	": do local stochastic mesh operations");
	OPTSDC(opts,lfit,	"ni nli : do ni iters, each nli local fits");
	OPTSDC(opts,four1split,	": do global four-to-one split");
	OPTSDC(opts,outmesh,	"filename : output current mesh to file");
	opts.c("",":");
	OPTSDC(opts,pclp,	": print projections onto mesh (lines)");
	OPTSDC(opts,record,	": print mesh changes on cout, -noout");
	OPTSDC(opts,spawn,	"'command': send record to popen");
	OPTSFC(opts,nooutput,	": don't print final mesh on stdout");
	OPTSPC(opts,verb,	"i : verbosity level (1=avg,2=more,3=lots)");
	opts.c("",":");
	OPTSPC(opts,crbf,	"ratio : set repr. energy boundary factor");
	OPTSPC(opts,spbf,	"ratio : set spring constant boundary factor");
	OPTSPC(opts,fliter,	"factor : modify # local iters done in stoc");
	OPTSPC(opts,feswaasym,	"f : set drss threshold (fraction of crep)");
	signal(SIGUSR1,HHSIG_PF(sigUSR1));
	signal(SIGUSR2,HHSIG_PF(sigUSR2));
	TIMER(Meshfit);
	opts.allmustparse();
	if (!opts.parse(argc,argv)) { opts.problem(argc,argv); return 1; }
	perhapsinitialize();
	SHOWDF("\n");
	analyzemesh("FINAL");
	ETIMER(Meshfit);
	CleanUp();
	SHOWDF("\n"),SHOWDF("EndMeshfit\n");
	if (!nooutput) { meshtransform(xformi); mesh.write(cout); }
	if (filespawn) { delete ospawn; pclose(filespawn); }
	mesh.clear();
	pt.clear();
	return 0;
}

static int domfilename(int argc, char const** argv)
{
	TIMER(_mfilename);
	assertx(argc>1);
	assertx(!mesh.numVertices());
	RFile is(argv[1]);
	mesh.read(is());
	return 2;
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
		pt.enter(el[0].p);
	}
	SHOWDF("%d points read\n",pt.n);
	return 2;
}

static void perhapsinitialize()
{
	assertx(pt.n && mesh.numVertices());
	assertw1(spring>0);	// just warn user
	if (pt.cmf[0]) return;	// already initialized
	computexform();
	initialprojection();
}

static void computexform()
{
	Bbox bb;
	bb.clear();
	NEST { for (int i=0;i<pt.n;i++) bb.takeunion(pt.co[i]); }
	ForMeshVertex(mesh,v) {
		bb.takeunion(mesh.point(v));
	} EndFor;
	xform=bb.getFrameToSmallCube();
	SHOWDF("Applying xform: %s",FrameIO::string(xform,1,0));
	xformi=~xform;
	NEST { for (int i=0;i<pt.n;i++) pt.co[i]*=xform; }
	meshtransform(xform);
}

static void meshtransform(const Frame& f)
{
	ForMeshVertex(mesh,v) {
		mesh.setPoint(v,mesh.point(v)*f);
	} EndFor;
}

static void initialprojection()
{
	ForMeshFace(mesh,f) {
		pt.mfpts.enter(f,new Set<int>);
	} EndFor;
	NEST {
		TIMER(_initialproj);
		GlobalProject(mesh,pt);
	}
	if (sdebug) pt.OK();
	analyzemesh("INITIAL");
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
	TIMER(reconstruct);
	if (!spring) spring=sprsche[0];
	perhapsinitialize();
	const char* argv[2]; argv[1]="20";
	dofgfit(2, argv); // -fgfit 20
	applyschedule();
}

static void dosimplify()
{
	TIMER(simplify);
	if (!spring) spring=sprsche[0];
	perhapsinitialize();
	applyschedule();
}

static int dogfit(int argc, char const** argv)
{
        int i;
	perhapsinitialize();
	TIMER(_gfit);
	assertx(argc>1); int niter=atoi(argv[1]);
	SHOWDF("\n");
	SHOWDF("Beginning gfit, %d iterations, spr=%g\n",niter,spring);
	double ecsc=getedis()+getespr();	// energy constant simplicial complex
//	for (int i=0;!niter || i<niter;) {
	for (i=0;!niter || i<niter;) {
		if (CheckSignal()) break;
		if (!niter && i>=MAXGFITITER) break;
		i++;
		ATIMER(__gfit_iter);
		SHOWF("iter %d/%d\n",i,niter),cout << flush;
		NEST { ATIMER(__glls); GlobalFit(mesh,pt,spring,spbf); }
		NEST { ATIMER(__gproject); GlobalProject(mesh,pt); }
		if (sdebug) pt.OK(),mesh.OK();
		double edis=getedis(), espr=getespr(), necsc=edis+espr;
		double echange=necsc-ecsc;
		assertw(echange<0);
		SHOWF(" edis=%g espr=%g  necsc=%g  echange=%g\n",
		      edis,espr,necsc,echange);
		if (!niter && echange>-1e-4) break;
		ecsc=necsc;
	}
	SHOWDF("Finished gfit, did %d iterations\n",i);
	analyzemesh("after_gfit");
	return 2;
}

static double dot(const Array<Vector>& a, const Array<Vector>&b)
{
	assertx(a.num() && a.num()==b.num());
	double v=0; ForIndex(i,a.num()) { v+=dot(a[i],b[i]); } EndFor; return v;
}

class FG {
  public:
	FG(GMesh& pmesh);
	~FG();
	void computelinedir();
	double computeonline(double t);
	double computedonline(); // assumes previous t!
  private:
	GMesh& mesh;
	int nv;
	int first;
	Map<Vertex,int> mvi; Array<Vertex> iv;
	Array<Point> opos;
	Array<Vector> dir, g, h;
	void computegrad(Array<Vector>& grad);
	DISABLECOPY(FG);	// for GNUG 2.5.8
};

static FG* fg;

FG::FG(GMesh& pmesh)
: mesh(pmesh), nv(mesh.numVertices()), first(1),
  opos(nv), dir(nv), g(nv), h(nv)
{
	ForMeshVertex(mesh,v) { mvi.enter(v,iv.num()); iv+=v; } EndFor;
}

FG::~FG() { }

void FG::computegrad(Array<Vector>& grad)
{
	assertx(grad.num()==nv);
	ForIndex(j,nv) { grad[j]=Vector(0,0,0); } EndFor;
	ForIndex(i,pt.n) {
		Vertex va[3]; mesh.vertices(pt.cmf[i],va);
		Vector vtop=pt.co[i]-pt.clp[i];
		ForIndex(k,3) {
			grad[mvi.get(va[k])]+=vtop*(-2*pt.bary[i][k]);
		} EndFor;
	} EndFor;
	ForMeshVertex(mesh,v) {
		ForVertexEdge(mesh,v,e) {
			Vertex vv=mesh.oppVertex(v,e);
			Vector vtovv=mesh.point(vv)-mesh.point(v);
			double sp=mesh.isBoundary(e)?spring*spbf:spring;
			grad[mvi.get(v)]+=vtovv*(-2*sp);
		} EndFor;
	} EndFor;
}

void FG::computelinedir()
{
	// Record current vertex positions
	ForIndex(j,nv) { opos[j]=mesh.point(iv[j]); } EndFor;
	// modeled after numerical recipes frprmn()
	Array<Vector> grad(nv); computegrad(grad);
	if (first) {
		first=0;
		ForIndex(j,nv) { dir[j]=h[j]=g[j]=-grad[j]; } EndFor;
	} else {
		double gg=dot(g,g), dgg=0;
		const int fletcher=0;
		if (fletcher) {	// Fletcher-Reeves
			dgg=dot(grad,grad);
		} else {	// Polak-Ribiere
			ForIndex(j,nv) {
				dgg+=dot(grad[j]+g[j],grad[j]);
			} EndFor;
		}
		double gam=dgg/assertv(gg);
		ForIndex(j,nv) {
			g[j]=-grad[j];
			dir[j]=h[j]=g[j]+gam*h[j];
		} EndFor;
	}
}

double FG::computeonline(double t)
{
	ForIndex(j,nv) { mesh.setPoint(iv[j],opos[j]+dir[j]*t); } EndFor;
	NEST { ATIMER(__fgproject); GlobalProject(mesh,pt); }
	double e=getedis()+getespr();
	// SHOWF("computeonline t=%g e=%g\n",t,e);
	return e;
}

double FG::computedonline()
{
	double v; {		// GNUG 2.5.8
		Array<Vector> grad(nv); computegrad(grad); v=dot(grad,dir);
	} return v;
}

static double fgeval(double t) { return fg->computeonline(t); }
static double fgevald(double) { return fg->computedonline(); }

static int dofgfit(int argc, char const** argv)
{
	perhapsinitialize();
	TIMER(_fgfit);
	assertx(argc>1); int niter=atoi(argv[1]);
	SHOWDF("\n");
	SHOWDF("Beginning fgfit, %d iterations, spr=%g\n",niter,spring);
	double ecsc=getedis()+getespr();	// energy constant simplicial complex
	SHOWDF(" initial ecsc=%g\n",ecsc);
	fg=new FG(mesh);
	double oldtmin=.05;
	assertx(niter);
	ForIndex(ni,niter) {
		if (CheckSignal()) break;
		ATIMER(__fgfit_iter);
		fg->computelinedir();
		const double ftol=.05;
		double t1=oldtmin*.9, t2=oldtmin*1.1+.001;
		double tmin, emin; int neval;
		LineOptimization(&fgeval,fgevald,t1,t2,ftol,tmin,emin,neval);
		SHOWF("iter %d/%d  neval=%d tmin=%9.6lf e=%lg\n",
		      ni+1,niter,neval,tmin,emin);
		oldtmin=tmin;
		SSTAT(Sneval,neval); SSTAT(Stmin,tmin);
	} EndFor;
	delete fg,fg=0;
	analyzemesh("after_fgfit");
	return 2;
}

static void dostoc()
{
	perhapsinitialize();
	TIMER(_stoc);
	SHOWDF("\n");
	SHOWDF("Beginning stoc, spring=%g, fliter=%g\n",spring,fliter);
	NEST { for (int c=0;c<OPNUM;c++) opstat.na[c]=opstat.ns[c]=0; }
	NEST { for (int c=0;c<ORNUM;c++) opstat.nor[c]=0; }
	opstat.notswaps=0;
	int ni=0, nbad=0, lasti=-99999;
	ForMeshEdge(mesh,e) {
		ecand.enter(e);
	} EndFor;
	for (;!ecand.empty();) {
		if (CheckSignal()) break;
		ni++;
		HashStructIter<MEdge> hi(ecand,Random::G);
		Edge e=hi();
		assertx(ecand.remove(e));
		mesh.valid(e); // optional make sdebug
		Operation op; OResult oor=OR_illegal; double edrss;
		if (oor!=OR_success) op=OP_ecol,oor=tryop(e,op,edrss);
		if (oor!=OR_success) op=OP_espl,oor=tryop(e,op,edrss);
		if (oor!=OR_success) op=OP_eswa,oor=tryop(e,op,edrss);
		if (oor!=OR_success) { nbad++; continue; }
		if (verb>=2 || verb>=1 && ni>=lasti+100) {
			SHOWF("it %5d, %s (after %3d) [%5d/%-5d] edrss=%e\n",
			       ni,opname[op],nbad,
			       ecand.num(),mesh.numEdges(),edrss);
			lasti=ni;
		}
		if (ospawn) *ospawn << flush;
		nbad=0;
	}
	SHOWDF("it %d, last search: %d wasted attempts\n",ni,nbad);
	int nat=0, nst=0;
	ForIndex(i,OPNUM) { nat+=opstat.na[i]; nst+=opstat.ns[i]; } EndFor;
	SHOWDF("%s(ecol=%d/%d, espl=%d/%d eswa=%d/%d tot=%d/%d)\n",
	       "Endstoc: ",
	       opstat.ns[OP_ecol],opstat.na[OP_ecol],
	       opstat.ns[OP_espl],opstat.na[OP_espl],
	       opstat.ns[OP_eswa],opstat.na[OP_eswa], nst,nat);
	SHOWDF("         (otswaps=%d)\n",opstat.notswaps);
	SHOWDF("Result of %d attempted operations:\n",nat);
	ForIndex(i,ORNUM) {
		SHOWDF("  %5d %s\n",opstat.nor[i],orname[i]);
	} EndFor;
	analyzemesh("after_stoc");
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
		if (CheckSignal()) break;
		ForMeshVertex(mesh,v) {
			if (CheckSignal()) break;
			FitRing(mesh,pt,spring,spbf,mincos,v,nli);
		} EndFor;
	}
	SHOWDF("Finished lfit\n");
	analyzemesh("after_lfit");
	return 3;
}

static void dofour1split()
{
        int i;
	perhapsinitialize();
	TIMER(_four1split);
	static int level=0;
	if (!level) {
		ForMeshVertex(mesh,v) {
			mesh.setInfo(v,new
				     MeshInfoString(hform(" %d ",level)));
		} EndFor;
	}
	level++;
	Stack<Edge> stacke;
	ForMeshEdge(mesh,e) { stacke.push(e); } EndFor;
	Stack<Face> stackf;
	ForMeshFace(mesh,f) { stackf.push(f); } EndFor;
	Map<Edge,Vertex> menewv; // old Edge -> Vertex
	// Create new vertices and compute their positions
	ForStack(stacke,Edge,e) {
		Vertex v=mesh.createVertex();
		menewv.enter(e,v);
		mesh.setPoint(v,interp(mesh.point(mesh.vertex1(e)),
				       mesh.point(mesh.vertex2(e))));
		mesh.setInfo(v,new MeshInfoString(hform(" %d ",level)));
	} EndFor;
	// Subdivide faces and update projections
	ForStack(stackf,Face,f) {
		Vertex va[3], vs[3];
		mesh.vertices(f,va);
//		for (int i=0;i<3;i++)
		for (i=0;i<3;i++)
			vs[i]=menewv.get(mesh.edge(va[i],va[(i+1)%3]));
		Stack<int> stpts;
		Stack<Face> stfaces;
		pushfacepoints(f,stpts);
		for (i=0;i<3;i++) {
			Face f=assertv(mesh.createFace(va[i],vs[i],
						       vs[(i+2)%3]));
			stfaces.push(f); pt.mfpts.enter(f,new Set<int>);
		}
		NEST {
			Face f=mesh.createFace(vs,3);
			assertx(f);
			stfaces.push(f); pt.mfpts.enter(f,new Set<int>);
		}
		ReprojectLocally(mesh,pt,stpts,stfaces);
		Set<int>* set=pt.mfpts.get(f);
		assertx(set->empty());
		delete set;
	} EndFor;
	ForStack(stackf,Face,f) {
		mesh.destroyFace(f);
	} EndFor;
	if (sdebug) pt.OK(),assertx(mesh.isNice());
}

static int dooutmesh(int argc, char const** argv)
{
	perhapsinitialize();
	assertx(argc>1);
	WFile os(argv[1]);
	meshtransform(xformi);
	mesh.write(os());
	meshtransform(xform);
	return 2;
}

static void dopclp()
{
	nooutput=1;
	perhapsinitialize();
	A3dElem el;
	for (int i=0;i<pt.n;i++) {
		el.init(A3dElem::TPolyline);
		Point pco=pt.co[i], pclp=pt.clp[i];
		pco*=xformi, pclp*=xformi;
		el+=A3dVertex(pco,Vector(0,0,0),A3dVertexColor::White);
		el+=A3dVertex(pclp,Vector(0,0,0),A3dVertexColor::White);
		oa3d.write(el);
	}
	oa3d.flush();
}

static void dorecord()
{
	// xform not undone!
	mesh.recordChanges(&cout);
	nooutput=1;
}

static int dospawn(int argc, char const** argv)
{
	// xform not undone!
	filespawn=assertv(popen(argv[1],"w"));
	ospawn=new boost::fdostream(fileno(filespawn));
	mesh.write(*ospawn);
	mesh.recordChanges(ospawn);
	return 2;
}

//*** stoc

static OResult tryop(Edge e, Operation op, double& edrss)
{
	if (CheckSignal()) return OR_illegal;
	ATIMER(__tryop);
	OResult oor;
	oor=(op==OP_ecol?tryecol(e,int(4*fliter+.5),int(2*fliter+.5),edrss):
	     op==OP_espl?tryespl(e,int(3*fliter+.5),int(4*fliter+.5),edrss):
	     op==OP_eswa?tryeswa(e,int(3*fliter+.5),int(2*fliter+.5),edrss):
	     (assertnever(""),OR_success));
	opstat.na[op]++;
	if (oor==OR_success) opstat.ns[op]++;
	opstat.nor[oor]++;
	return oor;
}

static OResult tryecol(Edge e, int ni, int nri, double& edrss)
{
	ATIMER(__tryecol);
	if (!mesh.niceEdgeCollapse(e)) return OR_illegal; // not a legal move
	Vertex v1=mesh.vertex1(e), v2=mesh.vertex2(e);
	Face f1=mesh.face1(e), f2=mesh.face2(e);
	Array<const Point*> wa;
	GatherEdgeRing(mesh,e,wa);
	int nbvb=mesh.isBoundary(v1)+mesh.isBoundary(v2);
	double minb=min(MinDihedralAboutVertex(mesh,v1),
		       MinDihedralAboutVertex(mesh,v2));
	double rssf=0;
	Stack<int> stpts;
	Stack<Face> stfaces;
	ForVertexFace(mesh,v1,f) {
		pushfacepoints(f,stpts);
		if (f==f1 || f==f2) continue;
		stfaces.push(f);
	} EndFor;
	ForVertexFace(mesh,v2,f) {
		if (f==f1 || f==f2) continue;
		pushfacepoints(f,stpts);
		stfaces.push(f);
	} EndFor;
	ForStack(stpts,int,pi) {
		rssf+=dist2(pt.co[pi],pt.clp[pi]);
	} EndFor;
	ForVertexVertex(mesh,v1,v) {
		rssf+=springenergy(v1,v);
	} EndFor;
	ForVertexVertex(mesh,v2,v) {
		if (v==v1) continue;
		rssf+=springenergy(v2,v);
	} EndFor;
	// Find the best starting location by exploring one iteration.
	double minrss1=1e30; int minii=-1;
	for (int ii=0;ii<3;ii++) {
		Point newp=interp(mesh.point(v1),mesh.point(v2),ii*.5);
		double rss0,rss1;
		LocalFit(pt,spring,spbf,stpts,wa,1,newp,rss0,rss1);
		double mina=MinLocalDihedral(wa,newp);
		if (mina<mincos && mina<minb) continue; // change disallowed
		if (rss1<minrss1) minrss1=rss1,minii=ii;
	}
	if (minii<0) return OR_dih; // no dihedrally admissible configuration
	// Then, explore ni iterations from that chosen starting point
	Point newp=interp(mesh.point(v1),mesh.point(v2),minii*.5);
	double rss0,rss1;
	LocalFit(pt,spring,spbf,stpts,wa,ni,newp,rss0,rss1);
	double mina=MinLocalDihedral(wa,newp);
	if (mina<mincos && mina<minb) return OR_dih; // change disallowed
	double drss=rss1-rssf-(nbvb==2?crbf:1)*double(crep);
	edrss=drss;
	if (verb>=3) SHOWF("# ecol: rssf=%lg rss1=%lg drss=%lg\n",
			   rssf,rss1,drss);
	if (drss>=0) return OR_energy; // energy function does not decrease
	// ALL SYSTEMS GO
	ATIMER(__doecol);
	ForEdgeVertex(mesh,e,v) {
		ForVertexEdge(mesh,v,e) { // one duplication
			ecand.remove(e);
		} EndFor;
	} EndFor;
	removeface(f1);
	if (f2) removeface(f2);
	mesh.collapseEdge(e);	// v1 kept
	// add about 12-16 edges
	ForVertexEdge(mesh,v1,e) {
		ecand.add(e);
	} EndFor;
	ForVertexFace(mesh,v1,f) {
		ecand.add(mesh.oppEdge(v1,f));
	} EndFor;
	if (sdebug>=2) assertx(mesh.isNice());
	mesh.setPoint(v1,newp);
	ReprojectLocally(mesh,pt,stpts,stfaces);
	cleanupneighborhood(v1,nri);
	return OR_success;
}

static OResult tryespl(Edge e, int ni, int nri, double& edrss)
{
	ATIMER(__tryespl);
	// always legal
	Vertex v1=mesh.vertex1(e), v2=mesh.vertex2(e);
	Face f1=mesh.face1(e), f2=mesh.face2(e);
	Vertex vo1=mesh.sideVertex1(e), vo2=mesh.sideVertex2(e);
	Array<const Point*> wa;
	wa+=&mesh.point(v2);
	wa+=&mesh.point(vo1);
	wa+=&mesh.point(v1);
	if (vo2) wa+=&mesh.point(vo2),wa+=wa[0];
	double minb=vo2?EdgeDihedralAngleCos(mesh,e):2;
	double rssf=springenergy(e);
	Stack<int> stpts;
	ForEdgeFace(mesh,e,f) {
		pushfacepoints(f,stpts);
	} EndFor;
	ForStack(stpts,int,pi) {
		rssf+=dist2(pt.co[pi],pt.clp[pi]);
	} EndFor;
	Point newp=interp(*wa[0],*wa[2]);
	double rss0,rss1;
	LocalFit(pt,spring,spbf,stpts,wa,ni,newp,rss0,rss1);
	double mina=MinLocalDihedral(wa,newp);
	if (mina<mincos && mina<minb) return OR_dih; // change disallowed
	double drss=rss1-rssf+(vo2?1:crbf)*double(crep);
	edrss=drss;
	if (verb>=3) SHOWF("# espl: rssf=%lg rss1=%lg drss=%lg\n",
			   rssf,rss1,drss);
	if (drss>=0) return OR_energy; // energy function does not decrease
	// ALL SYSTEMS GO
	ATIMER(__doespl);
	ForEdgeFace(mesh,e,f) {
		ForFaceEdge(mesh,f,e) { // one duplication
			ecand.remove(e);
		} EndFor;
	} EndFor;
	Vertex v=mesh.splitEdge(e);
	mesh.setPoint(v,newp);
	// add 8 edges (5 if boundary)
	ForVertexFace(mesh,v,f) {
		ForFaceEdge(mesh,f,e) { // four duplications
			ecand.add(e);
		} EndFor;
	} EndFor;
	ForVertexFace(mesh,v,f) {
		if (f!=f1 && f!=f2) pt.mfpts.enter(f,new Set<int>);
	} EndFor;
	// Since stpts project onto f1+f2 (which are still there),
	// it is easy to update the projections:
	FitRing(mesh,pt,spring,spbf,mincos,v,2);
	cleanupneighborhood(v,nri);
	return OR_success;
}

static OResult tryeswa(Edge e, int ni, int nri, double& edrss)
{
	ATIMER(__tryeswa);
	if (!mesh.legalEdgeSwap(e)) return OR_illegal; // not legal move
	Vertex v1=mesh.vertex1(e), v2=mesh.vertex2(e);
	Face f1=mesh.face1(e), f2=mesh.face2(e);
	Vertex vo1=mesh.sideVertex1(e), vo2=mesh.sideVertex2(e);
	// Compare angles immediately before and after swap
	double minb=EdgeDihedralAngleCos(mesh,e);
	double mina=DihedralAngleCos(mesh.point(vo1),mesh.point(vo2),
				    mesh.point(v1),mesh.point(v2));
	if (mina<mincos && mina<minb) return OR_dih;
	// Will try both cases, but randomly select which one to try first
	OResult oor;
	if (Random::G.getint()%2) {
		oor=checkhalfeswa(e,vo1,v1,vo2,v2,f2,f1,ni,nri,edrss);
		if (oor==OR_success) return oor;
		oor=checkhalfeswa(e,vo2,v2,vo1,v1,f1,f2,ni,nri,edrss);
	} else {
		oor=checkhalfeswa(e,vo2,v2,vo1,v1,f1,f2,ni,nri,edrss);
		if (oor==OR_success) return oor;
		oor=checkhalfeswa(e,vo1,v1,vo2,v2,f2,f1,ni,nri,edrss);
	}
	return oor;
}

// Note: correspondences on Edge e such as vertex1(e)==v1 may fail here!
// Try swapping edge (v1,v2) into edge (vo1,vo2),
//  allowing vertex vo1 to move.
// To do this, gather points in current ring of vo1 plus points on f2.
static OResult checkhalfeswa(Edge e, Vertex vo1, Vertex v1,
			     Vertex vo2, Vertex v2, Face f2, Face f1,
			     int ni, int nri, double& edrss)
{
	Array<const Point*> wa;
	NEST {
		assertx(mesh.ccwVertex(vo1,v1)==v2);
		Vertex w=mesh.mostClwVertex(vo1), wf=w;
		for (;;) {
			wa+=&mesh.point(w);
			if (w==v1) wa+=&mesh.point(vo2);
			w=mesh.ccwVertex(vo1,w);
			if (!w || w==wf) break;
		}
		if (w) wa+=wa[0];
	}
	Stack<int> stpts;
	Stack<Face> stfaces;
	ForVertexFace(mesh,vo1,f) {
		if (f==f1) continue;
		assertx(f!=f2);
		stfaces.push(f);
		pushfacepoints(f,stpts);
	} EndFor;
	pushfacepoints(f1,stpts);
	pushfacepoints(f2,stpts);
	double rssf=0;
	ForStack(stpts,int,pi) {
		rssf+=dist2(pt.co[pi],pt.clp[pi]);
	} EndFor;
	ForVertexVertex(mesh,vo1,v) {
		rssf+=springenergy(vo1,v);
	} EndFor;
	rssf+=springenergy(e);
	Point newp=mesh.point(vo1);
	double minb=min(MinDihedralAboutVertex(mesh,vo1),
		       EdgeDihedralAngleCos(mesh,e));
	double rss0,rss1;
	LocalFit(pt,spring,spbf,stpts,wa,ni,newp,rss0,rss1);
	double mina=MinLocalDihedral(wa,newp);
	if (mina<mincos && mina<minb) return OR_dih; // change disallowed
	double drss=rss1-rssf+crep*feswaasym;
	edrss=drss;
	if (verb>=3) SHOWF("# eswa: rssf=%lg rss1=%lg drss=%lg\n",
			   rssf,rss1,drss);
	if (drss>0) return OR_energy;
	// ALL SYSTEMS GO
	ATIMER(__doeswa);
	ForEdgeFace(mesh,e,f) {
		ForFaceEdge(mesh,f,e) {
			ecand.remove(e);
		} EndFor;
	} EndFor;
	removeface(f1);
	removeface(f2);
	Edge enew=assertv(mesh.swapEdge(e));
	mesh.setPoint(vo1,newp);
	// add about 9 edges
	ecand.add(mesh.edge(vo2,v1));
	ecand.add(mesh.edge(vo2,v2));
	ForVertexEdge(mesh,vo1,e) {
		ecand.add(e);
	} EndFor;
	ForEdgeFace(mesh,enew,f) {
		stfaces.push(f);
		pt.mfpts.enter(f,new Set<int>);
	} EndFor;
	ReprojectLocally(mesh,pt,stpts,stfaces);
	FitRing(mesh,pt,spring,spbf,mincos,vo2,nri);
	FitRing(mesh,pt,spring,spbf,mincos,v1,nri);
	FitRing(mesh,pt,spring,spbf,mincos,v2,nri);
	FitRing(mesh,pt,spring,spbf,mincos,vo1,nri);
	cleanupneighborhood(vo1,nri);
	FitRing(mesh,pt,spring,spbf,mincos,vo2,nri);
	FitRing(mesh,pt,spring,spbf,mincos,v1,nri);
	FitRing(mesh,pt,spring,spbf,mincos,v2,nri);
	FitRing(mesh,pt,spring,spbf,mincos,vo1,nri);
	return OR_success;
}

static void cleanupneighborhood(Vertex v, int nri)
{
	if (nri) {
		ForVertexVertex(mesh,v,w) {
			FitRing(mesh,pt,spring,spbf,mincos,w,nri);
		} EndFor;
		FitRing(mesh,pt,spring,spbf,mincos,v,nri);
	}
}

//*** other

static void analyzemesh(const char* s)
{
	int nv=mesh.numVertices(), nf=mesh.numFaces(), ne=mesh.numEdges();
	int nbv=getnbv();
	double edis=getedis(), espr=getespr(), erep=geterep();
	double etot=edis+espr+erep;
	SHOWDF("%-12s: mesh v=%d f=%d e=%d (nbv=%d)\n",s,nv,nf,ne,nbv);
	SHOWDF("  energies: edis=%g espr=%g erep=%g etot=%g\n",
	       edis,espr,erep,etot);
}

static double getedis()
{
	double edis=0;
	for (int i=0;i<pt.n;i++) edis+=dist2(pt.co[i],pt.clp[i]);
	return edis;
}

static double getespr()
{
	double espr=0;
	ForMeshEdge(mesh,e) {
		espr+=springenergy(e);
	} EndFor;
	return espr;
}

static double geterep()
{
	return crep*(mesh.numVertices()+getnbv()*(crbf-1));
}

static int getnbv()
{
	int nbv=0;
	ForMeshVertex(mesh,v) {
		if (mesh.isBoundary(v)) nbv++;
	} EndFor;
	return nbv;
}

static double springenergy(Vertex v1, Vertex v2)
{
	return dist2(mesh.point(v1),mesh.point(v2))*
		(mesh.isBoundary(mesh.edge(v1,v2))?spring*spbf:spring);
}

static double springenergy(Edge e)
{
	return dist2(mesh.point(mesh.vertex1(e)),mesh.point(mesh.vertex2(e)))*
		(mesh.isBoundary(e)?spring*spbf:spring);
}

static void pushfacepoints(Face f, Stack<int>& stpts)
{
	ForSet(*pt.mfpts.get(f),int,pi) {
		stpts.push(pi);
	} EndFor;
}

// Make face's points project nowhere and remove from pt
static void removeface(Face f)
{
	Set<int>& set=*pt.mfpts.get(f);
	for (;!set.empty();) PointChangeFace(pt,set.getone(),0);
	delete assertv(pt.mfpts.remove(f));
}

//*** DataPts

DataPts::DataPts() : n(0) { }

DataPts::~DataPts() { clear(); }

void DataPts::clear()
{
	ForMapValue(mfpts,Face,Set<int>*,set) { delete set; } EndFor;
	co.clear(); cmf.clear(); clp.clear(); bary.clear(); mfpts.clear();
}

void DataPts::enter(const Point& p)
{
	co += p;
	cmf += static_cast<Face>(0);
	clp.add(1);
	bary.add(1);
	n++;
}

void DataPts::OK()
{
	DTIMER(__ptOK);
	mfpts.OK();
	ForMapValue(mfpts,Face,Set<int>*,set) { set->OK(); } EndFor;
	// ::mesh needed for GNUG 2.5.8
	for (int i=0;i<n;i++) {
		Face f=assertv(cmf[i]); ::mesh.valid(f);
		assertx(mfpts.get(f)->contains(i));
		Point pa[3];
		::mesh.points(f,pa);
		Point tclp;
		ProjectPTri(co[i],pa[0],pa[1],pa[2],0,0,0,&tclp);
		assertx(dist2(tclp,clp[i])<1e-10);
	}
	ForMapKeyValue(mfpts,Face,f,Set<int>*,set) {
		::mesh.valid(f);
		ForSet(*set,int,pi) {
			assertx(pi>=0 && pi<n);
			assertx(cmf[pi]==f);
		} EndFor;
	} EndFor;
	ForMeshFace(::mesh,f) {
		assertx(mfpts.contains(f));
	} EndFor;
}

static int gotsignalUSR1=0;
static int gotsignalUSR2=0;

static void sigUSR1()
{
	signal(SIGUSR1,HHSIG_PF(sigUSR1)); // for ATT unix
	gotsignalUSR1=1;
}

static void sigUSR2() {
	signal(SIGUSR2,HHSIG_PF(sigUSR2)); // for ATT unix
	gotsignalUSR2=1;
}

int CheckSignal()
{
	if (gotsignalUSR2) { SHOWDF("*** sigUSR2 exit\n"); return 1; }
	if (!gotsignalUSR1) return 0;
	gotsignalUSR1=0;
	SHOWS("***starting outmesh");
	NEST {
		const char* argv[2]; argv[1]="v.meshfit.m.Z";
		dooutmesh(2, argv);
	}
	SHOWS("***done with outmesh");
	return 0;
}
