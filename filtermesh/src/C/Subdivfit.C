// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include <csignal>
#include "Options.h"
#include "A3dStream.h"
#include "FileIO.h"
#include "PolygonSpatial.h"
#include "SubMesh.h"
#include "Facedistance.h"
#include "MeshOp.h"
#include "GeomOp.h"
#include "Homogeneous.h"
#include "LLS.h"
#include "Set.h"
#include "Map.h"
#include "Array.h"
#include "Stat.h"
#include "Timer.h"
#include "HashStruct.h"
#include "Recipes.h"		// LineOptimization()
#include "Random.h"
#include "Bbox.h"
#include "FrameIO.h"

#include <iomanip>
#include <iostream>
using std::ostream;
using std::cout;
using std::flush;

// Meshinfo:
//  in SubMesh: at Vertex, store numsharpe
//  in fgfit: at Vertex, store gradient vector
//  in tryop: at Vertex of lmesh, store vi
//            at Vertex of smesh, store cvih

// SIGUSR1: dump current mesh out to file and continue
// SIGUSR2: stop optimization and exit gracefully

// weighta is subdivision mask weight, not averaging mask weight,
//  so it must be adjusted.
static int nsubdiv=2;
static int nolimit=0;
static double selective=0;	// set to 180 to refine near sharp only
static int s222=0;		// use s222 instead of n222 (not C1 for n=3)
static double weighta=0;		// for extraord. interior case, interior weight

static int outn;
static int nooutput=0;
static int ecol=0;
static int gecol=0;
static int esha=0;
static int eswa=0;
static int espl=0;
static double feshaasym=.05;
static double feswaasym=.05;
static double crep;
static double csharp;
static double wcrep;
static double wcsharp;

static WFile* wfrecord=0;
const double mincos=-1./3.;	// acos(109.471) == tetrahedron angle

static Array<Point> co;		// points
static GMesh gmesh;		// current control mesh
static Bbox gbbox;
static Frame xform;
static double xformscale;

// always updated
static Array<double> gdis2;	// squared distance associated with each point

// for stoc
static Map<Face,Set<int>*> mfpts; // Face in gmesh -> Set of point indices
static Array<Face> gcmf;	  // pi -> Face in gmesh
static Array<int> gscmfi;	  // smesh closest face index (in face gcmf)

// for general procedures
static Array<Face> gscmf;	// face point projects to in some_smesh
static Array<Bary> gbary;	// barycentric coordinates in some_smesh
static Array<Point> gclp;	// closest point on some_smesh

enum Operation { OP_ecol, OP_espl, OP_eswa, OP_esha, OPNUM};
const char* opname[OPNUM]={"ecol", "espl", "eswa", "esha"};
enum OResult { OR_success, OR_energy, OR_dih, OR_sharp, OR_illegal, ORNUM };
const char* orname[ORNUM]={"success", "positive_energy", "bad_dihedral",
			   "bad_sharp", "illegal_move"};
static int opstat[OPNUM][ORNUM];

static Univ hashe(const MEdge* e) {
	return Conv<int>::e(gmesh.vertexid(gmesh.vertex1(const_cast<Edge>(e)))+
			    gmesh.vertexid(gmesh.vertex2(const_cast<Edge>(e)))*76541);
}
static int cmpe(const MEdge* e1, const MEdge* e2) { return e1!=e2; }
// Set of candidate edges in stoc
// HashStruct so that hashing does not use pointer values (portable random)
static HashStruct<MEdge> ecand(hashe,cmpe);

// Declarations
static int domfilename(int argc, char const** argv);
static int dofilename(int argc, char const** argv);
static int dorecord(int argc, char const** argv);
static void doreconstruct();
static int dogfit(int argc, char const** argv);
static int dofgfit(int argc, char const** argv);
static void dointerp();
static void dostoc();
static int dooutmesh(int argc, char const** argv);
//
static int checksignal();
static void sigUSR1();
static void sigUSR2();
//
static void markmesh(GMesh& m);
static void subdivide(SubMesh& smesh);
// unordered

int main(int argc, char const** argv)
{
	SHOWDF("Created using:\n%s",CreateHeader(argc,argv));
	gbbox.clear();
	Options opts;
	OPTSDC(opts,mfilename,	"name : read initial mesh (can be -)");
	OPTSDC(opts,filename,	"name : read point data (can be -)");
	OPTSFC(opts,outn,	": at end, write meshn instead of mesh0");
	opts.c("",":");
	OPTSDC(opts,record,	"file : send orig+changes to file");
	OPTSPC(opts,crep,	"val : set repr. constant (def. 0)");
	OPTSPC(opts,csharp,	"val : penalize sharp edges (def. 0)");
	OPTSDC(opts,reconstruct,": run standard optimization");
	opts.c("",":");
	OPTSDC(opts,gfit,	"niter : run global optimization (fixed K)");
	OPTSDC(opts,fgfit,	"niter : conjugate gradient method");
	OPTSDC(opts,interp,	": find interpolating mesh");
	OPTSFC(opts,ecol,	":  try edge collapses, keep sharp topology");
	OPTSFC(opts,gecol,	":  try general edge collapses");
	OPTSFC(opts,esha,	":  try sharp changes");
	OPTSFC(opts,eswa,	":  try edge swaps");
	OPTSFC(opts,espl,	":  try edge splits");
	OPTSDC(opts,stoc,	": run local discrete optimization");
	opts.c("",":");
	OPTSDC(opts,outmesh,	"filename : output current mesh0 to file");
	OPTSFC(opts,nooutput,	": do not print final mesh0");
	opts.c("",":");
	OPTSPC(opts,feshaasym,	"f : set drss threshold (fraction of crep)");
	OPTSPC(opts,feswaasym,	"f : set drss threshold (fraction of crep)");
	opts.c("",":");
	OPTSPC(opts,nsubdiv,	"i : number of subdivision iters (def. 2)");
	OPTSFC(opts,nolimit,	": do not send final vertices to limit");
	opts.c("",":");
	OPTSPC(opts,selective,	"deg : refine sharpe edges and others>deg");
	OPTSFC(opts,s222,	": use s222 mask instead of n222");
	OPTSPC(opts,weighta,	"a : override interior extraord. weight");
	opts.c("",":");
	opts.c("    Examples"
	       ,  ": Subdivfit -mf - -outn");
	opts.c("",": Subdivfit -mf v.m -fi v.pts  -gfit 10  -ecol -stoc >v");
	opts.c("",":");
	signal(SIGUSR1,HHSIG_PF(sigUSR1));
	signal(SIGUSR2,HHSIG_PF(sigUSR2));
	TIMER(Subdivfit);
	opts.allmustparse();
	if (!opts.parse(argc,argv)) { opts.problem(argc,argv); return 1; }
	if (outn) {
		assertx(gmesh.numVertices());
		ForMeshVertex(gmesh,v) {
			gmesh.modFlag(v,SubMesh::VARIABLEV,0);
		} EndFor;
		SubMesh smesh(gmesh);
		subdivide(smesh);
		smesh.updatevertexpositions();
		ETIMER(Subdivfit);
		CleanUp();
		if (!nooutput) {
			GMesh& m=smesh.mesh();
			markmesh(m); m.write(cout);
		}
	} else {
		ETIMER(Subdivfit);
		CleanUp();
		if (!nooutput) { markmesh(gmesh); gmesh.write(cout); }
	}
	if (wfrecord) { gmesh.recordChanges(0); delete wfrecord,wfrecord=0; }
	gmesh.clear(); co.clear(); gbary.clear(); gclp.clear();
	return 0;
}

static int domfilename(int argc, char const** argv)
{
	assertx(argc>1);
	assertx(!gmesh.numVertices());
	RFile is(argv[1]);
	gmesh.read(is());
	SHOWDF("Initial mesh:\n");
	SHOWDF(MeshGenusString(gmesh));
	ForMeshVertex(gmesh,v) { gbbox.takeunion(gmesh.point(v)); } EndFor;
	return 2;
}

static int dofilename(int argc, char const** argv)
{
	assertx(argc>1);
	assertx(!co.num());
	RFile is(argv[1]);
	RSA3dStream ia3d(is());
	A3dElem el;
	for (;;) {
		ia3d.read(el);
		if (el.type()==A3dElem::TEndFile) break;
		if (el.type()==A3dElem::TComment) continue;
		if (assertw1(el.type()==A3dElem::TPoint)) continue;
		co+=el[0].p;
		gbbox.takeunion(el[0].p);
		gcmf+=Face(0); gscmfi+=0; gdis2+=double(0); gscmf+=Face(0);
		gbary+=Bary(0,0,0); gclp+=Point(0,0,0);
	}
	SHOWDF("%d points read\n",co.num());
	return 2;
}

static void initialize()
{
	static int done=0;
	if (!done++) {
		xform=gbbox.getFrameToSmallCube();
		SHOWDF(" internal xform: %s",FrameIO::string(xform,1,0));
		xformscale=assertv(xform[0][0]);
	}
	wcrep=crep/square(xformscale);
	wcsharp=csharp/square(xformscale);
}

//*** auxiliary

static void markmesh(GMesh& m)
{
	ForMeshVertex(m,v) {
		m.setInfo(v,m.flag(v,GMesh::CUSPV)?
			  new MeshInfoString("cusp"):0);
	} EndFor;
	ForMeshEdge(m,e) {
		m.setInfo(e,m.flag(e,GMesh::SHARPE)?
			  new MeshInfoString("sharp"):0);
	} EndFor;
}

static void wfframe()
{
	if (!wfrecord) return;
	(*wfrecord)() << "# frame\n" << flush;
}

static void subdivide(SubMesh& smesh)
{
	ATIMER(___submesh);
	smesh.maskparameters(s222,weighta);
	smesh.subdividen(nsubdiv,!nolimit,cos(torad(selective)));
}

// translate Vertex from one Mesh to another Mesh
static Vertex trvmm(Vertex v, const Mesh& mf, const Mesh& mt)
{
	return mt.idvertex(mf.vertexid(v));
}

// translate Face from one Mesh to another Mesh
static Face trfmm(Face f, const Mesh& mf, const Mesh& mt)
{
	return mt.idface(mf.faceid(f));
}

// translate Edge from one Mesh to another Mesh
static Edge tremm(Edge e, const Mesh& mf, const Mesh& mt)
{
	return mt.edge(trvmm(mf.vertex1(e),mf,mt),
		       trvmm(mf.vertex2(e),mf,mt));
}

static int vnumsharpe(const GMesh& mesh, Vertex v)
{
	int nsharpe=0;
	ForVertexEdge(mesh,v,e) {
		if (mesh.isBoundary(e) || mesh.flag(e,GMesh::SHARPE))
			nsharpe++;
	} EndFor;
	return nsharpe;
}

static double mindihedralaboutvertices(const GMesh& mesh, Vertex v)
{
	Set<Edge> sete;
	ForVertexVertex(mesh,v,vv) {
		ForVertexEdge(mesh,vv,e) {
			sete.add(e);
		} EndFor;
	} EndFor;
	double mindic=2;
	ForSet(sete,Edge,e) {
		if (mesh.isBoundary(e)) continue;
		double dic=EdgeDihedralAngleCos(mesh,e);
		mindic=min(mindic,dic);
	} EndFor;
	return mindic;
}

static int oksharpedgechange(Edge e)
{
	if (!gmesh.flag(e,GMesh::SHARPE)) {
		if (vnumsharpe(gmesh,gmesh.vertex1(e))>=1 &&
		    vnumsharpe(gmesh,gmesh.vertex2(e))>=1) return 0;
	} else {
		if (vnumsharpe(gmesh,gmesh.vertex1(e))>=3 &&
		    vnumsharpe(gmesh,gmesh.vertex2(e))>=3) return 0;
	}
	return 1;
}

static double getedis()
{
	double sum=0;
	ForIndex(i,co.num()) { sum+=gdis2[i]; } EndFor;
	return sum;
}

static int numsharpe()
{
	int nsharpe=0;
	// boundary edges are never sharp!
	ForMeshEdge(gmesh,e) {
		if (gmesh.flag(e,GMesh::SHARPE)) {
			assertx(!gmesh.isBoundary(e));
			nsharpe++;
		}
	} EndFor;
	return nsharpe;
}

static double getetot()
{
	return getedis()+gmesh.numVertices()*wcrep+numsharpe()*wcsharp;
}

static void analyzemesh(const char* s)
{
	SHOWDF("%s: v=%d nse=%d/%d  edis=%g etot=%g\n",
	       s,gmesh.numVertices(),numsharpe(),gmesh.numEdges(),
	       getedis(),getetot());
}

static void gallproject(const SubMesh& smesh)
{
	ATIMER(___gallproject);
	const GMesh& mesh=smesh.mesh();
	PolygonSpatial psp(35);
	Map<Polygon*,Face> mpf;
	NEST {
		ATIMER(____gmakeSpatial);
		ForMeshFace(mesh,f) {
			Polygon& spoly=*new Polygon(3);
			mesh.polygon(f,spoly);
			ForIndex(i,spoly.num()) { spoly[i]*=xform; } EndFor;
			psp.enter(&spoly);
			mpf.enter(&spoly,f);
		} EndFor;
	}
	NEST {
		ATIMER(____gspatialproject);
		ForIndex(i,co.num()) {
			SpatialSearch ss(psp,co[i]*xform);
			Polygon* spoly=Conv<Polygon*>::d(ss.next());
			gscmf[i]=mpf.get(spoly);
			static Polygon poly;
			mesh.polygon(gscmf[i],poly);
			ProjectPTri(co[i],poly[0],poly[1],poly[2],
				    &gdis2[i],0,&gbary[i],&gclp[i]);
		} EndFor;
	}
	ForMapKey(mpf,Polygon*,Face,poly) { delete poly; } EndFor;
}

static void gneighproject(const SubMesh& smesh)
{
	ATIMER(___gneighproject);
	ForIndex(i,co.num()) {
		ProjectPNeighb(smesh.mesh(),co[i],gscmf[i],
			       &gdis2[i],&gbary[i],&gclp[i],1);
	} EndFor;
}

static void globallls(SubMesh& smesh, double& rss0, double& rss1)
{
	ATIMER(___glls);
 	GMesh& omesh=smesh.origmesh();
 	GMesh& mesh=smesh.mesh();
	Map<Vertex,int> mvi; Array<Vertex> iv;
	ForMeshVertex(omesh,v) {
		mvi.enter(v,iv.num()); iv+=v;
	} EndFor;
	int m=co.num(), n=iv.num();
	SparseLLS* lls=new SparseLLS(m,n,3); // LLS::create(m,n,3,12./n);
	lls->maxiter(10);
	ForIndex(i,m) {
		Vertex va[3]; mesh.vertices(gscmf[i],va);
		Combvh tricomb, comb;
		for (int j=0;j<3;j++) tricomb.c[va[j]]=gbary[i][j];
		smesh.composecmvcvh(tricomb,comb);
		SSTAT(Scombnum,comb.c.num());
		Homogeneous h=co[i]-comb.h;
		lls->enterRh(i,&h[0]);
		ForCombination(comb.c,Vertex,v,val) {
			lls->enter(i,mvi.get(v),val);
		} EndFor;
	} EndFor;
	ForIndex(i,n) {
		lls->enterEstimate(i, static_cast<double const *>(omesh.point(iv[i])));
	} EndFor;
	NEST { ATIMER(____gsolve); lls->solve(&rss0,&rss1); }
	ForIndex(i,n) {
		Point p; lls->getValue(i,&p[0]);
		omesh.setPoint(iv[i],p);
	} EndFor;
	delete lls,lls=0;
}

//*** main procedures

static int dorecord(int argc, char const** argv)
{
	assertx(argc>1);
	wfrecord=new WFile(argv[1]);
	markmesh(gmesh);
	gmesh.write((*wfrecord)());
	gmesh.recordChanges(&(*wfrecord)());
	return 2;
}

static int dogfit(int argc, char const** argv)
{
	assertx(co.num() && gmesh.numVertices());
	initialize();
	ostream* os=gmesh.recordChanges(0);
	SHOWDF("\n");
	TIMER(_gfit);
	assertx(argc>1); int niter=atoi(argv[1]);
	ForMeshVertex(gmesh,v) {
		gmesh.modFlag(v,SubMesh::VARIABLEV,1);
	} EndFor;
	SubMesh smesh(gmesh);
	subdivide(smesh);
	smesh.updatevertexpositions();
	gallproject(smesh);
	analyzemesh("gfit_before");
	ForIndex(ni,niter) {
		if (checksignal()) break;
		ATIMER(__gfit_iter);
		double rss0, rss1;
		globallls(smesh,rss0,rss1);
		SHOWF(" gopt %d/%d lls rss0=%g rss1=%g\n",
		      ni+1,niter,rss0,rss1);
		smesh.updatevertexpositions();
		gneighproject(smesh);
	} EndFor;
	if (os) {
		gmesh.recordChanges(os);
		ForMeshVertex(gmesh,v) {
			gmesh.setPoint(v,gmesh.point(v));
		} EndFor
		wfframe();
	}
	analyzemesh("gfit_after ");
	return 2;
}

//*** fgfit

static double dot(const Array<Vector>& a, const Array<Vector>&b)
{
	assertx(a.num() && a.num()==b.num());
	register double sum=0;
	ForIndex(i,a.num()) { sum+=dot(a[i],b[i]); } EndFor;
	return sum;
}

class FG {
  public:
	FG(GMesh& pmesh, SubMesh& psmesh);
	~FG();
	void computelinedir();
	double computeonline(double t);
	double computedonline(); // assumes previous t!
  private:
	GMesh& mesh; SubMesh& smesh;
	int nv;
	int first;
	Map<Vertex,int> mvi; Array<Vertex> iv;
	Array<Point> opos;
	Array<Vector> dir, g, h;
	void computegrad(Array<Vector>& grad);
	DISABLECOPY(FG);	// GNUG 2.5.8
};

static FG* fg;

FG::FG(GMesh& pmesh, SubMesh& psmesh)
: mesh(pmesh), smesh(psmesh), nv(mesh.numVertices()), first(1),
  opos(nv), dir(nv), g(nv), h(nv)
{
	ForMeshVertex(mesh,v) { mvi.enter(v,iv.num()); iv+=v; } EndFor;
}

FG::~FG() { }

struct fgVInfo : MeshInfo {
	fgVInfo() : g(0,0,0) { }
	Vector g;
};

void FG::computegrad(Array<Vector>& grad)
{
	ATIMER(___computegrad);
	assertx(grad.num()==nv);
	ForMeshVertex(mesh,v) { mesh.setInfo(v,new fgVInfo); } EndFor;
	ForIndex(i,co.num()) {
		Vertex va[3]; smesh.mesh().vertices(gscmf[i],va);
		Vector vtop=co[i]-gclp[i];
		// this is faster than composecmvcv(tricomb,comb);
		ForIndex(j,3) {
			double a=-2*gbary[i][j];
			if (!a) continue;
			const Combvh& comb=smesh.combination(va[j]);
			assertx(comb.h.iszero());
			ForCombination(comb.c,Vertex,v,val) {
				reinterpret_cast<fgVInfo *>(mesh.info(v))->g += vtop*(a*val);
			} EndFor;
		} EndFor;
	} EndFor;
	ForIndex(j,nv) { grad[j] = reinterpret_cast<fgVInfo *>(mesh.info(iv[j]))->g; } EndFor;
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
	smesh.updatevertexpositions();
	gneighproject(smesh);
	double e=getedis();
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
	ostream* os=gmesh.recordChanges(0);
	assertx(co.num() && gmesh.numVertices());
	initialize();
	SHOWDF("\n");
	TIMER(_fgfit);
	assertx(argc>1); int niter=atoi(argv[1]);
	ForMeshVertex(gmesh,v) {
		gmesh.modFlag(v,SubMesh::VARIABLEV,1);
	} EndFor;
	SubMesh smesh(gmesh);
	subdivide(smesh);
	smesh.updatevertexpositions();
	gallproject(smesh);
	analyzemesh("fgfit_before");
	fg=new FG(gmesh,smesh);
	double oldtmin=.05;
	ForIndex(ni,niter) {
		if (checksignal()) break;
		ATIMER(__fgfit_iter);
		fg->computelinedir();
		const double ftol=.05;
		double t1=oldtmin*.9, t2=oldtmin*1.1+.001;
		double tmin, emin; int neval;
		LineOptimization(&fgeval,&fgevald,t1,t2,ftol,tmin,emin,neval);
		SHOWF("iter %d/%d  neval=%d tmin=%9.6lf e=%lg\n",
		      ni+1,niter,neval,tmin,emin);
		oldtmin=tmin;
		SSTAT(Sneval,neval); SSTAT(Stmin,tmin);
	} EndFor;
	delete fg,fg=0;
	if (os) {
		gmesh.recordChanges(os);
		ForMeshVertex(gmesh,v) {
			gmesh.setPoint(v,gmesh.point(v));
		} EndFor
		wfframe();
	}
	analyzemesh("fgfit_after ");
	return 2;
}

static void dointerp()
{
	assertx(gmesh.numVertices());
	initialize();
	SHOWDF("\n");
	TIMER(_interp);
	ForMeshVertex(gmesh,v) {
		gmesh.modFlag(v,SubMesh::VARIABLEV,1);
	} EndFor;
	SubMesh smesh(gmesh);
	subdivide(smesh);
	Map<Vertex,int> mvi; Array<Vertex> iv;
	ForMeshVertex(gmesh,v) {
		mvi.enter(v,iv.num()); iv+=v;
	} EndFor;
	int n=iv.num();
	LLS* lls=LLS::create(n,n,3,12./n);
	ForIndex(i,n) {
		Vertex v=iv[i];
		Point p=gmesh.point(v);
		lls->enterRh(i,&p[0]);
		Vertex vs=trvmm(v,gmesh,smesh.mesh());
		const Combvh& comb=smesh.combination(vs);
		assertx(comb.h.iszero());
		ForCombination(comb.c,Vertex,vv,val) {
			lls->enter(i,mvi.get(vv),val);
		} EndFor;
	} EndFor;
	ForIndex(i,n) {
		Point p=gmesh.point(iv[i]);
		lls->enterEstimate(i,&p[0]);
	} EndFor;
	double rss0, rss1;
	NEST { ATIMER(__isolve); lls->solve(&rss0,&rss1); }
	SHOWDF("lls rss0=%g, rss1=%g\n",rss0,rss1);
	ForIndex(i,n) {
		Point p; lls->getValue(i,&p[0]);
		gmesh.setPoint(iv[i],p);
	} EndFor;
	delete lls,lls=0;
	if (co.num()) {
		smesh.updatevertexpositions();
		gallproject(smesh);
		analyzemesh("interp_after ");
	}
}

//*** stoc

class Combvih {
  public:
	Combvih() {}
	Array<double> c;		// partial comb.: vertex index in iv -> double
	Homogeneous h;		// remainder of combination (constant term)
  private:
	DISABLECOPY(Combvih);	// GNUG 2.5.8
};

struct stocoVinfo : MeshInfo {
	int index;		// defined if vertex is in setmv
};

struct stocsVinfo : MeshInfo {
	stocsVinfo() {}
	Combvih c;
  private:
	DISABLECOPY(stocsVinfo); // GNUG 2.5.8
};

class Mvcvih {
  public:
	Mvcvih() {}
	Array<Vertex> iv;	// int -> mesh Vertex
  private:
	DISABLECOPY(Mvcvih);	// GNUG 2.5.8
};

static void createmvcvih(SubMesh& smesh, const Set<Vertex>& setmv,
			Mvcvih& mvcvih)
{
	GMesh& omesh=smesh.origmesh(); GMesh& m=smesh.mesh();
	ForMeshVertex(omesh,v) { omesh.setInfo(v,new stocoVinfo); } EndFor;
	ForMeshVertex(m,v) { m.setInfo(v,new stocsVinfo); } EndFor;
	ForSet(setmv,Vertex,v) {
		reinterpret_cast<stocoVinfo*>(m.info(v))->index=mvcvih.iv.num();
		mvcvih.iv+=v;
	} EndFor;
	int nv=setmv.num();
	ForMeshVertex(m,v) {
		Combvih& cvih=reinterpret_cast<stocsVinfo*>(m.info(v))->c;
		cvih.c.init(nv);
		ForIndex(j,nv) { cvih.c[j]=0; } EndFor;
		const Combvh& comb=smesh.combination(v);
		cvih.h=comb.h;
		ForCombination(comb.c,Vertex,vl,val) {
			int vi= reinterpret_cast<stocoVinfo*>(m.info(vl))->index;
			cvih.c[vi]=val;
		} EndFor;
	} EndFor;
}

static void updatelocal(SubMesh& smesh, const Mvcvih& mvcvih)
{
	GMesh& m=smesh.mesh();
	int nv=mvcvih.iv.num();
	ForMeshVertex(m,v) {
		const Combvih& cvih = reinterpret_cast<stocsVinfo*>(m.info(v))->c;
		Homogeneous h(cvih.h);
		ForIndex(j,nv) {
			double b=cvih.c[j];
			if (b) h+=smesh.origmesh().point(mvcvih.iv[j])*b;
		} EndFor;
		m.setPoint(v,toPoint(h));
	} EndFor;
}

static void localallproject(const SubMesh& smesh, const Set<Face>& setgoodf,
			    const Set<int>& setpts, const Set<int>& setbadpts)
{
	ATIMER(___lallproject);
	const GMesh& mesh=smesh.mesh();
	PolygonSpatial psp(35);
	Map<Polygon*,Face> mpf;
	NEST {
		ATIMER(____lmakespatial);
		ForMeshFace(mesh,f) {
			if (!setgoodf.contains(smesh.origface(f))) continue;
			Polygon& spoly=*new Polygon(3);
			mesh.polygon(f,spoly);
			ForIndex(i,spoly.num()) { spoly[i]*=xform; } EndFor;
			psp.enter(&spoly);
			mpf.enter(&spoly,f);
		} EndFor;
	}
	ATIMER(____lspatialproject);
	ForSet(setpts,int,i) {
		if (setbadpts.contains(i)) {
			SpatialSearch ss(psp,co[i]*xform);
			Polygon* spoly=Conv<Polygon*>::d(ss.next());
			gscmf[i]=mpf.get(spoly);
		} else {
			Face f=trfmm(gcmf[i],gmesh,smesh.origmesh());
			gscmf[i]=smesh.getface(f,gscmfi[i]);
		}
	} EndFor;
	ForMapKey(mpf,Polygon*,Face,poly) { delete poly; } EndFor;
}

// Optimization of setmv given setpts, gscmf, mvcvih
static void optimizelocal(SubMesh& smesh, const Set<Vertex>& setmv,
			  const Set<int>& setpts, const Mvcvih& mvcvih,
			  double& rss1)
{
	ATIMER(___loptimize);
	updatelocal(smesh,mvcvih);
	const GMesh& mesh=smesh.mesh();
	ATIMER(____lneighproject);
	ForSet(setpts,int,pi) {
		ProjectPNeighb(mesh,co[pi],gscmf[pi],0,&gbary[pi],0,1);
	} EndFor;
	ETIMER(____lneighproject);
	int m=setpts.num(), n=mvcvih.iv.num();
	assertx(n==setmv.num()); // optional
	LLS* lls=LLS::create(m,n,3,1.);
	int rowi=0;
	ATIMER(____lcombinations);
	ForSet(setpts,int,pi) {
		Vertex va[3]; mesh.vertices(gscmf[pi],va);
		const double* bary=&gbary[pi][0];
		register double b0=bary[0], b1=bary[1], b2=bary[2];
		const Combvih* cviha[3];
		const double* cvihcf[3]; const double* cvihhf[3];
		ForIndex(c,3) {
			cviha[c] = & reinterpret_cast<stocsVinfo*>(mesh.info(va[c]))->c;
		} EndFor;
		ForIndex(c,3) { cvihcf[c] = cviha[c]->c; } EndFor;
		ForIndex(c,3) { cvihhf[c] = static_cast<double const *>(cviha[c]->h); } EndFor;
		double h[3];
		h[0]=b0*cvihhf[0][0]+b1*cvihhf[1][0]+b2*cvihhf[2][0];
		h[1]=b0*cvihhf[0][1]+b1*cvihhf[1][1]+b2*cvihhf[2][1];
		h[2]=b0*cvihhf[0][2]+b1*cvihhf[1][2]+b2*cvihhf[2][2];
		ForIndex(j,n) {
			lls->enter(rowi,j,(b0*cvihcf[0][j]+b1*cvihcf[1][j]+
					   b2*cvihcf[2][j]));
		} EndFor;
		// Homogeneous hrh=co[pi]-h;
		Vector vrh; ForIndex(c,3) { vrh[c]=co[pi][c]-h[c]; } EndFor;
		lls->enterRh(rowi,&vrh[0]);
		rowi++;
	} EndFor;
	ETIMER(____lcombinations);
 	GMesh& omesh=smesh.origmesh();
	ForIndex(i,n) {
		const Point& p=omesh.point(mvcvih.iv[i]);
		lls->enterEstimate(i, static_cast<const double*>(p));
	} EndFor;
	NEST { ATIMER(____lsolve); lls->solve(0,&rss1); }
	if (0) {
		SHOWF("local optimization m=%d n=%d lls rss1=%g\n",
		      m,n,rss1);
	}
	ForIndex(i,n) {
		Point p; lls->getValue(i,&p[0]);
		omesh.setPoint(mvcvih.iv[i],p);
	} EndFor;
	delete lls;
}

static void buildlmesh1(const GMesh& gmesh, const Set<Vertex>& setgmv,
			const Set<Face>& setbadfg, GMesh& lmesh,
			Set<int>& setpts, Set<int>& setbadpts,
			double& rssf)
{
	Set<Vertex> setvg; Set<Face> setfg; Set<Face> setmfg;
	ForSet(setgmv,Vertex,v) {
		ForVertexVertex(gmesh,v,vv) {
			ForVertexVertex(gmesh,vv,vvv) {
				ForVertexVertex(gmesh,vvv,vvvv) {
					setvg.add(vvvv);
				} EndFor;
				ForVertexFace(gmesh,vvv,f) {
					setfg.add(f);
				} EndFor;
			} EndFor;
			ForVertexFace(gmesh,vv,f) {
				setmfg.add(f);
			} EndFor;
		} EndFor;
	} EndFor;
	ForSet(setvg,Vertex,v) {
		Vertex vn=lmesh.createVertexI(gmesh.vertexid(v));
		lmesh.setPoint(vn,gmesh.point(v));
	} EndFor;
	ForSet(setfg,Face,f) {
		Vertex va[3]; gmesh.vertices(f,va);
		for (int i=0;i<3;i++) va[i]=trvmm(va[i],gmesh,lmesh);
		assertx(lmesh.createFaceI(gmesh.faceid(f),va,3));
	} EndFor;
	// Add flags to lmesh
	ForMeshEdge(lmesh,e) {
		lmesh.modFlag(e,gmesh.flag(tremm(e,lmesh,gmesh),GMesh::ALL),1);
	} EndFor;
	ForMeshVertex(lmesh,v) {
		lmesh.modFlag(v,gmesh.flag(trvmm(v,lmesh,gmesh),GMesh::ALL),1);
	} EndFor;
	// Gather into setpts the points projecting on faces in setmfg
	ForSet(setmfg,Face,f) {
		ForSet(*mfpts.get(f),int,pi) { setpts.enter(pi); } EndFor;
	} EndFor;
	NEST {
		double sum=0;
		ForSet(setpts,int,pi) { sum+=gdis2[pi]; } EndFor;
		rssf=sum;
	}
	// Gather into setbadpts the points projecting on faces in setbadfg
	ForSet(setbadfg,Face,f) {
		ForSet(*mfpts.get(f),int,pi) { setbadpts.enter(pi); } EndFor;
	} EndFor;
}

static void buildlmesh2(GMesh& lmesh, const Set<Vertex>& setmv,
			Set<Face>& setmf)
{
	ForMeshVertex(lmesh,v) {
		lmesh.modFlag(v,SubMesh::VARIABLEV,setmv.contains(v));
	} EndFor;
	ForSet(setmv,Vertex,v) {
		ForVertexVertex(lmesh,v,vv) {
			ForVertexFace(lmesh,vv,f) {
				setmf.add(f);
			} EndFor;
		} EndFor;
	} EndFor;
}

// subdivide submesh; trim its outlying faces and vertices.
static void subdivtrim(SubMesh& smesh, const Set<int>& setpts,
		       const Set<Face>& setmf, const Set<Vertex>& setmv)
{
	subdivide(smesh);
	NEST {			// trim faces from smesh outside setmf
		Set<Face> setfrem;
		ForMeshFace(smesh.mesh(),f) {
			if (!setmf.contains(smesh.origface(f)))
				setfrem.enter(f);
		} EndFor;
		SSTAT(Ssetfrem,setfrem.num());
		ForSet(setfrem,Face,f) {
			smesh.mesh().destroyFace(f);
		} EndFor;
	}
	NEST {			// trim hanging vertices
		Set<Vertex> setvrem;
		ForMeshVertex(smesh.mesh(),v) {
			if (!smesh.mesh().degree(v)) setvrem.enter(v);
		} EndFor;
		SSTAT(Ssetvrem,setvrem.num());
		ForSet(setvrem,Vertex,v) {
			smesh.mesh().destroyVertex(v);
		} EndFor;
	}
}

static void lupdategmesh(const SubMesh& smesh, const Set<Vertex>& setmv,
			 const Set<int>& setpts,
			 Face f1l=0, Face f1g=0, Face f2l=0, Face f2g=0,
			 Vertex v1l=0, Vertex v1g=0)
{
	const GMesh& lmesh=smesh.origmesh();
	// Update positions of setmv in gmesh
	ForSet(setmv,Vertex,v) {
		Vertex vg= v==v1l?v1g: trvmm(v,lmesh,gmesh);
		gmesh.setPoint(vg,lmesh.point(v));
	} EndFor;
	ForSet(setpts,int,pi) {
		const GMesh& mesh=smesh.mesh();
		Vertex va[3]; mesh.vertices(gscmf[pi],va);
		ProjectPTri(co[pi],mesh.point(va[0]),mesh.point(va[1]),
			    mesh.point(va[2]),&gdis2[pi],0,0,0);
		Face fl; int index;
		smesh.origfaceindex(gscmf[pi],fl,index);
		Face fg= fl==f1l?f1g: fl==f2l?f2g: trfmm(fl,lmesh,gmesh);
		gscmfi[pi]=index;
		if (fg==gcmf[pi]) continue;
		assertx(mfpts.get(gcmf[pi])->remove(pi));
		gcmf[pi]=fg;
		mfpts.get(fg)->enter(pi);
	} EndFor;
}

static int tryopt(SubMesh& smesh, const Set<Vertex>& setmv,
		  const Set<int>& setpts, const Mvcvih& mvcvih,
		  double threshrss, double& edrss)
{
        int ni;
	double lrss=1e30, rss1;
	const double anticipation=3; const int maxtni=10;
	const int minni=6, maxni=12;
	// It is worthwhile investing into the fit since it may eliminate
	// future operations.
//	for (int ni=0;ni<maxtni;) {
	for (ni=0;ni<maxtni;) {
		optimizelocal(smesh,setmv,setpts,mvcvih,rss1);
		ni++;
		if (checksignal()) break;
		if (rss1-(lrss-rss1)*anticipation>threshrss) break;
		lrss=rss1;
	}
	SSTAT(Soptnit,ni);
	edrss=rss1-threshrss;
	if (edrss>=0) return 0;
	// do some more fitting
	for (;ni<maxni;) {
		optimizelocal(smesh,setmv,setpts,mvcvih,rss1);
		ni++;
		if (checksignal()) break;
		if (ni<minni) continue;
		if (lrss-rss1<wcrep*.1) break;
		lrss=rss1;
	}
	SSTAT(Soptnig,ni);
	edrss=rss1-threshrss;
	if (edrss>=0) { Warning("tryopt strange"); return 0; }
	updatelocal(smesh,mvcvih);
	return 1;
}

static void removeface(Face f)
{
	assertx(mfpts.contains(f));
	assertx(mfpts.get(f)->empty());
	delete mfpts.remove(f);
}

static OResult tryecol(Edge eg, double& edrss)
{
	// SHOWS("tryecol");
	if (!gmesh.niceEdgeCollapse(eg)) return OR_illegal;
	if (!gecol && !oksharpedgechange(eg)) return OR_sharp;
	ATIMER(__tryecol);
	GMesh lmesh; Set<int> setpts, setbadpts; double rssf;
	NEST {
		Set<Vertex> setgmv; Set<Face> setbadfg;
		ForEdgeVertex(gmesh,eg,v) {
			ForVertexVertex(gmesh,v,vv) { setgmv.add(vv); } EndFor;
			ForVertexFace(gmesh,v,f) { setbadfg.add(f); } EndFor;
		} EndFor;
		buildlmesh1(gmesh,setgmv,setbadfg,lmesh,setpts,setbadpts,rssf);
	}
	double minb=1e30;
	ForEdgeVertex(gmesh,eg,v) {
		minb=min(minb,mindihedralaboutvertices(gmesh,v));
	} EndFor;
	Edge e=tremm(eg,gmesh,lmesh);
	Vertex v1=lmesh.vertex1(e), v2=lmesh.vertex2(e);
	Point p1=lmesh.point(v1), p2=lmesh.point(v2);
	int ndsharp=lmesh.flag(e,GMesh::SHARPE)?1:0;
	NEST {
		Vertex vo1=lmesh.sideVertex1(e), vo2=lmesh.sideVertex2(e);
		if (lmesh.flag(lmesh.edge(v1,vo1),GMesh::SHARPE) &&
		    lmesh.flag(lmesh.edge(v2,vo1),GMesh::SHARPE)) ndsharp++;
		if (vo2 &&
		    lmesh.flag(lmesh.edge(v1,vo2),GMesh::SHARPE) &&
		    lmesh.flag(lmesh.edge(v2,vo2),GMesh::SHARPE)) ndsharp++;
	}
	lmesh.collapseEdge(e);	// keep v1
	Set<Vertex> setmv; setmv.enter(v1);
	ForVertexVertex(lmesh,v1,vv) { setmv.enter(vv); } EndFor;
	Set<Face> setmf; buildlmesh2(lmesh,setmv,setmf);
	Set<Face> setgoodf;
	ForVertexFace(lmesh,v1,f) { setgoodf.enter(f); } EndFor;
	SubMesh smesh(lmesh); subdivtrim(smesh,setpts,setmf,setmv);
	Mvcvih mvcvih; createmvcvih(smesh,setmv,mvcvih);
	SSTAT(Secolpts,setpts.num());
	SSTAT(Secolmf,setmf.num()); SSTAT(Secolmv,setmv.num());
	SSTAT(Secolsmv,smesh.mesh().numVertices());
	lmesh.setPoint(v1,interp(p1,p2));
	updatelocal(smesh,mvcvih);
	localallproject(smesh,setgoodf,setpts,setbadpts); // get guess
	double minrss=1e30; int minii=-1;
	ForIndex(ii,3) {
		lmesh.setPoint(v1,interp(p1,p2,ii*.5));
		double mina=mindihedralaboutvertices(lmesh,v1);
		if (mina<mincos && mina<minb) continue;	// disallow
		double rss1;
		optimizelocal(smesh,setmv,setpts,mvcvih,rss1);
		mina=mindihedralaboutvertices(lmesh,v1);
		// Return to initial state.
		ForSet(setmv,Vertex,v) {
			lmesh.setPoint(v,gmesh.point(trvmm(v,lmesh,gmesh)));
		} EndFor;
		if (mina<mincos && mina<minb) continue;	// disallow
		if (rss1<minrss) minrss=rss1,minii=ii;
	} EndFor;
	if (minii<0) return OR_dih;
	double threshrss=rssf+wcrep+ndsharp*wcsharp;
	// SHOWF("minii=%d threshrss=%g\n",minii,threshrss);
	lmesh.setPoint(v1,interp(p1,p2,minii*.5));
	if (!tryopt(smesh,setmv,setpts,mvcvih,threshrss,edrss))
		return OR_energy;
	double mina=mindihedralaboutvertices(lmesh,v1);
	if (mina<mincos && mina<minb) return OR_dih;
	// ALL SYSTEMS GO
	Face f1g=gmesh.face1(eg), f2g=gmesh.face2(eg);
	ForEdgeVertex(gmesh,eg,v) {
		ForVertexEdge(gmesh,v,e) { ecand.remove(e); } EndFor;
	} EndFor;
	Vertex v1g=gmesh.vertex1(eg);
	gmesh.collapseEdge(eg);	// keep v1g
	// add about 12-16 edges
	ForVertexEdge(gmesh,v1g,e) { ecand.add(e); } EndFor;
	ForVertexFace(gmesh,v1g,f) { ecand.add(gmesh.oppEdge(v1g,f)); } EndFor;
	lupdategmesh(smesh,setmv,setpts);
	removeface(f1g); if (f2g) removeface(f2g);
	return OR_success;
}

static OResult tryesha(Edge eg, double& edrss)
{
	// SHOWS("tryesha");
//	const testdih=0;
	const int testdih=0;
	if (gmesh.isBoundary(eg)) return OR_illegal;
	int issharp=gmesh.flag(eg,GMesh::SHARPE);
	double vcos=EdgeDihedralAngleCos(gmesh,eg);
	static const double cos30d=cos(torad(30));
	if (!issharp && vcos>cos30d) return OR_sharp; // quick culling
	// if issharp then always consider smoothing it
	ATIMER(__tryesha);
	// setbadpts and setgoodf are empty
	GMesh lmesh; Set<int> setpts, setbadpts; double rssf;
	NEST {
		Set<Vertex> setgmv; Set<Face> setbadfg;
		ForEdgeVertex(gmesh,eg,v) {
			ForVertexVertex(gmesh,v,vv) { setgmv.add(vv); } EndFor;
		} EndFor;
		buildlmesh1(gmesh,setgmv,setbadfg,lmesh,setpts,setbadpts,rssf);
	}
	double minb=1e30;
	if (testdih) {
		ForEdgeVertex(gmesh,eg,v) {
			minb=min(minb,mindihedralaboutvertices(gmesh,v));
		} EndFor;
	}
	Edge e=tremm(eg,gmesh,lmesh);
	lmesh.modFlag(e,GMesh::SHARPE,!issharp);
	Set<Vertex> setmv;
	ForEdgeVertex(lmesh,e,v) {
		ForVertexVertex(lmesh,v,vv) { setmv.add(vv); } EndFor;
	} EndFor;
	Set<Face> setmf; buildlmesh2(lmesh,setmv,setmf);
	Set<Face> setgoodf;
	SubMesh smesh(lmesh); subdivtrim(smesh,setpts,setmf,setmv);
	Mvcvih mvcvih; createmvcvih(smesh,setmv,mvcvih);
	SSTAT(Seshapts,setpts.num());
	SSTAT(Seshamf,setmf.num()); SSTAT(Seshamv,setmv.num());
	SSTAT(Seshasmv,smesh.mesh().numVertices());
	updatelocal(smesh,mvcvih);
	localallproject(smesh,setgoodf,setpts,setbadpts);
	double threshrss=rssf-wcrep*feshaasym+(issharp?1:-1)*wcsharp;
	if (!tryopt(smesh,setmv,setpts,mvcvih,threshrss,edrss))
		return OR_energy;
	if (testdih) {
		double mina=1e30;
		ForEdgeVertex(lmesh,e,v) {
			mina=min(mina,mindihedralaboutvertices(lmesh,v));
		} EndFor;
		if (mina<mincos && mina<minb) return OR_dih;
	}
	// ALL SYSTEMS GO
	issharp=!issharp;
	gmesh.modFlag(eg,GMesh::SHARPE,issharp);
	if (wfrecord) {
		(*wfrecord)() << "Edge " << gmesh.vertexid(gmesh.vertex1(e)) <<
			" " << gmesh.vertexid(gmesh.vertex2(e)) << " {" <<
			(issharp?"sharp":"") << "}\n";
	}
	ForEdgeVertex(gmesh,eg,v) {
		ForVertexEdge(gmesh,v,e) { ecand.add(e); } EndFor;
		ForVertexFace(gmesh,v,f) {
			ecand.add(gmesh.oppEdge(v,f));
		} EndFor;
	} EndFor;
	lupdategmesh(smesh,setmv,setpts);
	return OR_success;
}

static OResult tryeswa(Edge eg, double& edrss)
{
	// SHOWS("tryeswa");
	if (!gmesh.legalEdgeSwap(eg)) return OR_illegal;
	int issharp=gmesh.flag(eg,GMesh::SHARPE);
	Vertex v1g=gmesh.vertex1(eg), v2g=gmesh.vertex2(eg);
	Vertex vo1g=gmesh.sideVertex1(eg), vo2g=gmesh.sideVertex2(eg);
	Face of1g=gmesh.face1(eg), of2g=gmesh.face2(eg);
	double minb=EdgeDihedralAngleCos(gmesh,eg);
	double mina=DihedralAngleCos(gmesh.point(vo1g),gmesh.point(vo2g),
				    gmesh.point(v1g),gmesh.point(v2g));
	if (mina<mincos && mina<minb) return OR_dih;
	// could do culling check if mina>cos5 && minb>cos5 ??
	ATIMER(__tryeswa);
	GMesh lmesh; Set<int> setpts, setbadpts; double rssf;
	NEST {
		Set<Vertex> setgmv; Set<Face> setbadfg;
		setgmv.enter(v1g); setgmv.enter(v2g);
		setgmv.enter(vo1g); setgmv.enter(vo2g);
		setbadfg.enter(of1g); setbadfg.enter(of2g);
		buildlmesh1(gmesh,setgmv,setbadfg,lmesh,setpts,setbadpts,rssf);
	}
	Edge oe=tremm(eg,gmesh,lmesh);
	Vertex v1=lmesh.vertex1(oe), v2=lmesh.vertex2(oe);
	Vertex vo1=lmesh.sideVertex1(oe), vo2=lmesh.sideVertex2(oe);
	Edge ne=lmesh.swapEdge(oe); // new edge is not sharp.  make sharp??
	Face nf1l=lmesh.face1(ne), nf2l=lmesh.face2(ne);
	Set<Vertex> setmv;
	setmv.enter(v1); setmv.enter(v2); setmv.enter(vo1); setmv.enter(vo2);
	Set<Face> setmf; buildlmesh2(lmesh,setmv,setmf);
	Set<Face> setgoodf; setgoodf.enter(nf1l); setgoodf.enter(nf2l);
	SubMesh smesh(lmesh); subdivtrim(smesh,setpts,setmf,setmv);
	Mvcvih mvcvih; createmvcvih(smesh,setmv,mvcvih);
	SSTAT(Seswapts,setpts.num());
	SSTAT(Seswamf,setmf.num()); SSTAT(Seswamv,setmv.num());
	SSTAT(Seswasmv,smesh.mesh().numVertices());
	updatelocal(smesh,mvcvih);
	localallproject(smesh,setgoodf,setpts,setbadpts);
	double threshrss=rssf-wcrep*feswaasym+(issharp?1:0)*wcsharp;
	if (!tryopt(smesh,setmv,setpts,mvcvih,threshrss,edrss))
		return OR_energy;
	// ALL SYSTEMS GO
	ecand.remove(eg);
	ecand.remove(gmesh.edge(v1g,vo1g)); ecand.remove(gmesh.edge(v2g,vo1g));
	ecand.remove(gmesh.edge(v1g,vo2g)); ecand.remove(gmesh.edge(v2g,vo2g));
	Edge neg=gmesh.swapEdge(eg); // new edge is not sharp!!
	Face nf1g=gmesh.face1(neg), nf2g=gmesh.face2(neg);
	ecand.add(neg);
	ecand.add(gmesh.edge(v1g,vo1g)); ecand.add(gmesh.edge(v2g,vo1g));
	ecand.add(gmesh.edge(v1g,vo2g)); ecand.add(gmesh.edge(v2g,vo2g));
	if (nf1g!=of1g && nf1g!=of2g) mfpts.enter(nf1g,new Set<int>);
	if (nf2g!=of1g && nf2g!=of2g) mfpts.enter(nf2g,new Set<int>);
	lupdategmesh(smesh,setmv,setpts,nf1l,nf1g,nf2l,nf2g);
	if (of1g!=nf1g && of1g!=nf2g) removeface(of1g);
	if (of2g!=nf1g && of2g!=nf2g) removeface(of2g);
	return OR_success;
}

static OResult tryespl(Edge eg, double& edrss)
{
	// SHOWS("tryespl");
	// always legal
	int issharp=gmesh.flag(eg,GMesh::SHARPE);
	Vertex v1g=gmesh.vertex1(eg), v2g=gmesh.vertex2(eg);
	Vertex vo1g=gmesh.sideVertex1(eg), vo2g=gmesh.sideVertex2(eg);
	Face f1g=gmesh.face1(eg), f2g=gmesh.face2(eg);
	// vo2g and f2g may be zero
	ATIMER(__tryespl);
	GMesh lmesh; Set<int> setpts, setbadpts; double rssf;
	NEST {
		Set<Vertex> setgmv; Set<Face> setbadfg;
		setgmv.enter(v1g); setgmv.enter(v2g);
		setgmv.enter(vo1g); if (vo2g) setgmv.enter(vo2g);
		setbadfg.enter(f1g); if (f2g) setbadfg.enter(f2g);
		buildlmesh1(gmesh,setgmv,setbadfg,lmesh,setpts,setbadpts,rssf);
	}
	Edge e=tremm(eg,gmesh,lmesh); Vertex v2=lmesh.vertex2(e);
	Vertex vn=lmesh.splitEdge(e);
	Edge eul=lmesh.edge(vn,v2);
	Face nf1l=lmesh.face1(eul), nf2l=lmesh.face2(eul);
	// edges (vn,vo1) [and (vn,vo2)] are not sharp
	Set<Vertex> setmv; setmv.enter(vn);
	ForVertexVertex(lmesh,vn,v) { setmv.enter(v); } EndFor;
	Set<Face> setmf; buildlmesh2(lmesh,setmv,setmf);
	Set<Face> setgoodf;
	ForVertexFace(lmesh,vn,f) { setgoodf.enter(f); } EndFor;
	SubMesh smesh(lmesh); subdivtrim(smesh,setpts,setmf,setmv);
	Mvcvih mvcvih; createmvcvih(smesh,setmv,mvcvih);
	SSTAT(Sesplpts,setpts.num());
	SSTAT(Sesplmf,setmf.num()); SSTAT(Sesplmv,setmv.num());
	SSTAT(Sesplsmv,smesh.mesh().numVertices());
	updatelocal(smesh,mvcvih);
	localallproject(smesh,setgoodf,setpts,setbadpts);
	double threshrss=rssf-wcrep-(issharp?1:0)*wcsharp;
	if (!tryopt(smesh,setmv,setpts,mvcvih,threshrss,edrss))
		return OR_energy;
	// ALL SYSTEMS GO
	ForEdgeFace(gmesh,eg,f) {
		ForFaceEdge(gmesh,f,e) { ecand.remove(e); } EndFor;
	} EndFor;
	Vertex vng=gmesh.splitEdge(eg);
	ForVertexFace(gmesh,vng,f) {
		ForFaceEdge(gmesh,f,e) { ecand.add(e); } EndFor;
	} EndFor;
	Edge eug=gmesh.edge(vng,v2g);
	Face nf1g=gmesh.face1(eug), nf2g=gmesh.face2(eug);
	mfpts.enter(nf1g,new Set<int>);
	if (nf2g) mfpts.enter(nf2g,new Set<int>);
	lupdategmesh(smesh,setmv,setpts,nf1l,nf1g,nf2l,nf2g,vn,vng);
	return OR_success;
}

static OResult tryop(Edge e, Operation op, double& edrss)
{
	if (checksignal()) return OR_illegal;
	OResult oor;
	oor=(op==OP_ecol?tryecol(e,edrss):
	     op==OP_esha?tryesha(e,edrss):
	     op==OP_eswa?tryeswa(e,edrss):
	     op==OP_espl?tryespl(e,edrss):
	     (assertnever(""),OR_success));
	opstat[op][oor]++;
	return oor;
}

static void stocinit()
{
	DTIMER(_initial_fit);
	ForMeshFace(gmesh,f) { mfpts.enter(f,new Set<int>); } EndFor;
	NEST {
		ForMeshVertex(gmesh,v) {
			gmesh.modFlag(v,SubMesh::VARIABLEV,0);
		} EndFor;
		SubMesh smesh(gmesh);
		subdivide(smesh);
		smesh.updatevertexpositions();
		gallproject(smesh);
		// Build up global projection information
		const GMesh& mesh=smesh.mesh();
		ForIndex(i,co.num()) {
			smesh.origfaceindex(gscmf[i],gcmf[i],gscmfi[i]);
			mfpts.get(gcmf[i])->enter(i);
			Vertex va[3]; mesh.vertices(gscmf[i],va);
			ProjectPTri(co[i],mesh.point(va[0]),mesh.point(va[1]),
				    mesh.point(va[2]),&gdis2[i],0,0,0);
		} EndFor;
	}
	ForIndex(i,OPNUM) ForIndex(j,ORNUM) opstat[i][j]=0; EndFor; EndFor;
}

static void stocend()
{
	SHOWDF("Summary of attempts and results:\n");
	NEST {
		const char* s=hform("%20s","");
		ForIndex(i,OPNUM) { s=hform("%s%10s",s,opname[i]); } EndFor;
		SHOWDF("%s\n",s);
	}
	NEST {
		const char* s=hform("%20s"," total_attempts");
		ForIndex(i,OPNUM) {
			int sum=0;
			ForIndex(j,ORNUM) { sum+=opstat[i][j]; } EndFor;
			s=hform("%s%10d",s,sum);
		} EndFor;
		SHOWDF("%s\n",s);
	}
	ForIndex(j,ORNUM) {
		const char* s=hform("%20s",orname[j]);
		ForIndex(i,OPNUM) { s=hform("%s%10d",s,opstat[i][j]); } EndFor;
		SHOWDF("%s\n",s);
	} EndFor;
	NEST {
		ForMeshVertex(gmesh,v) {
			gmesh.modFlag(v,SubMesh::VARIABLEV,0);
		} EndFor;
		SubMesh smesh(gmesh);
		subdivide(smesh);
		smesh.updatevertexpositions();
		gallproject(smesh);
	}
	ForMeshFace(gmesh,f) { delete assertv(mfpts.remove(f)); } EndFor;
	assertx(!mfpts.num());
}

static void dostoc()
{
	SHOWDF("\n");
	TIMER(_stoc);
	assertx(co.num() && gmesh.numVertices());
	initialize();
	SHOWDF("Stoc, crep=%g csharp=%g wcrep=%g wcsharp=%g\n",
	       crep,csharp,wcrep,wcsharp);
	stocinit();
	analyzemesh("stoc_before");
	if (gecol) ecol=1;
	ForMeshEdge(gmesh,e) { ecand.enter(e); } EndFor;
	double cedis=getedis(), cetot=getetot();
	int i=0, nbad=0;
	for (;!ecand.empty();) {
		if (checksignal()) { ecand.clear(); break; }
		cout << flush;
		ATIMER(__lattempt);
		i++;
		HashStructIter<MEdge> hi(ecand,Random::G);
		Edge e=hi();
		assertx(ecand.remove(e));
		gmesh.valid(e); // optional
		Operation op=OP_ecol; OResult oor=OR_illegal; double edrss=0;
		if (oor!=OR_success && ecol) op=OP_ecol,oor=tryop(e,op,edrss);
		if (oor!=OR_success && esha) op=OP_esha,oor=tryop(e,op,edrss);
		if (oor!=OR_success && eswa) op=OP_eswa,oor=tryop(e,op,edrss);
		if (oor!=OR_success && espl) op=OP_espl,oor=tryop(e,op,edrss);
		if (oor==OR_success) wfframe();
		SHOWF("# it %5d, %s (after %3d) [%5d/%-5d] %s\n",
		      i,opname[op],nbad,ecand.num(),gmesh.numEdges(),
		      (oor==OR_success?hform("* success e=%e",edrss):
		       oor==OR_energy?hform("positive e=%e",edrss):
		       orname[oor]));
		if (oor==OR_success) nbad=0; else nbad++;
		if (oor==OR_success) {
			double nedis=getedis(), netot=getetot();
			// SHOWF("edis:%g->%g  etot:%g->%g\n",
			// cedis,nedis,cetot,netot);
			SSTAT(Sechange,netot-cetot);
			if (assertw(netot<=cetot)) {
				SSTAT(HHH_PECHANGE,netot-cetot);
			}
			cedis=nedis,cetot=netot;
		}
	}
	SHOWDF("it %d, last search: %d wasted attempts\n",i,nbad);
	SHOWDF("New mesh:\n");
	SHOWDF(MeshGenusString(gmesh));
	stocend();
	SHOWDF(" last cedis=%g cetot=%g\n",cedis,cetot);
	analyzemesh("stoc_after ");
	ecol=gecol=esha=eswa=espl=0; // reset flags (for perhaps other stoc)
}

static void doreconstruct()
{
	SHOWDF("Starting reconstruction sequence\n");
	SHOWDF(" crep=%g, csharp=%g\n",crep,csharp);
	initialize();
	const char* argv[10];
	argv[1]="40"; dofgfit(2, argv); // -fgfit 40
	ecol=1; dostoc();		       // -ecol -stoc
	argv[1]="10"; dofgfit(2, argv); // -fgfit 10
	esha=1; dostoc();		       // -esha -stoc
	argv[1]="10"; dofgfit(2, argv); // -fgfit 10
	gecol=1; esha=1; eswa=1; dostoc();     // -gecol -esha -eswa -stoc
	argv[1]="10"; dofgfit(2, argv); // -fgfit 10
}

static int dooutmesh(int argc, char const** argv)
{
	assertx(argc>1);
	WFile os(argv[1]);
	markmesh(gmesh);
	gmesh.write(os());
	return 2;
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

static int checksignal()
{
	if (gotsignalUSR2) { SHOWDF("*** sigUSR2 exit\n"); return 1; }
	if (!gotsignalUSR1) return 0;
	gotsignalUSR1=0;
	SHOWS("***starting outmesh");
	NEST {
		WFile os("v.subdivfit.m.Z");
		markmesh(gmesh);
		gmesh.write(os());
	}
	SHOWS("***done with outmesh");
	return 0;
}
