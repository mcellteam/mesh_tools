// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include <cstdarg>
#include "Options.h"
#include "GMesh.h"
#include "MeshOp.h"
#include "A3dStream.h"
#include "FileIO.h"
#include "HashPoint.h"
#include "Homogeneous.h"
#include "GeomOp.h"		// for DihedralAngleCos()
#include "Polygon.h"
#include "Array.h"
#include "Set.h"
#include "Map.h"
#include "Stack.h"
#include "Pqueue.h"
#include "Stat.h"
#include "Random.h"
#include "Timer.h"
#include "Principal.h"		// for PrincipalComponents()

#include <iostream>

using std::cin;
using std::cout;

static WSA3dStream oa3d(cout);
static GMesh mesh;
const A3dVertexColor SURFCOL(A3dColor(.8,.5,.4),A3dColor(.5,.5,.5),
			     A3dColor(5,0,0));
const double UNDEF=-1e31;

static double cosangle=UNDEF;
static double solidangle=0;
static int nooutput=0;
static int sloannormal=0;

// reduce options;
static int nfaces=0;
static double maxcrit=1e20;
static enum { Length=1, Inscribed=2} reducecrit;

// Declarations
static void dofroma3d();
static int doangle(int argc, char const** argv);
static int docosangle(int argc, char const** argv);
static int dosolidangle(int argc, char const** argv);
static int doexpsubdiv(int argc, char const** argv);
static void dosetb3d();
static void dotoa3d();
static void dotob3d();
static void doendobject();
static void dorenumber();
static void dooutmesh();
static void doaddmesh();
static void dorecord();
static void domark();
static void dodelaunay();
static void dodiagonal();
static void dosegment();
static void doselsmooth();
static void doshowsharpe();
static void dotrisubdiv();
static void dosilsubdiv();
static int dofillholes(int argc, char const** argv);
static void dotriangulate();
static int dormcomp(int argc, char const** argv);
static int docoalesce(int argc, char const** argv);
static void donice();
static void doflip();
static void dofixvertices();
static void dofixfaces();
static void dosmootha3d();
static void dogenus();
static void doremoveinfo();
static void dostat();
static void doreduce();
static void lengthc();
static void inscribedc();
static int dorandpts(int argc, char const** argv);
static int mybsearch(const double* ar, int l, int h, double v);
static Point trianglerandompoint(const Polygon& poly);
static void outputpoint(const Point& p);
static void dovertexpts();
static int dobndpts(int argc, char const** argv);
// auxiliary
static int sharp(Edge e);
static void gathersegments(Map<Face,int>& mfseg, Stack<Face>& stackrepf);
static void gatherfollowseg(Face f, Map<Face,int>& mfseg,
			    int segnum, int& nf, Point& pcentroid);
static void edgeanglestats(const Map<Face,int>& mfseg);
static void printsegment(Face f, Map<Face,int>& mfseg);
static void selsmooth();
static void meshsubdiv();
static int removeboundary(Queue<Edge>& queuee);
static void removesimpleboundary(const Queue<Edge>& queuee);
static int removecomponent(const Set<Face>& setf);
static int fixvertex(Vertex v);
static double reducecriterion(Edge e);
static double trycoalesce(Edge e);
static void gatheredgecoalescevertices(Edge e, Array<Vertex>& va);
static Face coalesce(Edge e);
// multi-use
static void outputedge(Edge e);
static void outputface(Face f);
static void outputmesh();
// unordered

int main(int argc, char const** argv)
{
	int froma3d;		// dummy
	Options opts;
	
	OPTSFC(opts,froma3d,	": build mesh from a3d input");
	opts.c("",":");
	OPTSDC(opts,renumber,	": renumber vertices and faces");
	OPTSDC(opts,outmesh,	": output mesh now");
	OPTSDC(opts,record,	": record mesh changes on stdout");
	OPTSFC(opts,nooutput,	": do not print mesh at program end");
	OPTSDC(opts,setb3d,	": set a3d output to binary");
	OPTSDC(opts,toa3d,	": output a3d version of mesh");
	OPTSDC(opts,tob3d,	": output binary a3d version of mesh");
	OPTSDC(opts,endobject,	": output EndObject marker");
	OPTSFC(opts,sloannormal,": use Sloan normals instead of limit subdiv");
	opts.c("",":");
	OPTSDC(opts,angle,	"deg : tag sharp edges");
	OPTSDC(opts,cosangle,	"fcos : tag sharp edges, acos(fcos)");
	OPTSDC(opts,solidangle,	"sterad. : tag cusp vertices (def. 0)");
	opts.c("",":");
	OPTSDC(opts,delaunay,	": retriangulate based on circumradii");
	OPTSDC(opts,diagonal,	": retriangulate based on diagonal lengths");
	OPTSDC(opts,segment,	": segment mesh and output a3d objects");
	OPTSDC(opts,selsmooth,	": selectively smooth, output a3d object");
	OPTSDC(opts,showsharpe,	": output sharp edges");
	opts.c("",": (next two are obsolete, use Subdivfit)");
	OPTSDC(opts,trisubdiv,	": 1 iter of triangular subdivision");
	OPTSDC(opts,silsubdiv,	": 1 iter of silhouette subdivision");
	OPTSDC(opts,expsubdiv,	"n : n iter of experimental subdivision");
	opts.c("",":");
	OPTSDC(opts,mark,	": mark tagged elements on output");
	OPTSDC(opts,removeinfo,	": remove tagged information at {v,f,e}");
	OPTSDC(opts,genus,	": print mesh genus");
	OPTSDC(opts,stat,	": print statistics about mesh");
	OPTSDC(opts,fixvertices,": disconnect non-nice vertices");
	OPTSDC(opts,fixfaces,	": remove contradictory faces");
	OPTSDC(opts,nice,	": assert mesh is nice");
	OPTSDC(opts,flip,	": flip orientation of faces");
	OPTSDC(opts,smootha3d,	": == -angle 180 -selsmooth");
	OPTSDC(opts,fillholes,	"nedges : fill holes with <= nedges");
	OPTSDC(opts,triangulate,": triangulate large faces (centersplit)");
	OPTSDC(opts,rmcomp,	"nfaces : remove components with <=nfaces");
	OPTSDC(opts,coalesce,	"fcrit : coalesce planar faces into polygons");
	opts.c("",":");
	OPTSPC(opts,nfaces,	"n :  stop when mesh has <=n faces");
	OPTSPC(opts,maxcrit,	"f :  stop when edge crit >=f");
	OPTSPC(opts,lengthc,	":  use edge length as reduction criterion");
	OPTSPC(opts,inscribedc,	":  use face inscribed radius as criterion");
	OPTSDC(opts,reduce,	": reduce mesh (must follow these args)");
	opts.c("",":");
	OPTSDC(opts,randpts,	"n : print n random points on mesh");
	OPTSDC(opts,vertexpts,	": print mesh vertices as points");
	OPTSDC(opts,bndpts,	"n : print n points on each boundary edge");
	OPTSDC(opts,addmesh,	": output a3d endfile + mesh now");
	opts.c("",":");
	
	TIMER(Filtermesh);
	if (argc>1 && !strcmp(argv[1],"-froma3d")) {
		argc--; argv++;
		dofroma3d();
	} else if (argc==1 || strcmp(argv[1],"-h")) {
		const char* s="-";
		if (argc>1 && argv[1][0]!='-') {
			s=argv[1];
			argc--; argv++;
		}
		RFile is(s);
		TIMER(_readmesh);
		for (;is().peek()=='#';) {
			char str[1000];
			is().get(str,sizeof(str));
			is().ignore();
			if (!str[1]) continue;
			SHOWFF(hform("%s\n",&str[2]));
		}
		mesh.read(is());
	}
	opts.allmustparse();
	if (!opts.parse(argc,argv) || !opts.noargs(argc,argv)) {
		opts.problem(argc,argv);
		mesh.clear();
		return 1;
	}
	ETIMER(Filtermesh);
	CleanUp();
	if (!nooutput) mesh.write(cout);
	mesh.clear();		// clear Pools
	return 0;
}

//*** froma3d

static void dofroma3d()
{
	TIMER(_froma3d);
	STAT(Sppdist2);
	RSA3dStream ia3d(cin);
	HashPoint hp;
	int ci=-1;
	Array<Vertex> gva;
	Array<Vertex> va;
	A3dElem el;
	for (;;) {
		ia3d.read(el);
		if (el.type()==A3dElem::TEndFile) break;
		if (el.type()==A3dElem::TComment) continue;
		if (assertw1(el.type()==A3dElem::TPolygon)) continue;
		va.init(el.num());
		for (int i=0;i<el.num();i++) {
			int k=hp.enter(el[i].p);
			if (k>ci) {
				ci=k;
				va[i]=mesh.createVertex();
				mesh.setPoint(va[i],el[i].p);
				gva.need(k+1);
				gva[k]=va[i];
			} else {
				va[i]=gva[k];
				Sppdist2+=dist2(el[i].p,mesh.point(va[i]));
			}
		}
		assertw1(mesh.createFace(va,el.num()));
	}
}

static void recordsharpe()
{
	assertx(cosangle!=UNDEF);
	assertw(cosangle>=-1 && cosangle<=1);
	STAT(Ssharp); STAT(Ssmooth); STAT(Sang);
	ForMeshEdge(mesh,e) {
		if (mesh.isBoundary(e)) continue;
		double angcos=EdgeDihedralAngleCos(mesh,e);
		double ang=todeg(myacos(angcos));
		Sang+=ang;
		if (angcos>cosangle) {
			Ssmooth+=ang;
			mesh.modFlag(e,GMesh::SHARPE,0);
		} else {
			Ssharp+=ang;
			mesh.modFlag(e,GMesh::SHARPE,1);
		}
	} EndFor;
}

static void recordcuspv()
{
	ForMeshVertex(mesh,v) {
		if (mesh.isBoundary(v)) continue;
		int iscusp=VertexSolidAngle(mesh,v)<solidangle;
		mesh.modFlag(v,GMesh::CUSPV,iscusp);
	} EndFor;
}

static int sharp(Edge e)
{
	return mesh.isBoundary(e) || mesh.flag(e,GMesh::SHARPE);
}

static int doangle(int argc, char const** argv)
{
	assertx(argc>1);
	double angle=atof(argv[1]);
	assertx(angle>=0 && angle<=180);
	cosangle=cos(torad(angle)),recordsharpe();
	return 2;
}

static int docosangle(int argc, char const** argv)
{
	assertx(argc>1);
	cosangle=atof(argv[1]);
	recordsharpe();
	return 2;
}

static int dosolidangle(int argc, char const** argv)
{
	assertx(argc>1);
	solidangle=atof(argv[1]);
	assertx(solidangle>=0 && solidangle<=4*PI);
	recordcuspv();
	return 2;
}

static void dosetb3d()
{
	setenv("A3DBINARY","",1);
}


//*** toa3d, tob3d

static void dotoa3d()
{
	TIMER(_toa3d);
	nooutput=1;
	outputmesh();
}

static void dotob3d()
{
	TIMER(_tob3d);
	setenv("A3DBINARY","",1);
	nooutput=1;
	outputmesh();
}

static void doendobject()
{
	oa3d.writeEndObject();
}

//*** other

static void dorenumber()
{
	TIMER(_renumber);
	mesh.renumber();
}

static void dooutmesh()
{
	TIMER(_outmesh);
	mesh.write(cout);
}

static void doaddmesh()
{
	TIMER(_addmesh);
	A3dElem el(A3dElem::TEndFile);
	oa3d.write(el);
	mesh.write(cout);
}

static void dorecord()
{
	nooutput=1;
	mesh.recordChanges(&cout);
}

static void domark()
{
	ForMeshVertex(mesh,v) {
		mesh.setInfo(v,mesh.flag(v,GMesh::CUSPV)?
			     new MeshInfoString("cusp"):0);
	} EndFor;
	ForMeshEdge(mesh,e) {
		mesh.setInfo(e,mesh.flag(e,GMesh::SHARPE)?
			     new MeshInfoString("sharp"):0);
	} EndFor;
}

//*** delaunay

static void dodelaunay()
{
	TIMER(_delaunay);
	assertx(cosangle!=UNDEF);
	int ns=RetriangulateAll(mesh,cosangle,CircumRadiusSwapCrit,0,0);
	SHOWDF("Swapped %d edges\n",ns);
}

//*** diagonal

static void dodiagonal()
{
	TIMER(_diagonal);
	assertx(cosangle!=UNDEF);
	int ns=RetriangulateAll(mesh,cosangle,DiagonalDistanceSwapCrit,0,0);
	SHOWDF("Swapped %d edges\n",ns);
}

//*** segment

static void dosegment()
{
	TIMER(_segment);
	nooutput=1;
	Map<Face,int> mfseg;	// Face -> segment#
	Stack<Face> stackrepf;
	gathersegments(mfseg,stackrepf);
	edgeanglestats(mfseg);
	SHOWDF("Detected %d segments\n",stackrepf.height());
	while (!stackrepf.empty()) {
		printsegment(stackrepf.pop(),mfseg);
		oa3d.writeEndObject();
	}
	assertx(mfseg.empty());
}

static void gathersegments(Map<Face,int>& mfseg, Stack<Face>& stackrepf)
{
	HPqueue<Face> pq;
	NEST {
		STAT(Sseg);
		int segnum=0;
		ForMeshFace(mesh,f) {
			if (mfseg.contains(f)) continue;
			segnum++;
			int nf; Point pc;
			gatherfollowseg(f,mfseg,segnum,nf,pc);
			Sseg+=nf;
			double pri=fabs(pc[0]*37+pc[1]*17+pc[2]);
			pq.enter(f,pri);
		} EndFor;
	}
	for (;!pq.empty();) {
		Face f=pq.removemin();
		stackrepf.push(f);
	}
}

static void gatherfollowseg(Face f, Map<Face,int>& mfseg,
			    int segnum, int& nf, Point& pcentroid)
{
	int nfaces=0;
	Homogeneous h;
	Queue<Face> queue;
	mfseg.enter(f,segnum);
	Polygon poly;
	for (;;) {
		nfaces++;
		mesh.polygon(f,poly);
		h+=centroid(poly)*poly.getarea();
		ForFaceEdge(mesh,f,e) {
			if (sharp(e)) continue;
			Face f2=mesh.oppFace(f,e);
			if (!mfseg.contains(f2)) {
				mfseg.enter(f2,segnum);
				queue.enqueue(f2);
			}
		} EndFor;
		if (queue.empty()) break;
		f=queue.dequeue();
	}
	nf=nfaces;
	pcentroid=toPoint(normalize(h));
}

static void edgeanglestats(const Map<Face,int>& mfseg)
{
	STAT(Sbndang); STAT(Sintang);
	ForMeshEdge(mesh,e) {
		if (mesh.isBoundary(e)) continue;
		double angcos=EdgeDihedralAngleCos(mesh,e);
		double ang=todeg(myacos(angcos));
		Face f1=mesh.face1(e), f2=mesh.face2(e);
		if (mfseg.get(f1)==mfseg.get(f2))
			Sintang+=ang;
		else
			Sbndang+=ang;
	} EndFor;
}

static void printsegment(Face f, Map<Face,int>& mfseg)
{
	Queue<Face> queue;
	assertx(mfseg.remove(f));
	for (;;) {
		outputface(f);
		ForFaceEdge(mesh,f,e) {
			if (sharp(e)) { outputedge(e); continue; }
			Face f2=mesh.oppFace(f,e);
			if (mfseg.remove(f2)) queue.enqueue(f2);
		} EndFor;
		if (queue.empty()) break;
		f=queue.dequeue();
	}
}

//*** selsmooth

static void doselsmooth()
{
	TIMER(_selsmooth);
	nooutput=1;
	selsmooth();
}

static void selsmooth()
{
	Map<Vertex,Vnors*> mvnors;
	NEST {
		STAT(Svncomp);
		ForMeshVertex(mesh,v) {
			Vnors* vnors=new Vnors;
			mvnors.enter(v,vnors);
			Svncomp+=ComputeVnors(mesh,v,*vnors,!sloannormal);
		} EndFor;
	}
	Set<Face> fvisited;
	Queue<Face> queue;
	ForMeshFace(mesh,ff) {
		Face f=ff;
		if (!fvisited.add(f)) continue;
		oa3d.writeComment(" beginning of new segment");
		for (;;) {
			static A3dElem el;
			el.init(A3dElem::TPolygon);
			ForFaceVertex(mesh,f,v) {
				const Vector& nor=mvnors.get(v)->getnor(f);
				el+=A3dVertex(mesh.point(v),nor,SURFCOL);
			} EndFor;
			oa3d.write(el);
			ForFaceEdge(mesh,f,e) {
				if (sharp(e)) continue;
				Face f2=mesh.oppFace(f,e);
				if (fvisited.add(f2)) queue.enqueue(f2);
			} EndFor;
			if (queue.empty()) break;
			f=queue.dequeue();
		}
	} EndFor;
	oa3d.flush();
	ForMapValue(mvnors,Vertex,Vnors*,vnors) { delete vnors; } EndFor;
}

static void doshowsharpe()
{
	nooutput=1;
	int ne=0;
	ForMeshEdge(mesh,e) {
		if (sharp(e)) outputedge(e),ne++;
	} EndFor;
	SHOWDF("Output %d sharp edges\n",ne);
}

//*** trisubdiv

static void dotrisubdiv()
{
	TIMER(_trisubdiv);
	meshsubdiv();
}

// Create a new mesh in situ with the old vertices
// (each old vertex will temporarily have 2 rings of faces!)
static void meshsubdiv()
{
        int i;
	Stack<Vertex> stackv; Stack<Face> stackf; Stack<Edge> stacke;
	ForMeshVertex(mesh,v) { stackv.push(v); } EndFor;
	ForMeshFace(mesh,f) { stackf.push(f); } EndFor;
	ForMeshEdge(mesh,e) { stacke.push(e); } EndFor;
	Map<Edge,Vertex> menewv;
	// Create new vertices and compute their positions
	ForStack(stacke,Edge,e) {
		Vertex v=mesh.createVertex();
		menewv.enter(e,v);
		Homogeneous h=(3*mesh.point(mesh.vertex1(e))+
			       3*mesh.point(mesh.vertex2(e)));
		if (!sharp(e))
			h+=(mesh.point(mesh.sideVertex1(e))+
			    mesh.point(mesh.sideVertex2(e)));
		mesh.setPoint(v,toPoint(h.normalize()));
	} EndFor;
	// Update old vertex positions
	Map<Vertex,Point*> mapvp;
	ForStack(stackv,Vertex,v) {
		int ne=0, nesharp=0;
		Homogeneous h, hsharp;
		ForVertexEdge(mesh,v,e) {
			const Point& p=mesh.point(mesh.oppVertex(v,e));
			ne++,h+=p;
			if (sharp(e)) nesharp++,hsharp+=p;
		} EndFor;
		if (nesharp>2) { // vertex position held constant
			continue;
		} else if (nesharp==2) { // cubic bspline boundary
			h=hsharp+6*mesh.point(v);
		} else {	// quartic bspline interior
			// (if nesharp==1, treat as interior vertex)
			int n=ne;
			double a=5./8.-square((3+2*cos(2*PI/n))/8);
			double centerw=n*(1-a)/a;
			h+=centerw*mesh.point(v);
		}
		mapvp.enter(v,new Point(toPoint(h.normalize())));
	} EndFor;
	ForMapKeyValue(mapvp,Vertex,v,Point*,p) {
		mesh.setPoint(v,*p);
		delete p;
	} EndFor;
	// Create new triangulation
	ForStack(stackf,Face,f) {
		Vertex va[3], vs[3];
		mesh.vertices(f,va);
//		for (int i=0;i<3;i++) {
		for (i=0;i<3;i++) {
			Edge e=mesh.edge(va[i],va[(i+1)%3]);
			vs[i]=menewv.get(e);
		}
		for (i=0;i<3;i++)
			assertx(mesh.createFace(va[i],vs[i],vs[(i+2)%3]));
		assertx(mesh.createFace(vs,3));
	} EndFor;
	// Update sharp edges
	ForStack(stacke,Edge,e) {
		if (!mesh.flag(e,GMesh::SHARPE)) continue;
		Vertex vnew=menewv.get(e);
		mesh.modFlag(mesh.edge(vnew,mesh.vertex1(e)),GMesh::SHARPE,1);
		mesh.modFlag(mesh.edge(vnew,mesh.vertex2(e)),GMesh::SHARPE,1);
	} EndFor;
	// Remove old triangulation
	ForStack(stackf,Face,f) { mesh.destroyFace(f); } EndFor;
}

//*** silsubdiv

struct vvnewv {
	vvnewv(Vertex pv1, Vertex pv2, Vertex pvnew, int psharp)
		: vnew(pvnew), sharp(psharp)
	{ if (pv1<pv2) v1=pv1,v2=pv2; else v1=pv2,v2=pv1; }
	Vertex v1, v2, vnew;
	int sharp;
};

static Univ hashvvnewv(const vvnewv* s)
{
	return Conv<int>::e(mesh.vertexid(s->v1)+mesh.vertexid(s->v2)*761);
}

static int cmpvvnewv(const vvnewv* s1, const vvnewv* s2)
{
	return s1->v1!=s2->v1 || s1->v2!=s2->v2;
}

static void dosilsubdiv()
{
        int i;
	TIMER(_silsubdiv);
	Stack<Face> stackf;
	ForMeshFace(mesh,f) { stackf.push(f); } EndFor;
	// Determine which edges will be subdivided
	Set<Edge> subde;	// edges to subdivide
	ForMeshEdge(mesh,e) {
		if (sharp(e)) subde.enter(e);
	} EndFor;
	Queue<Face> queuef;
	ForStack(stackf,Face,f) { queuef.enqueue(f); } EndFor;
	while (!queuef.empty()) {
		Face f=queuef.dequeue();
		int nnew=0;
		ForFaceEdge(mesh,f,e) {
			if (subde.contains(e)) nnew++;
		} EndFor;
		if (nnew!=2) continue; // ok, no propagating changes
		ForFaceEdge(mesh,f,e) {
			if (!subde.add(e)) continue;
			Face f2=mesh.oppFace(f,e);
			if (f2) queuef.enqueue(f2);
		} EndFor;
	}
	// Introduce new vertices at midpoints
	HashStruct<vvnewv> mvvnewv(hashvvnewv,cmpvvnewv);
	ForSet(subde,Edge,e) {
		Vertex v=mesh.createVertex();
		mvvnewv.enter(new vvnewv(mesh.vertex1(e),mesh.vertex2(e),v,
					 mesh.flag(e,GMesh::SHARPE)));
		mesh.setPoint(v,interp(mesh.point(mesh.vertex1(e)),
				       mesh.point(mesh.vertex2(e))));
	} EndFor;
	// Subdivide faces, destroys validity of flags(e)
	ForStack(stackf,Face,f) {
		Vertex va[3], vs[3]; int nnew=0, i0=-1;
		mesh.vertices(f,va);
//		for (int i=0;i<3;i++) {
		for (i=0;i<3;i++) {
			Edge e=mesh.edge(va[i],va[(i+1)%3]);
			vvnewv si(mesh.vertex1(e),mesh.vertex2(e),0,0);
			vvnewv* s=mvvnewv.retrieve(&si);
			vs[i]=s?s->vnew:0;
			if (vs[i]) nnew++,i0=i;
		}
		if (!nnew) continue;
		assertx(nnew==1 || nnew==3);
		mesh.destroyFace(f);
		if (nnew==1) {
			int i1=(i0+1)%3, i2=(i0+2)%3;
			assertx(mesh.createFace(va[i0],vs[i0],va[i2]));
			assertx(mesh.createFace(vs[i0],va[i1],va[i2]));
			continue;
		}
		for (i=0;i<3;i++)
			assertx(mesh.createFace(va[i],vs[i],vs[(i+2)%3]));
		assertx(mesh.createFace(vs,3));
	} EndFor;
	ForHashStruct(mvvnewv,vvnewv,s) {
		if (s->sharp) {
			mesh.modFlag(mesh.edge(s->vnew,s->v1),GMesh::SHARPE,1);
			mesh.modFlag(mesh.edge(s->vnew,s->v2),GMesh::SHARPE,1);
		}
		delete s;
	} EndFor;
	mvvnewv.clear();	// avoid warning message
	// Update vertex positions
	ForMeshVertex(mesh,v) {
		int nesharp=0;
		Homogeneous h;
		ForVertexEdge(mesh,v,e) {
			if (sharp(e))
				nesharp++,h+=mesh.point(mesh.oppVertex(v,e));
		} EndFor;
		if (nesharp!=2) continue;
		// important note: here order does not matter because newly
		// inserted midpoints remain midpoints.
		// Homogeneous() cast for DECCXX
		mesh.setPoint(v,toPoint(Homogeneous(.5*mesh.point(v))+
					Homogeneous(.25*h)));
	} EndFor;
}

//*** expsubdiv

struct faceinfo : MeshInfo {
	faceinfo(const Point& pp) : p(pp) { }
	Point p;
	Point newp;
};

static void expsubdiv()
{
	Stack<Face> stackf; ForMeshFace(mesh,f) { stackf.push(f); } EndFor;
	Stack<Edge> stacke; ForMeshEdge(mesh,e) { stacke.push(e); } EndFor;
	Map<Edge,Vertex> menewv;
	// Create new vertices
	ForStack(stacke,Edge,e) {
		menewv.enter(e,mesh.createVertex());
	} EndFor;
	int try1=GetenvValue("TRY1");
	// Create new triangulation
	ForStack(stackf,Face,f) {
		const Point& p=reinterpret_cast<faceinfo*>(mesh.info(f))->p;
		Vertex va[3], vs[3]; mesh.vertices(f,va);
		Face fs[3];
		ForIndex(i,3) {
			Edge e=mesh.edge(va[i],va[(i+1)%3]);
			vs[i]=menewv.get(e);
			fs[i]=mesh.oppFace(f,e); // may be 0
		} EndFor;
		ForIndex(i,3) {
			Face f=assertv(mesh.createFace(va[i],
						       vs[i],vs[(i+2)%3]));
			if (try1 && fs[0] && fs[1] && fs[2]) {
				Homogeneous h(p*.75);
				h+=reinterpret_cast<faceinfo*>(mesh.info(fs[(i+1)%3]))->p*-.25;
				h+=reinterpret_cast<faceinfo*>(mesh.info(fs[(i+0)%3]))->p*.25;
				h+=reinterpret_cast<faceinfo*>(mesh.info(fs[(i+2)%3]))->p*.25;
				mesh.setInfo(f,new faceinfo(toPoint(h)));
			} else {
				mesh.setInfo(f,new faceinfo(p));
			}
		} EndFor;
		Face f=assertv(mesh.createFace(vs,3));
		mesh.setInfo(f,new faceinfo(p));
	} EndFor;
	// Update sharp edges
	ForStack(stacke,Edge,e) {
		if (!mesh.flag(e,GMesh::SHARPE)) continue;
		Vertex vnew=menewv.get(e);
		mesh.modFlag(mesh.edge(vnew,mesh.vertex1(e)),GMesh::SHARPE,1);
		mesh.modFlag(mesh.edge(vnew,mesh.vertex2(e)),GMesh::SHARPE,1);
	} EndFor;
	// Remove old triangulation
	ForStack(stackf,Face,f) { mesh.destroyFace(f); } EndFor;
	// Update positions
	double wcenter=1;
	if (getenv("WCENTER")) wcenter=atof(getenv("WCENTER"));
	double woutside=1.;
	int nconv=2;
	if (getenv("NCONV")) nconv=atoi(getenv("NCONV"));
	if (try1) nconv=0;
	ForIndex(ni,nconv) {
		ForMeshFace(mesh,f) {
			faceinfo* fi=reinterpret_cast<faceinfo*>(mesh.info(f));
			Homogeneous h(fi->p*wcenter);
			ForFaceFace(mesh,f,ff) {
				h+=reinterpret_cast<faceinfo*>(mesh.info(ff))->p*woutside;
			} EndFor;
			fi->newp=toPoint(normalize(h));
		} EndFor;
		ForMeshFace(mesh,f) {
			faceinfo* fi=reinterpret_cast<faceinfo*>(mesh.info(f));
			fi->p=fi->newp;
		} EndFor;
	} EndFor;
}

static int doexpsubdiv(int argc, char const** argv)
{
	assertx(argc>1); int niter=atoi(argv[1]);
	ForMeshFace(mesh,f) {
		Homogeneous h;
		ForFaceVertex(mesh,f,v) { h+=mesh.point(v); } EndFor;
		mesh.setInfo(f,new faceinfo(toPoint(normalize(h))));
	} EndFor;
	ForIndex(ni,niter) { expsubdiv(); } EndFor;
	ForMeshVertex(mesh,v) {
		Homogeneous h;
		ForVertexFace(mesh,v,f) {
			h+=reinterpret_cast<faceinfo*>(mesh.info(f))->p;
		} EndFor;
		mesh.setPoint(v,toPoint(normalize(h)));
	} EndFor;
	if (0) {
		// this will provide an approximation to the refined mesh,
		// with vertices on a subdivision of the original mesh.
	} else if (0) {
		Stack<Face> stackf;
		ForMeshFace(mesh,f) { stackf.push(f); } EndFor;
		ForStack(stackf,Face,f) {
			const Point p=reinterpret_cast<faceinfo*>(mesh.info(f))->p;
			Vertex va[3]; mesh.vertices(f,va);
			mesh.destroyFace(f);
			Vertex v=mesh.createVertex();
			mesh.setPoint(v,p);
			ForIndex(i,3) {
				assertx(mesh.createFace(v,va[i],va[(i+1)%3]));
			} EndFor;
		} EndFor;
	} else {
		ForMeshVertex(mesh,v) {
			static A3dElem el;
			el.init(A3dElem::TPolygon,0);
			ForVertexCcwVertex(mesh,v,vv) {
				Vertex vvv=mesh.ccwVertex(v,vv);
				if (!vvv) continue;
				Face f=mesh.face(v,vv,vvv);
				Point p=reinterpret_cast<faceinfo*>(mesh.info(f))->p;
				el+=A3dVertex(p,Vector(0,0,0),SURFCOL);
			} EndFor;
			oa3d.write(el);
		} EndFor;
		nooutput=1;
	}
	return 2;
}

//*** fillholes

static int dofillholes(int argc, char const** argv)
{
	TIMER(_fillholes);
	assertx(argc>1);
	int maxnume=atoi(argv[1]); // ==maxnumv
	Set<Edge> setbe;	   // boundary edges
	ForMeshEdge(mesh,e) {
		if (mesh.isBoundary(e)) setbe.enter(e);
	} EndFor;
	STAT(Sbndlen); STAT(Sbndsub);
	for (;!setbe.empty();) {
		Edge e=setbe.getone();
		Queue<Edge> queuee;
		GatherBoundary(mesh,e,queuee);
		ForQueue(queuee,Edge,e) {
			assertx(setbe.remove(e));
		} EndFor;
		int ne=queuee.length();
		if (ne>maxnume) continue;
		Sbndlen+=ne;
		Sbndsub+=removeboundary(queuee);
	}
	SHOWDF("Filled in %d holes\n",Sbndlen.num());
	return 2;
}

static int removeboundary(Queue<Edge>& queuee)
{
	int nsub=0;
	for (;!queuee.empty();) {
		nsub++;
		Set<Vertex> setv;
		Vertex v=0;
		ForQueue(queuee,Edge,e) {
			Vertex vv=mesh.vertex2(e);
			if (!setv.add(vv)) { v=vv; break; }
		} EndFor;
		Queue<Edge> qc;
		if (!v) {
			qc.addtoend(queuee);
		} else {
			// Rotate queuee to put v at front
			while (mesh.vertex2(queuee.front())!=v)
				queuee.enqueue(queuee.dequeue());
			// Extract loop from queuee into qc
			for (;;) {
				Edge e=queuee.dequeue();
				qc.enqueue(e);
				if (mesh.vertex1(e)==v) break;
			}
		}
		removesimpleboundary(qc);
	}
	return nsub;
}

static void removesimpleboundary(const Queue<Edge>& queuee)
{
	int nv=queuee.length();
	Array<Vertex> va;
	ForQueue(queuee,Edge,e) {
		va+=mesh.vertex2(e); // vertex1(e) ok but slower
	} EndFor;
	SHOWDF(" filling in hole with %d sides\n",nv);
	assertx(mesh.createFace(va,nv));
}

//*** triangulate

static void dotriangulate()
{
	STAT(Sfverts);
	Stack<Face> stackf;
	ForMeshFace(mesh,f) {
		int nv=mesh.numVertices(f);
		if (nv==3) continue;
		Sfverts+=nv;
		stackf.push(f);
	} EndFor;
	SHOWDF("Found %d faces to triangulate\n",
	       stackf.height());
	while (!stackf.empty()) CenterSplitFace(mesh,stackf.pop());
	if (!Sfverts.num()) Sfverts.setprint(0);
}

//*** rmcomponents

static int dormcomp(int argc, char const** argv)
{
	// component is assumed face-face connected (do not jump bowtie)
	TIMER(_rmcomponents);
	assertx(argc>1);
	int maxnumf=atoi(argv[1]);
	Set<Face> setfvis;	// faces already considered
	STAT(Svertsrem); STAT(Sfacesrem);
	for (;;) {
		int found=0;
		ForMeshFace(mesh,f) {
			if (setfvis.contains(f)) continue;
			Set<Face> setf;
			GatherComponent(mesh,f,setf);
			ForSet(setf,Face,ff) {
				setfvis.enter(ff);
			} EndFor;
			int nf=setf.num();
			if (nf>maxnumf) continue;
			Sfacesrem+=nf;
			Svertsrem+=removecomponent(setf);
			found=1;
			break;
		} EndFor;
		if (!found) break;
	}
	SHOWDF("Removed %d mesh components\n",Sfacesrem.num());
	return 2;
}

static int removecomponent(const Set<Face>& setf)
{
	Set<Vertex> setv;
	ForSet(setf,Face,f) {
		ForFaceVertex(mesh,f,v) {
			setv.add(v); // may already be there
		} EndFor;
		mesh.destroyFace(f);
	} EndFor;
	int nvrem=0;
	ForSet(setv,Vertex,v) {
		if (mesh.degree(v)) continue;
		nvrem++;
		mesh.destroyVertex(v);
	} EndFor;
	return nvrem;
}

//*** coalesce

// This version of coalesce has a weakness.  It only attempts to merge 2 faces
// by removing one edge common to them.  But, after several coalescences, 2
// faces may share more than one edge.
static int docoalesce(int argc, char const** argv)
{
	TIMER(_coalesce);
	assertx(argc>1);
	double fcrit=atof(argv[1]);
	int nerem=0;
	Set<Edge> sete;
	ForMeshEdge(mesh,e) {
		sete.enter(e);
	} EndFor;
	for (;!sete.empty();) {
		Edge e=sete.removeone();
		double f=trycoalesce(e);
		if (f>fcrit) continue;
		// all systems go
		// neighboring edges may change, so remove them from sete
		ForEdgeFace(mesh,e,f) {
			ForFaceEdge(mesh,f,ee) {
				(void)sete.remove(ee);
			} EndFor;
		} EndFor;
		// coalesce faces into new face fnew
		Face fnew=coalesce(e);
		nerem++;
		// reenter affected edges
		ForFaceEdge(mesh,fnew,ee) {
			sete.enter(ee);	// cannot still be there
		} EndFor;
	}
	SHOWDF("Removed %d edges\n",nerem);
	return 2;
}

static void gatheredgecoalescevertices(Edge e, Array<Vertex>& va)
{
	static Array<Vertex> va1, va2;
	mesh.vertices(mesh.face1(e),va1);
	mesh.vertices(mesh.face2(e),va2);
	int nv1=va1.num(), nv2=va2.num();
	int i1, i2, ic, i;
	// Find one vertex common to both faces (mesh.vertex1(e))
	for (i1=0;i1<nv1;i1++) if (va1[i1]==mesh.vertex1(e)) break;
	for (i2=0;i2<nv2;i2++) if (va2[i2]==mesh.vertex1(e)) break;
	assertx(i1<nv1 && i2<nv2);
	// Find most clw vertex on face1 common to both
	for (;va1[(i1-1+nv1)%nv1]==va2[(i2+1)%nv2];)
		i1=(i1-1+nv1)%nv1,i2=(i2+1)%nv2;
	// Let ic be the number of vertices common to both faces
	for (ic=1;;ic++) if (va1[(i1+ic)%nv1]!=va2[(i2-ic+nv2)%nv2]) break;
	va.init(0);
	for (i=ic-1;i<nv1;i++) va+=va1[(i1+i)%nv1];
	for (i=0;i<nv2-ic+1;i++) va+=va2[(i2+i)%nv2];
}

// return 1e30 if not legal
static double trycoalesce(Edge e)
{
	if (mesh.isBoundary(e)) return 1e30;
	static Array<Vertex> va;
	gatheredgecoalescevertices(e,va);
	NEST {			// check for duplicate vertices
		Set<Vertex> vset;
		for (int i=0;i<va.num();i++)
			if (!vset.add(va[i])) return 1e30;
	}
	static Array<Point> pa; pa.init(va.num());
	for (int i=0;i<va.num();i++) pa[i]=mesh.point(va[i]);
	Frame frame; double eimag[3];
	PrincipalComponents(pa,va.num(),frame,eimag);
	// ratio of thickness to width
	// width judged better than length (two adjacent skinny triangles)
	return eimag[2]/eimag[1];
}

static Face coalesce(Edge e)
{
	static Array<Vertex> va;
	gatheredgecoalescevertices(e,va);
	Face f1=mesh.face1(e), f2=mesh.face2(e);
	Set<Vertex> vbefore;
	ForFaceVertex(mesh,f1,v) { vbefore.enter(v); } EndFor;
	ForFaceVertex(mesh,f2,v) { vbefore.add(v); } EndFor;
	mesh.destroyFace(f1);
	mesh.destroyFace(f2);
	for (int i=0;i<va.num();i++) (void)vbefore.remove(va[i]);
	ForSet(vbefore,Vertex,v) { mesh.destroyVertex(v); } EndFor;
	return assertv(mesh.createFace(va,va.num()));
}

//*** other

static void donice()
{
	assertx(mesh.isNice());
}

static void doflip()
{
	Queue<Array<Vertex>*> qva;
	ForMeshOrderedFace(mesh,f) {
		Array<Vertex>& va=*new Array<Vertex>;
		mesh.vertices(f,va);
		qva.enqueue(&va);
	} EndFor;
	Stack<Face> sf; ForMeshFace(mesh,f) { sf.push(f); } EndFor;
	ForStack(sf,Face,f) { mesh.destroyFace(f); } EndFor;
	Array<Vertex> vnew;
	ForQueue(qva,Array<Vertex>*,pva) {
		Array<Vertex>& va=*pva;
		int nv=va.num();
		vnew.init(nv);
		ForIndex(i,nv) { vnew[i]=va[nv-1-i]; } EndFor;
		assertx(mesh.createFace(vnew,nv));
		delete pva;
	} EndFor;
	mesh.renumber();
}

static void dofixvertices()
{
	STAT(Svnrings);
	Queue<Vertex> queuev;
	ForMeshVertex(mesh,v) {
		if (!mesh.isNice(v)) queuev.enqueue(v);
	} EndFor;
	while (!queuev.empty())
		Svnrings+=fixvertex(queuev.dequeue());
	SHOWDF("Fixed %d vertices\n",Svnrings.num());
}

// Return the number of rings the old vertex had.
static int fixvertex(Vertex v)
{
	int nrings=0;
	Set<Face> setallf;
	ForVertexFace(mesh,v,f) {
		setallf.enter(f);
	} EndFor;
	for (;;) {
		nrings++;
		Face frep=setallf.getone();
		Set<Face> setf;
		setf.enter(frep);
		For (Face f=frep;;) {
			f=mesh.clwFace(v,f);
			if (!f || !setf.add(f)) break;
		} EndFor;
		For (Face f=frep;;) {
			f=mesh.ccwFace(v,f);
			if (!f || !setf.add(f)) break;
		} EndFor;
		ForSet(setf,Face,f) {
			assertx(setallf.remove(f));
		} EndFor;
		if (setallf.empty()) break; // is now a nice vertex
		Vertex vnew=mesh.createVertex();
		mesh.setPoint(vnew,mesh.point(v)); // same geometry
		ForSet(setf,Face,f) {
			assertx(mesh.substituteFaceVertex(f,v,vnew));
		} EndFor;
	}
	return nrings;
}

static void dofixfaces()
{
	int nf=0;
	Queue<Face> queuef;
	ForMeshFace(mesh,f) {
		if (!mesh.isNice(f)) queuef.enqueue(f);
	} EndFor;
	ForQueue(queuef,Face,f) {
		if (mesh.isNice(f)) continue; // other face may be gone now
		nf++;
		mesh.destroyFace(f);
	} EndFor;
	SHOWDF("Fixed %d faces\n",nf);
}

static void dosmootha3d()
{
	TIMER(_smootha3d);
	nooutput=1;
	// Not too efficient because normals are completely shared at vertices.
	selsmooth();
}

static void dogenus()
{
	SHOWDF(MeshGenusString(mesh));
}

static void doremoveinfo()
{
	ForMeshVertex(mesh,v) {
		mesh.setInfo(v,0); mesh.modFlag(v,GMesh::ALL,0);
	} EndFor;
	ForMeshFace(mesh,f) {
		mesh.setInfo(f,0); mesh.modFlag(f,GMesh::ALL,0);
	} EndFor;
	ForMeshEdge(mesh,e) {
		mesh.setInfo(e,0); mesh.modFlag(e,GMesh::ALL,0);
	} EndFor;
}

static void dostat()
{
	SHOWDF(MeshGenusString(mesh));
	NEST {
		STAT(Sbound);
		MeshStatBoundaries(mesh,Sbound);
	}
	NEST {
		STAT(Scompf);
		MeshStatComponents(mesh,Scompf);
	}
	NEST {
		STAT(Snormsolida); STAT(Ssolidang);
		ForMeshVertex(mesh,v) {
			if (mesh.numBoundaries(v)) continue;
			double solidang=VertexSolidAngle(mesh,v);
			Ssolidang+=solidang;
			Snormsolida+=1-solidang/(2*PI);
		} EndFor;
	}
	NEST {
		STAT(Sfacearea);
		Polygon poly;
		ForMeshFace(mesh,f) {
			mesh.polygon(f,poly);
			Sfacearea+=poly.getarea();
		} EndFor;
	}
	NEST {
		STAT(Sfvertices);
		ForMeshFace(mesh,f) {
			Sfvertices+=mesh.numVertices(f);
		} EndFor;
	}
	NEST {
		STAT(Sbvalence); STAT(Sivalence); STAT(Svalence);
		ForMeshVertex(mesh,v) {
			Svalence+=mesh.degree(v);
			if (!mesh.isNice(v)) continue;
			if (mesh.isBoundary(v)) Sbvalence+=mesh.degree(v);
			else Sivalence+=mesh.degree(v);
		} EndFor;
	}
}

//*** reduce

static void doreduce()
{
	TIMER(_reduce);
	assertx(reducecrit);
	// also use: nfaces, maxcrit
	HPqueue<Edge> pqe;
	NEST {
		TIMER(__initpq);
		STAT(Sred);	// optional
		ForMeshEdge(mesh,e) {
			pqe.enter(e,Sred.enter(reducecriterion(e)));
		} EndFor;
	}
	int nf=mesh.numFaces(), orignf=nf;
	int ne=mesh.numEdges(), origne=ne;
	int ncol=0;
	for (;;) {
		if (nf<=nfaces) break;
		if (pqe.empty()) break;
		double crit=pqe.minpriority();
		Edge e=pqe.removemin();
		if (crit>maxcrit) break;
		if (!mesh.niceEdgeCollapse(e)) continue;
		// Do edge collapse
		ForEdgeVertex(mesh,e,v) {
			ForVertexEdge(mesh,v,e2) {
				(void)pqe.remove(e2);
			} EndFor;
		} EndFor;
		int nfcol=1+(mesh.face2(e)?1:0);
		nf-=nfcol;
		ne-=1+nfcol;
		ncol++;
		Vertex vkept=mesh.vertex1(e);
		mesh.collapseEdge(e);
		ForVertexEdge(mesh,vkept,e) {
			pqe.enterupdate(e,reducecriterion(e));
		} EndFor;
		ForVertexFace(mesh,vkept,f) {
			Edge e=mesh.oppEdge(vkept,f);
			pqe.enterupdate(e,reducecriterion(e));
		} EndFor;
	}
	SHOWDF("Reduced %d times, deleted %d edges, %d faces\n",
	       ncol,origne-ne,orignf-nf);
}

static double reducecriterion(Edge e)
{
	if (reducecrit==Length) return mesh.length(e);
	if (reducecrit==Inscribed) return CollapseEdgeCriterion(mesh,e);
	assertnever(""); return 0;
}

static void lengthc() { reducecrit=Length; }

static void inscribedc() { reducecrit=Inscribed; }

//*** point sampling

static int dorandpts(int argc, char const** argv)
{
	TIMER(_randpts);
	nooutput=1;
	assertx(argc>1);
	int npoints=atoi(argv[1]);
	int nf=mesh.numFaces();
	SArray<Face> fface(nf);	// Face of this index
	SArray<double> farea(nf); // area of face
	SArray<double> fcarea(nf+1); // cumulative area
	double sumarea=0;	    // for accuracy
	int i=0;
	Polygon poly;
	ForMeshFace(mesh,f) {
		mesh.polygon(f,poly);
		if (assertw(poly.num()==3)) continue;
		fface[i]=f;
		sumarea+=(farea[i]=poly.getarea());
		i++;
	} EndFor;
	int nfaces=i;
	SHOWDF("Total area %g over %d faces\n",sumarea,nfaces);
	double area=0;		// for accuracy
	for (i=0;i<nfaces;i++) {
		fcarea[i]=area;
		area+=farea[i]/sumarea;
	}
	assertx(fabs(area-1)<1e-5);
	fcarea[nfaces]=1.00001;
	for (i=0;i<npoints;i++) {
		int fi=mybsearch(fcarea,0,nfaces,Random::G.unif());
		assertx(fi>=0 && fi<nfaces); // optional
		Face f=fface[fi];
		mesh.polygon(f,poly);
		outputpoint(trianglerandompoint(poly));
	}
	SHOWDF("Printed %d random points\n",npoints);
	return 2;
}

// given ar[l]<=v<ar[h]
static int mybsearch(const double* ar, int l, int h, double v)
{
	assertx(l<h);		// optional 3
	assertx(ar[l]<=v);
	assertx(v<ar[h]);
	if (h-l==1) return l;
	int m=(l+h)/2;
	if (v>=ar[m]) return mybsearch(ar,m,h,v);
	else return mybsearch(ar,l,m,v);
}

static Point trianglerandompoint(const Polygon& poly)
{
	assertx(poly.num()==3);
	double a=Random::G.unif();
	double b=Random::G.unif();
	if (a+b>1) a=1-a,b=1-b;
	return interp(poly[0],poly[1],poly[2],a,b);
}

static void outputpoint(const Point& p)
{
	static A3dElem el;
	el.init(A3dElem::TPoint,0);
	el+=A3dVertex(p,Vector(0,0,0),A3dVertexColor::White);
	oa3d.write(el);
}

static void dovertexpts()
{
	TIMER(_vertexpts);
	nooutput=1;
	ForMeshVertex(mesh,v) {
		outputpoint(mesh.point(v));
	} EndFor;
	SHOWDF("Printed %d vertex points\n",mesh.numVertices());
}

static int dobndpts(int argc, char const** argv)
{
	TIMER(_bndpts);
	nooutput=1;
	assertx(argc>1);
	int nperbnd=atoi(argv[1]);
	assertx(nperbnd>0);
	int noutput=0;
	ForMeshEdge(mesh,e) {
		if (!mesh.isBoundary(e)) continue;
		const Point& p1=mesh.point(mesh.vertex1(e));
		const Point& p2=mesh.point(mesh.vertex2(e));
		for (int i=0;i<nperbnd;i++)
			outputpoint(interp(p1,p2,(i+.5)/nperbnd));
		noutput+=nperbnd;
	} EndFor;
	SHOWDF("Printed %d points on boundary edges\n",noutput);
	return 2;
}

//*** helper

static void outputedge(Edge e)
{
	static A3dElem el;
	el.init(A3dElem::TPolyline);
	el+=A3dVertex(mesh.point(mesh.vertex1(e)),Vector(0,0,0),
		      A3dVertexColor::White);
	el+=A3dVertex(mesh.point(mesh.vertex2(e)),Vector(0,0,0),
		      A3dVertexColor::White);
	oa3d.write(el);
}

static void outputface(Face f)
{
	static A3dElem el;
	el.init(A3dElem::TPolygon);
	ForFaceVertex(mesh,f,v) {
		el+=A3dVertex(mesh.point(v),Vector(0,0,0),SURFCOL);
	} EndFor;
	oa3d.write(el);
}

static void outputmesh()
{
	mesh.write(oa3d,SURFCOL);
}
