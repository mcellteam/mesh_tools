// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "MeshOp.h"
#include "GeomOp.h"
#include "Set.h"
#include "Array.h"
#include "HashStruct.h"
#include "Polygon.h"
#include "Homogeneous.h"
#include "Facedistance.h"	// for LBDistPTri(), ProjectPTri()

//*** helper

static const GMesh* gmesh;
// HashStruct so that hashing does not use pointer values (portable random)
static Univ hashe(const MEdge* e) {
	return Conv<int>::e(gmesh->vertexid(gmesh->vertex1(const_cast<Edge>(e)))+
			    gmesh->vertexid(gmesh->vertex2(const_cast<Edge>(e)))*76541);
}
static int cmpe(const MEdge* e1, const MEdge* e2) { return e1!=e2; }
typedef HashStruct<MEdge> SetEdge;

static int retriangulate(GMesh& mesh, SetEdge& sete, int recurse,
			 Set<Vertex>* setvr,
			 double mincos, EDGEF fdoswap, EDGEF fdel, EDGEF fadd)
{
	assertx(fdoswap);
	int neswapped=0;
	for (;!sete.empty();) {
		Edge e=sete.removeone();
		assertx(!mesh.isBoundary(e));
		const Point& p1=mesh.point(mesh.vertex1(e));
		const Point& p2=mesh.point(mesh.vertex2(e));
		const Point& po1=mesh.point(mesh.sideVertex1(e));
		const Point& po2=mesh.point(mesh.sideVertex2(e));
		if (DihedralAngleCos(p1,p2,po1,po2)<mincos) continue;
		if (DihedralAngleCos(po1,po2,p2,p1)<mincos) continue;
		if (!mesh.legalEdgeSwap(e)) continue;
		if (!fdoswap(mesh,e)) continue;
		if (fdel) fdel(mesh,e);
		ForEdgeFace(mesh,e,f) {
			ForFaceEdge(mesh,f,e) {
				sete.remove(e);
			} EndFor;
		} EndFor;
		Edge ne=mesh.swapEdge(e);
		neswapped++;
		if (fadd) fadd(mesh,ne);
		if (!recurse) continue;
		ForEdgeFace(mesh,ne,f) {
			ForFaceEdge(mesh,f,e) {
				if (e==ne || mesh.isBoundary(e) ||
				    setvr &&
				    (!setvr->contains(mesh.sideVertex1(e)) ||
				     !setvr->contains(mesh.sideVertex2(e))))
					continue;
				sete.add(e);
			} EndFor;
		} EndFor;
	}
	return neswapped;
}

//*** Misc

void GatherBoundary(const Mesh& mesh, Edge e, Queue<Edge>& queuee)
{
	Edge ef=e;
	for (;;) {
		queuee.enqueue(e);
		e=mesh.ccwBoundary(e);
		if (e==ef) break;
	}
}

void GatherComponent(const Mesh& mesh, Face f, Set<Face>& setf)
{
	assertx(setf.empty());
	Queue<Face> queue;
	setf.enter(f);
	for (;;) {
		ForFaceFace(mesh,f,f2) {
			if (setf.add(f2)) queue.enqueue(f2);
		} EndFor;
		if (queue.empty()) break;
		f=queue.dequeue();
	}
}

void MeshStatBoundaries(const Mesh& mesh, Stat& Sbound)
{
	Set<Edge> setevis; // boundary edges already considered
	ForMeshEdge(mesh,e) {
		if (!mesh.isBoundary(e)) continue;
		if (setevis.contains(e)) continue;
		Queue<Edge> queuee;
		GatherBoundary(mesh,e,queuee);
		ForQueue(queuee,Edge,e) {
			setevis.enter(e);
		} EndFor;
		Sbound+=queuee.length();
	} EndFor;
}

void MeshStatComponents(const Mesh& mesh, Stat& Scompf)
{
	Set<Face> setfvis; // faces already considered
	ForMeshFace(mesh,f) {
		if (setfvis.contains(f)) continue;
		Set<Face> setf;
		GatherComponent(mesh,f,setf);
		ForSet(setf,Face,ff) {
			setfvis.enter(ff);
		} EndFor;
		Scompf+=setf.num();
	} EndFor;
}

const char* MeshGenusString(const Mesh& mesh)
{
	int nv=mesh.numVertices();
	int nf=mesh.numFaces();
	int ne=mesh.numEdges();
	Stat Sbound;
	MeshStatBoundaries(mesh,Sbound);
	int nb=Sbound.num();
	Stat Scompf;
	MeshStatComponents(mesh,Scompf);
	int nc=Scompf.num();
	int ec=nv-ne+nf;	// euler characteristic
	double genus=int(nc*2-ec-nb)/2;
	int nse=0, ncv=0;
	ForMeshEdge(mesh,e) {
		if (mesh.flag(e,GMesh::SHARPE)) nse++;
	} EndFor;
	ForMeshVertex(mesh,v) {
		if (mesh.flag(v,GMesh::CUSPV)) ncv++;
	} EndFor;
	return hform("Genus: c=%d b=%d  v=%d f=%d e=%d  genus=%g%s\n",
		     nc,nb,nv,nf,ne,genus,
		     nse+ncv?hform("  sharpe=%d cuspv=%d",nse,ncv):"");
}

Vertex CenterSplitFace(GMesh& mesh, Face f)
{
	Array<Vertex> va; mesh.vertices(f,va);
	Polygon poly; mesh.polygon(f,poly);
	mesh.destroyFace(f);
	Vertex vc=mesh.createVertex();
	mesh.setPoint(vc,centroid(poly));
	for (int i=0;i<va.num();i++)
		assertx(mesh.createFace(va[i],va[(i+1)%va.num()],vc));
	return vc;
}

int TriangulateFace(GMesh& mesh, Face f)
{
        int i;
	Array<Vertex> va;
	mesh.vertices(f,va);
	int nv=va.num();
	if (assertw1(nv>3)) return 1;
	NEST {
		for (int i=2;i<nv-1;i++)
			if (assertw1(!mesh.queryEdge(va[0],va[i]))) return 0;
	}
	mesh.destroyFace(f);
//	for (int i=0;i<nv-2;i++)
	for (i=0;i<nv-2;i++)
		assertx(mesh.createFace(va[0],va[i+1],va[i+2]));
	Set<Vertex> setvr;	// vertices on ring of original face
	for (i=0;i<nv;i++) setvr.enter(va[i]);
	gmesh=&mesh;
	SetEdge sete(hashe,cmpe); // initially, inner edges
	for (i=2;i<nv-1;i++) sete.enter(mesh.edge(va[0],va[i]));
	retriangulate(mesh,sete,1,&setvr,-2,CircumRadiusSwapCrit,0,0);
	return 1;
}

// Returns cos(angle) of non-boundary edge.  Ranges from -1 to 1.
// Returns -2 if either face is degenerate.
double EdgeDihedralAngleCos(const GMesh& mesh, Edge e)
{
	Vertex v1=mesh.vertex1(e), v2=mesh.vertex2(e);
	Vertex vo1=mesh.sideVertex1(e), vo2=mesh.sideVertex2(e);
	assertx(vo2);
	return DihedralAngleCos(mesh.point(v1),mesh.point(v2),
				mesh.point(vo1),mesh.point(vo2));
}

double VertexSolidAngle(const GMesh& mesh, Vertex v)
{
	assertx(!mesh.numBoundaries(v));
	int np=mesh.degree(v), i=np;
	double ang; {		// for GNUG 2.5.8, ~SArray<Point> problem
		SArray<Point> pa(np);
		// Really want clockwise order so that solid angle points
		// toward inside of mesh.
		ForVertexCcwVertex(mesh,v,vv) {
			pa[--i]=mesh.point(vv);
		} EndFor;
		ang=SolidAngle(mesh.point(v),pa,np);
	} return ang;
}

double CollapseEdgeCriterion(const GMesh& mesh, Edge e)
{
	const Point& p1=mesh.point(mesh.vertex1(e));
	const Point& p2=mesh.point(mesh.vertex2(e));
	const Point& po1=mesh.point(mesh.sideVertex1(e));
	Vertex vo2=mesh.sideVertex2(e);
	double rc=InRadius(p1,p2,po1);
	if (vo2) {
		const Point& po2=mesh.point(vo2);
		rc=min(rc,InRadius(p1,po2,p2));
	}
	return dist(p1,p2)*rc;
}

//*** Retriangulate

int RetriangulateAll(GMesh& mesh, double mincos,
		     EDGEF fdoswap, EDGEF fdel, EDGEF fadd)
{
	gmesh=&mesh;
	SetEdge sete(hashe,cmpe);
	ForMeshEdge(mesh,e) {
		if (!mesh.isBoundary(e)) sete.enter(e);
	} EndFor;
	return retriangulate(mesh,sete,1,0,mincos,fdoswap,fdel,fadd);
}

int RetriangulateFromEdge(GMesh& mesh, Edge e, double mincos,
			  EDGEF fdoswap, EDGEF fdel, EDGEF fadd)
{
	gmesh=&mesh;
	SetEdge sete(hashe,cmpe);
	sete.enter(e);
	return retriangulate(mesh,sete,1,0,mincos,fdoswap,fdel,fadd);
}

int RetriangulateOneEdge(GMesh& mesh, Edge e, double mincos,
			 EDGEF fdoswap, EDGEF fdel, EDGEF fadd)
{
	gmesh=&mesh;
	SetEdge sete(hashe,cmpe);
	sete.enter(e);
	return retriangulate(mesh,sete,0,0,mincos,fdoswap,fdel,fadd);
}

int CircumRadiusSwapCrit(const GMesh& mesh, Edge e)
{
	assertx(!mesh.isBoundary(e));
	const Point& p1=mesh.point(mesh.vertex1(e));
	const Point& p2=mesh.point(mesh.vertex2(e));
	const Point& po1=mesh.point(mesh.sideVertex1(e));
	const Point& po2=mesh.point(mesh.sideVertex2(e));
	double rc1=CircumRadius(p1,p2,po1);
	double rc2=CircumRadius(p1,po2,p2);
	double rs1=CircumRadius(p1,po2,po1);
	double rs2=CircumRadius(p2,po1,po2);
	return max(rs1,rs2)<max(rc1,rc2);
}

int DiagonalDistanceSwapCrit(const GMesh& mesh, Edge e)
{
	assertx(!mesh.isBoundary(e));
	const Point& p1=mesh.point(mesh.vertex1(e));
	const Point& p2=mesh.point(mesh.vertex2(e));
	const Point& po1=mesh.point(mesh.sideVertex1(e));
	const Point& po2=mesh.point(mesh.sideVertex2(e));
	return dist2(po1,po2)<dist2(p1,p2);
}

//*** Normal estimation

Vnors::Vnors() : mfnor(0) { }

Vnors::~Vnors() { clear(); }

void Vnors::clear()
{
	if (mfnor) ForMapValue(*mfnor,Face,Vector*,v) { delete v; } EndFor;
	delete mfnor,mfnor=0;
}

// Cache cos() and sin() calls for low valences.
const int sizecs=12;
static double costable[sizecs+1][12];
static double sintable[sizecs+1][12];

static void initcstable()
{
	for (int i=1;i<=sizecs;i++)
		for (int j=0;j<i;j++) {
			costable[i][j]=cos(j*2*PI/i);
			sintable[i][j]=sin(j*2*PI/i);
		}
}

inline int sharp(const GMesh& mesh, Edge e)
{
	return mesh.isBoundary(e) || mesh.flag(e,GMesh::SHARPE);
}

static int extraordinarycreasev(const GMesh& mesh, Vertex v)
{
	int ne=mesh.degree(v);
	if (mesh.isBoundary(v)) return ne!=4;
	if (ne!=6) return 1;
	int nside=0, sharpef=0;
	ForVertexCcwVertex(mesh,v,vv) {
		if (sharpef==1) nside++;
		if (sharp(mesh,mesh.edge(v,vv))) sharpef++;
	} EndFor;
	assertx(sharpef==2);
	return nside!=3;
}

int ComputeVnors(const GMesh& mesh, Vertex v, Vnors& vnors, int subdivnormal)
{
        Face f;
	static int first=1; if (first) { first=0; initcstable(); }
	if (subdivnormal && assertw1(mesh.isNice(v))) subdivnormal=0;
	const Point& vp=mesh.point(v);
	int ncomp=0, nfaces=0, nsharpe=0;
	int iscusp=mesh.flag(v,GMesh::CUSPV);
	Set<Face> setfvis;
	ForVertexFace(mesh,v,f) {
		nfaces++;
		if (vnors.mfnor && !vnors.mfnor->contains(f))
			delete vnors.mfnor,vnors.mfnor=0;
	} EndFor;
	if (vnors.mfnor && vnors.mfnor->num()!=nfaces)
		delete vnors.mfnor,vnors.mfnor=0;
	ForVertexEdge(mesh,v,e) {
		if (sharp(mesh,e)) nsharpe++;
	} EndFor;
	ForVertexFace(mesh,v,frep) {
		if (setfvis.contains(frep)) continue;
		ncomp++;	// new component
		int closed=0;
		static Array<Vertex> av; av.init(0);
		static Array<Face> af; af.init(0);
//		for (Face f=frep;;) { // find f: mostClw, or frep if closed
		for (f=frep;;) { // find f: mostClw, or frep if closed
			Edge e=mesh.ccwEdge(v,f); if (sharp(mesh,e)) break;
			f=mesh.oppFace(f,e); if (f==frep) { closed=1; break; }
		}
		for (;;) {
			av+=mesh.oppVertex(v,mesh.ccwEdge(v,f));
			af+=f;
			Edge e=mesh.clwEdge(v,f);
			if (iscusp && nsharpe<2) {
				closed=0;
				av+=mesh.oppVertex(v,e);
				break;
			}
			if (!closed && sharp(mesh,e)) {
				if (mesh.oppVertex(v,e)!=av[0])
					av+=mesh.oppVertex(v,e);
				break;
			}
			f=mesh.oppFace(f,e); if (closed && f==frep) break;
		}
		int avn=av.num();
		Vector vec(0,0,0);
		if (!subdivnormal) {
			static Polygon poly; poly.init(0);
			for (int i=0;i<af.num();i++) {
				mesh.polygon(af[i],poly);
				double area=poly.getarea();
				if (assertw1(area)) continue;
				vec+=poly.getnormal()/square(area);
			}
		} else if (nsharpe>2 || iscusp) { // corner
			vec=cross(vp,mesh.point(av[0]),mesh.point(av[avn-1]));
			if (avn>2) { // from Polygon::getnormaldir()
				Vector pnor(0,0,0);
				for (int i=0;i<avn-1;i++)
					pnor+=cross(vp,mesh.point(av[i]),
						    mesh.point(av[i+1]));
				if (dot(vec,pnor)<0) vec*=-1;
			}
		} else if (nsharpe==2 && !extraordinarycreasev(mesh,v)) {
			// regular crease vertex 1,.5,-1,-1,.5
			assertx(avn==4);
			Vector v1=mesh.point(av[3])-mesh.point(av[0]);
			Vector v2=toVector(vp+mesh.point(av[0])*.5-
					   mesh.point(av[1])-mesh.point(av[2])+
					   mesh.point(av[3])*.5);
			vec=cross(v1,v2);
		} else if (nsharpe==2) { // non-regular crease vertex
			Vector v1=mesh.point(av[avn-1])-mesh.point(av[0]);
			Homogeneous h;
			if (avn==2) { // 2,-1,-1
				h=vp*2-mesh.point(av[0])-mesh.point(av[1]);
			} else if (avn==3) { // 1,0,-1,0
				h=vp-mesh.point(av[1]);
			} else {
				double theta=PI/(avn-1.);
				double w1=1./(2.-2.*cos(theta));
				h=(mesh.point(av[0])+mesh.point(av[avn-1]))*w1;
				for (int i=1;i<avn-1;i++)
					h-=mesh.point(av[i])*sin(i*theta)/
						sin(theta);
			}
			vec=cross(v1,toVector(h));
		} else {	// interior or dart
			Homogeneous h1, h2;
			for (int i=0;i<avn;i++) {
				double vcos=(avn<=sizecs?costable[avn][i]:
					    cos(i*2*PI/avn));
				double vsin=(avn<=sizecs?sintable[avn][i]:
					    sin(i*2*PI/avn));
				const Point& p=mesh.point(av[i]);
				h1+=p*vcos; h2+=p*vsin;
			}
			vec=cross(toVector(h1),toVector(h2));
		}
		assertw1(vec.normalize());
		if (af.num()==nfaces) {
			if (vnors.mfnor) delete vnors.mfnor,vnors.mfnor=0;
			vnors.nor=vec;
			break; // quick end
		}
		if (!vnors.mfnor) {
			vnors.mfnor=new Map<Face,Vector*>;
			ForVertexFace(mesh,v,f) {
				vnors.mfnor->enter(f,new Vector);
			} EndFor;
		}
		for (int i=0;i<af.num();i++) {
			setfvis.enter(af[i]);
			*assertv(vnors.mfnor->get(af[i]))=vec;
		}
	} EndFor;
	return ncomp;
}

//*** Project Point near face

void
ProjectPNeighb(const GMesh& mesh, const Point& p, Face& pf,
	       double* pdis2, Bary* pbary, Point* pclp, int fast)
{
	static int slowproj=GetenvValue("SLOWPROJECT");
        int ni;
	if (slowproj) fast=0;
	int pfsmooth=fast;
	ForFaceEdge(mesh,pf,e) {
		if (mesh.flag(e,GMesh::SHARPE)) pfsmooth=0;
	} EndFor;
	const double bnearedge=.08;
	Set<Face> setfvis;
	Point pa[3]; mesh.points(pf,pa);
	double mind2; Bary minbary;
	ProjectPTri(p,pa[0],pa[1],pa[2],&mind2,0,&minbary,pclp);
	double nearestedge=min(min(minbary[0],minbary[1]),minbary[2]);
	assertx(nearestedge>=0 && nearestedge<.34); // optional
	int nearedge=nearestedge<bnearedge;
	int projquick=pfsmooth && !nearedge; SSTAT(Sprojquick,projquick);
	if (projquick) {
		if (pdis2) *pdis2=mind2;
		if (pbary) *pbary=minbary;
		return;
	}
	int nvis=1; setfvis.enter(pf);
//	for (int ni=0;;ni++) {
	for (ni=0;;ni++) {
		// Look at faces adjacent to vertices with minbary>bnearvertex
		Set<Face> setf;
		Vertex va[3]; mesh.vertices(pf,va);
		ForIndex(j,3) {
			if (pfsmooth && minbary[(j+1)%3]>bnearedge &&
			    minbary[(j+2)%3]>bnearedge) continue;
			ForVertexFace(mesh,va[j],f) {
				if (!setfvis.contains(f)) setf.add(f);
			} EndFor;
		} EndFor;
		static Pqueue<Face> pq; pq.clear();
		ForSet(setf,Face,f) {
			Point pa[3];
			mesh.points(f,pa);
			double d2=square(LBDistPTri(p,pa[0],pa[1],pa[2]));
			pq.enterUnsorted(f,d2);
		} EndFor;
		nvis+=pq.num();
		pq.sort();
		for (;!pq.empty();) {
			if (pq.minpriority()>=mind2) break;
			Face f=pq.removemin();
			Point pa[3];
			mesh.points(f,pa);
			double d2; Bary bary; Point clp;
			ProjectPTri(p,pa[0],pa[1],pa[2],&d2,0,&bary,&clp);
			if (d2>=mind2) continue;
			pf=f; mind2=d2; minbary=bary;
			if (pclp) *pclp=clp;
			pfsmooth=fast;
			ForFaceEdge(mesh,pf,e) {
				if (mesh.flag(e,GMesh::SHARPE)) pfsmooth=0;
			} EndFor;
		}
		if (setfvis.contains(pf)) break;
		ForSet(setf,Face,f) { setfvis.enter(f); } EndFor;
	}
	SSTAT(Sprojnei,ni); SSTAT(Sprojf,nvis);
	if (ni>0) { SSTAT(Sprojunexp,!nearedge); }
	if (pdis2) *pdis2=mind2;
	if (pbary) *pbary=minbary;
}
