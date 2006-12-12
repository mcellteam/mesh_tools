// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Meshfit.h"
#include "PolygonSpatial.h"
#include "Facedistance.h"
#include "MeshOp.h"

// Declarations
static void localproject(const GMesh& mesh, DataPts& pt);
static void globalproject(const GMesh& mesh, DataPts& pt);

void GlobalProject(const GMesh& mesh, DataPts& pt)
{
	int nzcmf=0;
	for (int i=0;i<pt.n;i++) nzcmf+=!pt.cmf[i];
	if (nzcmf==pt.n) globalproject(mesh,pt);
	else if (!nzcmf) localproject(mesh,pt);
	else assertnever("");
}

static void globalproject(const GMesh& mesh, DataPts& pt)
{
	const int nv=mesh.numVertices();
	PolygonSpatial psp(nv<10000?15:nv<30000?25:35);
	Map<Polygon*,Face> mpf;
	ForMeshFace(mesh,f) {
		Polygon* p=new Polygon(3);
		mesh.polygon(f,*p);
		psp.enter(p);
		mpf.enter(p,f);
	} EndFor;
	for (int i=0;i<pt.n;i++) {
		if (CheckSignal()) break;
		assertx(!pt.cmf[i]);
		SpatialSearch ss(psp,pt.co[i]);
		Polygon* p=Conv<Polygon*>::d(ss.next());
		Face f=mpf.get(p);
		PointChangeFace(pt,i,f);
		Point pa[3];
		mesh.points(f,pa);
		ProjectPTri(pt.co[i],pa[0],pa[1],pa[2],0,0,
			    &pt.bary[i],&pt.clp[i]);
	}
	ForMapKey(mpf,Polygon*,Face,poly) { delete poly; } EndFor;
}

static void localproject(const GMesh& mesh, DataPts& pt)
{
	for (int i=0;i<pt.n;i++) {
		if (CheckSignal()) break;
		Face cf=pt.cmf[i];
		ProjectPNeighb(mesh,pt.co[i],cf,0,&pt.bary[i],&pt.clp[i],0);
		PointChangeFace(pt,i,cf);
	}
}

void ReprojectLocally(const GMesh& mesh, DataPts& pt,
		      const Stack<int>& stpts,
		      const Stack<Face>& stfaces)
{
	ForStack(stpts,int,pi) {
		ProjectPointFaces(mesh,pt,pi,stfaces);
	} EndFor;
}

void ProjectPointFaces(const GMesh& mesh, DataPts& pt,
		       int i, const Stack<Face>& stfaces)
{
	static Pqueue<Face> pq;
	pq.clear();
	const Point& p=pt.co[i];
	ForStack(stfaces,Face,f) {
		Point pa[3];
		mesh.points(f,pa);
		double d2=square(LBDistPTri(p,pa[0],pa[1],pa[2]));
		pq.enterUnsorted(f,d2);
	} EndFor;
	pq.sort();
	double mind2=1e30;
	Face minf=0;
	for (;!pq.empty();) {
		if (pq.minpriority()>=mind2) break;
		Face f=pq.removemin();
		Point pa[3];
		mesh.points(f,pa);
		double d2; Bary b; Point clp;
		ProjectPTri(p,pa[0],pa[1],pa[2],&d2,0,&b,&clp);
		if (d2<mind2) {
			mind2=d2,minf=f;
			pt.bary[i]=b;
			pt.clp[i]=clp;
		}
	}
	assertx(minf);
	PointChangeFace(pt,i,minf);
}

void PointChangeFace(DataPts& pt, int i, Face newf)
{
	Face oldf=pt.cmf[i];
	if (newf==oldf) return;
	if (oldf) assertx(pt.mfpts.get(oldf)->remove(i));
	pt.cmf[i]=newf;
	if (newf) pt.mfpts.get(newf)->enter(i);
}
