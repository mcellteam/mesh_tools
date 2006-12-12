// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Meshfit.h"
#include "LLS.h"

void GlobalFit(GMesh& wmesh, const DataPts& pt,
	       double spring, double spbf)
{
	const GMesh& mesh=wmesh;
	Map<Vertex,int> mvi;
	Array<Vertex> va;
	ForMeshVertex(mesh,v) {
		mvi.enter(v,va.num());
		va+=v;
	} EndFor;
	int m=pt.n, n=mesh.numVertices();
	if (spring) m+=mesh.numEdges();
	SHOWF("GlobalFit: about to solve a %dx%d LLS system\n",m,n);
	LLS* lls=LLS::create(m,n,3,3./n);
	// Add point constraints
	For (int i=0;i<pt.n;i++) {
		Face cmf=assertv(pt.cmf[i]);
		Vertex va[3]; mesh.vertices(cmf,va);
		for (int j=0;j<3;j++)
			lls->enter(i,mvi.get(va[j]),pt.bary[i][j]);
		lls->enterRh(i, static_cast<double const *>(pt.co[i]));
	} EndFor;
	// Add spring constraints
	if (spring) {
		double sqrtit=mysqrt(spring), sqrtbt=mysqrt(spring*spbf);
		Vector vzero(0,0,0);
		int ri=pt.n;
		ForMeshEdge(mesh,e) {
			double sqrttension=mesh.isBoundary(e)?sqrtbt:sqrtit;
			lls->enter(ri,mvi.get(mesh.vertex1(e)),sqrttension);
			lls->enter(ri,mvi.get(mesh.vertex2(e)),-sqrttension);
			lls->enterRh(ri,&vzero[0]);
			ri++;
		} EndFor;
		assertx(ri==m);
	}
	// Suggest current solution
	For (int i=0;i<n;i++) {
		lls->enterEstimate(i, static_cast<double const *>(mesh.point(va[i])));
	} EndFor;
	// Solve
	lls->solve();
	// Update solution
	For (int i=0;i<n;i++) {
		Point p; lls->getValue(i,&p[0]);
		wmesh.setPoint(va[i],p);
	} EndFor;
	delete lls,lls=0;
}
