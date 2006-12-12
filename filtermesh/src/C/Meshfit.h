// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Meshfit_h
#define Meshfit_h

#include "GMesh.h"
#include "Map.h"
#include "Set.h"
#include "Array.h"

class DataPts {
  public:
	DataPts();
	~DataPts();
	void clear();
	void enter(const Point& p);
	void OK();
	int n;
	Array<Point> co;
	Array<Face> cmf;
	Array<Point> clp;	// clp[i] only defined if cmf[i]
	Array<Bary> bary;	// bary[i] only defined if cmf[i]
	Map<Face,Set<int>*> mfpts; // Face -> Set of point indices
  private:
	DISABLECOPY(DataPts);
};

//*** Meshfit

extern int CheckSignal();

//*** MFFit

extern void GlobalFit(GMesh& mesh, const DataPts& pt,
		      double spring, double spbf);

//*** MFProject

extern void GlobalProject(const GMesh& mesh, DataPts& pt);

// Project point i onto stfaces and update pt.
extern void ProjectPointFaces(const GMesh& mesh, DataPts& pt,
			      int i, const Stack<Face>& stfaces);

// Reproject the points stpts on the faces in stfaces, then update the global
// data structure pt.
extern void ReprojectLocally(const GMesh& mesh, DataPts& pt,
			     const Stack<int>& stpts,
			     const Stack<Face>& stfaces);

// Update pt to reflect that point i now projects on newf.
extern void PointChangeFace(DataPts& pt, int i, Face newf);

//*** MFLocalFit

extern void GatherVertexRing(const GMesh& mesh, Vertex v,
			     Array<const Point*>& wa);

extern void GatherEdgeRing(const GMesh& mesh, Edge e,
			   Array<const Point*>& wa);

extern void LocalFit(const DataPts& pt, double spring, double spbf,
		     const Stack<int>& stpts,
		     const Array<const Point*>& wa, int niter,
		     Point& newp, double& prss0, double& prss1);

extern void FitRing(GMesh& mesh, DataPts& pt, double spring, double spbf,
		    double mincos, Vertex v, int niter);

extern double MinLocalDihedral(const Array<const Point*>& wa,
			      const Point& newp);

extern double MinDihedralAboutVertex(const GMesh& mesh, Vertex v);

#endif
