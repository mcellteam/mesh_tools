// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#ifndef MeshOp_h
#define MeshOp_h

#include "GMesh.h"
#include "Stat.h"
#include "Set.h"
#include "Queue.h"
#include "Map.h"

//*** Misc

extern void GatherBoundary(const Mesh& mesh, Edge e, Queue<Edge>& queuee);
extern void GatherComponent(const Mesh& mesh, Face f, Set<Face>& setf);

extern void MeshStatComponents(const Mesh& mesh, Stat& Scompf);
extern void MeshStatBoundaries(const Mesh& mesh, Stat& Sbound);

extern const char* MeshGenusString(const Mesh& mesh);

// Introduce centroid vertex and split face.
extern Vertex CenterSplitFace(GMesh& mesh, Face f);

// For faces with >3 sides, find a good triangulation of the vertices.
// Return: success (may fail if some edges already exist).
extern int TriangulateFace(GMesh& mesh, Face f);

extern double EdgeDihedralAngleCos(const GMesh& mesh, Edge e);

// must be a nice interior vertex
extern double VertexSolidAngle(const GMesh& mesh, Vertex v);

// Return a criterion given an edge e that is low if the edge should be
// collapsed.  Use the product of the edge length with the smallest inscribed
// radius of the two adjacent faces (its dimension is area).
extern double CollapseEdgeCriterion(const GMesh& mesh, Edge e);

//*** Retriangulate

typedef int (*EDGEF)(const GMesh& m, Edge e);

// For all Mesh Edge e,
//  if dihedral angle cos of faces both before and after is >mincos,
//   and if (fdoswap(e)) then
//    call fdel(e), swap the edge, and call fadd(newedge).
// Consider all affect edges again.  Return number of edges swapped.
extern int RetriangulateAll(GMesh& mesh, double mincos,
			    EDGEF fdoswap, EDGEF fdel=0, EDGEF fadd=0);

// Consider only Edge e and recursively, all affected edges.
// e cannot be boundary edge!
extern int RetriangulateFromEdge(GMesh& mesh, Edge e, double mincos,
				 EDGEF fdoswap, EDGEF fdel=0, EDGEF fadd=0);

// Consider swapping Edge e.  Return new_edge or 0 if not swapped.
// e cannot be boundary edge!
extern int RetriangulateOneEdge(GMesh& mesh, Edge e, double mincos,
				EDGEF fdoswap, EDGEF fdel=0, EDGEF fadd=0);

extern int CircumRadiusSwapCrit(const GMesh& mesh, Edge e);
extern int DiagonalDistanceSwapCrit(const GMesh& mesh, Edge e);

//*** Normal estimation

class Vnors {
  public:
	Vnors();
	~Vnors();
	const Vector& getnor(Face f) const;
	void clear();
  private:
	friend int ComputeVnors(const GMesh& mesh, Vertex v,
				Vnors& vnors, int subdivnormal);
	Map<Face,Vector*>* mfnor; // if zero, common normal stored in nor
	Vector nor;
	DISABLECOPY(Vnors);
};

// Compute the normal(s) of the faces around vertex v.
// Overwrite given vnors with the information.
// Normals are obtained from the limit surface of a variant of Loop
// subdivision scheme.  If subdivnormal==0, use SloanNormals instead.
// Returns the number of unique normals computed.
extern int ComputeVnors(const GMesh& mesh, Vertex v, Vnors& vnors,
			int subdivnormal=1);

//*** Projection onto mesh

// If fast!=0 and point p projects within interior of face and edges of face
// are not sharp, do not consider neighboring faces.
extern void
ProjectPNeighb(const GMesh& mesh, const Point& p, Face& pf,
	       double* pdis2, Bary* pbary, Point* pclp, int fast);

//----------------------------------------------------------------------

inline const Vector& Vnors::getnor(Face f) const
{ return mfnor?*mfnor->get(f):nor; }

#endif
