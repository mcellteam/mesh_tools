// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef GraphOp_h
#define GraphOp_h

#include "Graph.h"
#include "Pqueue.h"
#include "Set.h"

#include "Geometry.h"		// Point pa[] requires size

class PointSpatial;
class Stat;

extern void GraphSymmetricClosure(BGraph& g);

// Return vertices in order of increasing graph distance from pvs.
// Vertex pvs itself is returned on first invocation of next().
// Graph may be directed.
class Dijkstra {
  public:
	Dijkstra(const BGraph& pg, Univ pvs,
		 double (*pfdist)(Univ v1, Univ v2));
	~Dijkstra();
	int done();
	Univ next(double& dis);	// ret vertex, or die
  private:
	const BGraph& g;
	Univ vs;
	double (*fdist)(Univ v1, Univ v2);
	HPqueue<Univ> pq;
	Set<Univ> set;
	DISABLECOPY(Dijkstra);
};

// Given a graph gnew consisting solely of vertices,
// computes the MST of undirectedg over the vertices in gnew under the cost
// metric fdist.
// Returns is_connected.
// Implementation: Kruskal's algorithm, O(e log(e))
//  Prim's algorithm is recommended when e=~n^2
extern int GraphMst(const BGraph& undirectedg,
		    double (*fdist)(Univ v1, Univ v2),
		    BGraph& gnew);

// Returns a newly allocated undirected graph that is the MST of undirectedg.
// Returns 0 if g is not connected.
extern BGraph* newGraphMst(const BGraph& undirectedg,
			   double (*fdist)(Univ v1, Univ v2));

// Returns a newly allocated undirected graph that is the minimum spanning
// tree of the full graph between the num points, where the cost metric
// between two points v1 and v2 is fdist(v1,v2).
// Implementation: Prim's algorithm, O(n^2)!!
extern Graph<int>* newGraphMst(int num, double (*fdist)(int v1, int v2));

// Same as GraphMst() but works specifically on an array of points and tries
// to do it more quickly by making use of a spatial data structure.
// Implementation: Prim's MST on series of subgraphs.
extern Graph<int>* newGraphQuickEmst(const Point pa[], int num,
				     const PointSpatial& sp);

// Return statistics about graph edge lengths.
// If undirected, edges stats are duplicated.
extern Stat* newGraphEdgeStats(const BGraph& g,
			       double (*fdist)(Univ v1, Univ v2));

// Returns a newly allocated directed graph that connects each vertex to its
// kcl closest neighbors (based on Euclidean distance).
// Consider applying GraphSymmetricClosure() !
extern Graph<int>* newGraphEKClosest(const Point pa[], int num, int kcl,
				     const PointSpatial& sp);

// Access each connected component of a graph.
// next() returns a representative vertex of each component.
class GraphComponent {
  public:
	GraphComponent(const BGraph& pg);
	~GraphComponent();
	operator bool() const;
	void next();
	Univ operator()() const;
  private:
	const BGraph& g;
	BGraphVertices gv;
	Set<Univ> set;
	DISABLECOPY(GraphComponent);
};

extern int GraphNumComponents(const BGraph& g);

//----------------------------------------------------------------------------

inline GraphComponent::GraphComponent(const BGraph& pg) : g(pg), gv(g) { }
inline GraphComponent::~GraphComponent() { }
inline GraphComponent::operator bool() const { return gv; }
inline Univ GraphComponent::operator()() const { return gv(); }

#endif
