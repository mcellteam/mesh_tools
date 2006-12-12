// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Graph_h
#define Graph_h

#include "Map.h"
#include "Stack.h"

#if 0
{
	Graph<int> g;
	for (GraphVertices<int> gv(g);gv;gv.next())
		for (GraphEdges<int> ge(g,gv());ge;ge.next())
			doedge(gv(),ge());
}
#endif

// A BGraph represents a set of tuples.
// Each tuple consists of (Univ, Univ).
// The BGraph allows quick iteration over outgoing edges.
// Duplicate edges are not allowed.
// Graph should be used to encode a relation from a set of elements to itself.
class BGraph {
  public:
	BGraph();
	~BGraph();
	void clear();
	// enter and remove domain vertices
	void enter(Univ v);	// must be new
	int contains(Univ v) const;
	int remove(Univ v);	// must have 0 outdegree, ret: was_there
	// enter an edge
	void enter(Univ v1, Univ v2); // must be new (v1 must be present)
	// enter undirected edge
	void enteru(Univ v1, Univ v2); // must be new (v1 & v2 must be present)
	int contains(Univ v1, Univ v2) const; // O(n) slow
	int remove(Univ v1, Univ v2);	      // O(n) , ret: was_there
	int removeu(Univ v1, Univ v2);	      // O(n) , ret: was_there
	int outdegree(Univ v) const;
	void add(const BGraph& g);
  private:
    friend class BGraphVertices;
    friend class BGraphEdges;
	Map<Univ,Stack<Univ>*> m;
	DISABLECOPY(BGraph);
};

// Iterate over the first vertices of the edges
class BGraphVertices {
  public:
	BGraphVertices(const BGraph& g);
	~BGraphVertices() { }
	operator bool() const;
	void next();
	Univ operator()() const;
  private:
	MapIter<Univ,Stack<Univ>*> mi;
	DISABLECOPY(BGraphVertices);
};

// Given a vertex, iterate over all outgoing edges
class BGraphEdges : public StackIter<Univ> {
  public:
	BGraphEdges(const BGraph& pg, Univ v);
	~BGraphEdges() { }
};

//----------------------------------------------------------------------------

inline BGraphVertices::BGraphVertices(const BGraph& pg) : mi(pg.m) { }
inline BGraphVertices::operator bool() const { return mi; }
inline void BGraphVertices::next() { mi.next(); }
inline Univ BGraphVertices::operator()() const { return mi.key(); }

inline BGraphEdges::BGraphEdges(const BGraph& pg, Univ v)
: StackIter<Univ>(*pg.m.get(v)) { }

//----------------------------------------------------------------------------

template<class T>
class Graph : public BGraph {
  public:
	Graph() { }
	~Graph() { }
	inline void enter(T v)
	{ BGraph::enter(Conv<T>::e(v)); }
	inline int contains(T v) const
	{ return BGraph::contains(Conv<T>::e(v)); }
	inline int remove(T v)
	{ return BGraph::remove(Conv<T>::e(v)); }
	inline void enter(T v1, T v2)
	{ BGraph::enter(Conv<T>::e(v1),Conv<T>::e(v2)); }
	inline void enteru(T v1, T v2)
	{ BGraph::enteru(Conv<T>::e(v1),Conv<T>::e(v2)); }
	inline int contains(T v1, T v2) const
	{ return BGraph::contains(Conv<T>::e(v1),Conv<T>::e(v2)); }
	inline int remove(T v1, T v2)
	{ return BGraph::remove(Conv<T>::e(v1),Conv<T>::e(v2)); }
	inline int removeu(T v1, T v2)
	{ return BGraph::removeu(Conv<T>::e(v1),Conv<T>::e(v2)); }
	inline int outdegree(T v) const
	{ return BGraph::outdegree(Conv<T>::e(v)); }
	inline void add(const Graph<T>& g) { BGraph::add(g); }
};

template<class T>
class GraphVertices : public BGraphVertices {
  public:
	GraphVertices(const Graph<T>& g) : BGraphVertices(g) { }
	~GraphVertices() { }
	inline T operator()() const
	{ return Conv<T>::d(BGraphVertices::operator()()); }
};

template<class T>
class GraphEdges : public BGraphEdges {
  public:
	GraphEdges(const Graph<T>& g, T v) : BGraphEdges(g,Conv<T>::e(v)) { }
	~GraphEdges() { }
	inline T next() { return Conv<T>::d(BGraphEdges::next()); }
	inline T operator()() const
	{return Conv<T>::d(BGraphEdges::operator()());}
};

#endif
