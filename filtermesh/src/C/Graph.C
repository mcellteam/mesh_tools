// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Graph.h"

BGraph::BGraph() { }

BGraph::~BGraph()
{
	clear();
}

void BGraph::clear()
{
	ForMapValue(m,Univ,Stack<Univ>*,stack) { delete stack; } EndFor;
	m.clear();
}

void BGraph::enter(Univ v)
{
	m.enter(v,new Stack<Univ>);
}

int BGraph::contains(Univ v) const
{
	return m.contains(v);
}

int BGraph::remove(Univ v)
{
	Stack<Univ>* stack=m.remove(v);
	if (stack) { assertx(stack->empty()); delete stack; }
	return stack?1:0;
}

void BGraph::enter(Univ v1, Univ v2)
{
	m.get(v1)->push(v2);
}

void BGraph::enteru(Univ v1, Univ v2)
{
	enter(v1,v2);
	enter(v2,v1);
}

int BGraph::contains(Univ v1, Univ v2) const
{
	return m.get(v1)->contains(v2);
}

int BGraph::remove(Univ v1, Univ v2)
{
	return m.get(v1)->remove(v2);
}

int BGraph::removeu(Univ v1, Univ v2)
{
	int r1=remove(v1,v2);
	int r2=remove(v2,v1);
	assertx(r1==r2);
	return r1;
}

int BGraph::outdegree(Univ v1) const
{
	return m.get(v1)->height();
}

void BGraph::add(const BGraph& g)
{
	for (BGraphVertices gv(g);gv;gv.next())
		for (BGraphEdges ge(g,gv());ge;ge.next())
			if (!contains(gv(),ge())) enter(gv(),ge());
}
