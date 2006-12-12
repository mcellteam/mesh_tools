// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "UnionFind.h"

// union-find: with path compression but without weight-balancing
//  -> worst case O(nlogn), good case O(n)

BUnionFind::BUnionFind() { }

BUnionFind::~BUnionFind() { }

void BUnionFind::clear()
{
	m.clear();
}

int BUnionFind::unify(Univ e1, Univ e2)
{
	if (e1==e2) return 0;
	int present;
	Univ r1=rep(e1,present); if (!present) m.enter(r1=e1,e1);
	Univ r2=rep(e2,present); if (!present) m.enter(r2=e2,e2);
	if (r1==r2) return 0;
	m.replace(r1,r2);
	return 1;
}

int BUnionFind::equal(Univ e1, Univ e2)
{
	if (e1==e2) return 1;
	int present;
	Univ r1=rep(e1,present); if (!present) return 0;
	Univ r2=rep(e2,present); if (!present) return 0;
	return r1==r2;
}

Univ BUnionFind::rep(Univ e, int& present)
{
        Univ t;
	Univ r=m.retrieve(e, present);
	if (!present) return 0;	// alone, ret anything
//	for (t=e;r!=t;) r=m.get(t=r);
	for (Univ t=e;r!=t;) r=m.get(t=r);
	for (;e!=r;e=t) {
		t=m.get(e);
		m.replace(e,r);
	}
	return r;
}
