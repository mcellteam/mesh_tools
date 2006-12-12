// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "HashPoint.h"

// ALLOCATEPOOL(HashPointEntry);
ALLOCATEPOOL(HashPoint::Entry);

Univ HashPoint::hashEntry(const Entry* e)
{
	const Point& p=e->p;
	uint64_t hv=0;
	for (int c = 0; c < 3; c++)
		hv = hv*17 + floattoi64(p[c]);
	return Conv<uint64_t>::e(hv);
}

int HashPoint::cmpEntry(const Entry* e1, const Entry* e2)
{
	return compare(e1->p,e2->p);
}

HashPoint::HashPoint(int pnignorebits, double psmall)
: hf(pnignorebits,psmall), hs(hashEntry,cmpEntry), index(0) { }

HashPoint::~HashPoint()
{
	ForHashStruct(hs,HashPoint::Entry,e) { delete e; } EndFor;
	hs.clear();		// avoid warning message
}

int HashPoint::enter(const Point& p)
{
	static Entry e;
	for (int c=0;c<3;c++) e.p[c]=hf.enter(p[c]);
	Entry* er=hs.retrieve(&e);
	if (er) return er->index;
	er=new Entry;
	er->p=e.p;
	er->index=index;
	hs.enter(er);
	return index++;
}
