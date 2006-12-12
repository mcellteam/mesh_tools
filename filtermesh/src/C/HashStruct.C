// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "HashStruct.h"

BHashStruct::BHashStruct(HASHF phashf, CMPF pcmpf)
: hashf(phashf), cmpf(pcmpf) { }

BHashStruct::~BHashStruct() { assertw1(!num()); }

void BHashStruct::enter(Univ s)
{
	assertx(s);
	Univ k=hashf(s);
	BMap::enter(k,s);
}

Univ BHashStruct::retrieve(Univ s) const
{
	Node* n;
	assertx(s);
	if (!b) return 0;
	Univ k=hashf(s);
//	for (Node* n=b[hashk(k)];n;n=n->n)
	for (n=b[hashk(k)];n;n=n->n)
		if (n->k==k && !cmpf(s,n->v)) break;
	if (!n) return 0;
	return n->v;
}

int BHashStruct::add(Univ s)
{
	assertx(s);
	if (retrieve(s)) return 0;
	return enter(s),1;
}

int BHashStruct::remove(Univ s)
{
	assertx(s);
	if (!b) return 0;
	Node* n; Node* last;
	Univ k=hashf(s);
	int buckn=hashk(k);
	for (last=0,n=b[buckn];n;last=n,n=n->n)
		if (n->k==k && !cmpf(s,n->v)) break;
	if (!n) return 0;
	removeaux(buckn,n,last);
	return 1;
}
