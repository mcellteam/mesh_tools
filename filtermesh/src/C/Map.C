// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Map.h"
#include "Set.h"
#include "Random.h"

ALLOCATEPOOL(BMap);
ALLOCATEPOOL(BMap::Node);

const BMap BMap::EMPTY;

BMap::BMap() : b(0), isize(0), inum(0), fbuckn(0) { }

void BMap::clear()
{
	for (int i=0;i<isize;i++) {
		for (Node* n=b[i];n;) {
			Node* last=n->n;
			delete n;
			n=last;
		}
	}
	delete[] b;
	b=0;
	isize=0;
	inum=0;
	fbuckn=0;
}

BMap::~BMap()
{
	if (b) clear();
}

void BMap::resize(int newsize)
{
	BMap hn;
	typedef BMap::Node* pNode;
	hn.b=new pNode[hn.isize=newsize];

	// This shouldn't be necessary, but pool allocator doesn't return
	// zeroed memory!
	memset(hn.b, 0, sizeof(pNode) * newsize);

	for (BMapIter mi(*this);mi;mi.next())
		hn.quickenter(mi.key(),mi.value());
	hn.inum=inum;
	clear();
	// want to simulate *this=*hn;
	b=hn.b; isize=hn.isize; inum=hn.inum; fbuckn=0;
	// pretend hn has no data so destructor does nothing
	hn.b=0; hn.isize=0; hn.inum=0;
}

Univ BMap::remove(Univ k)
{
	if (!b) return 0;
	Node* n; Node* last;
	int buckn=hashk(k);
	for (last=0,n=b[buckn];n;last=n,n=n->n)
		if (n->k==k) break;
	if (!n) return 0;
	Univ e=n->v;
	removeaux(buckn,n,last);
	return e;
}

void BMap::removeaux(int buckn, Node* n, Node* last)
{
	assertx(n);
	if (last) {
		last->n=n->n;
	} else {
		b[buckn]=n->n;
	}
	delete n;
	if (--inum<=isize/2 && isize>5) resize((isize-3)/5+1);
	if (!inum) clear();
}

void BMap::OK() const
{
        int i;
	Set<Univ> set;
	For (BMapIter mi(*this);mi;mi.next()) {
		assertx(set.add(mi.key()));
	} EndFor;
	assertx(set.num()==inum);
	assertx((b?1:0)==(isize?1:0));
	assertx((b?1:0)==(inum?1:0));
//	for (int i=0;i<isize;i++) if (b[i]) break;
	for (i=0;i<isize;i++) if (b[i]) break;
	assertx(fbuckn>=0 && fbuckn<=i);
}

//*** BMapIter

// update h.fbuckn if necessary
inline void BMapIter::findrealfbuckn()
{
	// register BMap* vthat=const_cast<BMap*>(&h);	// virtual const
	register int ofbuckn=h.fbuckn;
	bn=ofbuckn-1; n=0;
	advance();
	if (bn>ofbuckn) h.fbuckn=bn;
}

void BMapIter::advancewrap()
{
	advance();
	if (n) return;
	assertx(h.inum);
	bn=h.fbuckn-1;
	advance();
	assertx(n);
}

BMapIter::BMapIter(const BMap& hh)
: h(hh)
{
	findrealfbuckn();
}

BMapIter::BMapIter(const BMap& hh, Random& r)
: h(hh)
{
	findrealfbuckn();
	if (!h.b) return;
	bn=h.fbuckn+r.getint()%(h.isize-h.fbuckn);
	int ne=0;
	for (n=h.b[bn];n;n=n->n) ne++;
	int nskip=r.getint()%(20+ne);
	if (!(n=h.b[bn])) advancewrap();
	for (int i=0;i<nskip;i++) {
		if (!(n=n->n)) advancewrap();
	}
}
