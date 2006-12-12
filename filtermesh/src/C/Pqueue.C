// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Pqueue.h"

//*** BPqueue

BPqueue::BPqueue(int size)
: isize(size?size:1), inum(0), ar(new Node[isize]) { }

BPqueue::~BPqueue() { delete[] ar; }

void BPqueue::clear() { inum=0; }

inline void BPqueue::nswitch(int n1, int n2)
{
	swap(&ar[n1],&ar[n2]);
}

void BPqueue::padjust(int n, int uptoo) // copied in BHPqueue
{
	const double cp=ar[n].pri;
	int nn;
	for (;;nswitch(n,nn),n=nn) {
		if (uptoo) {
			nn=(n-1)/2;
			if (n && cp<ar[nn].pri) continue;
		}
		int ln=n*2+1;	// left child
		if (ln>=inum) break; // no children
		double lp=ar[ln].pri;
		int rn=n*2+2;	// right child
		if (rn>=inum) { // no right child
			if (lp<cp) { nn=ln; continue; }
			break;
		}
		double rp=ar[rn].pri;
		if (lp<cp && lp<rp) { nn=ln; continue; }
		if (rp<cp) { nn=rn; continue; }
		break;
	}
}

void BPqueue::vnswitch(int n1, int n2) { nswitch(n1,n2); }
void BPqueue::vadjust(int n) { padjust(n,1); }

Univ BPqueue::removeminI()
{
	assertx(inum);
	Univ e=ar[0].e;
	if (inum>1) vnswitch(0,inum-1);	// virtual call (called from BHPqueue)
	if (--inum>1) vadjust(0);	// virtual call
	return e;
}

void BPqueue::resize(int newsize)
{
	assertx(newsize>=inum);
	Node* nar=new Node[isize=newsize];
	for (int i=0;i<inum;i++) nar[i]=ar[i];
	delete[] ar;
	ar=nar;
}

void BPqueue::sort()		// virtual for speed
{
	for (int i=(inum-2)/2;i>=0;i--) padjust(i,0);
}

//*** BHPqueue

BHPqueue::BHPqueue(int size) : BPqueue(size) { }

BHPqueue::~BHPqueue() { }

void BHPqueue::clear()
{
	BPqueue::clear();
	m.clear();
}

inline void BHPqueue::nswitch(int n1, int n2)
{
	BPqueue::nswitch(n1,n2);
	assertx(m.replace(ar[n1].e,n1)==n2);
	assertx(m.replace(ar[n2].e,n2)==n1);
}

void BHPqueue::padjust(int n, int uptoo) // copied from BPqueue
{
	const double cp=ar[n].pri;
	int nn;
	for (;;nswitch(n,nn),n=nn) {
		if (uptoo) {
			nn=(n-1)/2;
			if (n && cp<ar[nn].pri) continue;
		}
		int ln=n*2+1;	// left child
		if (ln>=inum) break; // no children
		double lp=ar[ln].pri;
		int rn=n*2+2;	// right child
		if (rn>=inum) { // no right child
			if (lp<cp) { nn=ln; continue; }
			break;
		}
		double rp=ar[rn].pri;
		if (lp<cp && lp<rp) { nn=ln; continue; }
		if (rp<cp) { nn=rn; continue; }
		break;
	}
}

inline void BHPqueue::adjust(int n) { padjust(n,1); }
void BHPqueue::vnswitch(int n1, int n2) { nswitch(n1,n2); }
void BHPqueue::vadjust(int n) { adjust(n); }

void BHPqueue::enterUnsorted(Univ e, double pri)
{
	m.enter(e,inum);
	BPqueue::enterUnsorted(e,pri);
}

Univ BHPqueue::removeminI()
{
	Univ e=BPqueue::removeminI();
	(void)m.remove(e);	// may return index 0
	return e;
}

inline int BHPqueue::find(Univ e) const
{
	int present, i=m.retrieve(e,present);
	return present?i:-1;
}

double BHPqueue::retrieve(Univ e) const
{
	int i=find(e);
	return i>=0?ar[i].pri:-1;
}

double BHPqueue::remove(Univ e)
{
	int i=find(e);
	if (i<0) return -1;
	double ppri=ar[i].pri;
	if (inum-1!=i) nswitch(inum-1,i);
	if (i<--inum) adjust(i);
	(void)m.remove(e);	// may return index 0
	if (inum<isize*.4 && isize>100) resize(isize/2);
	return ppri;
}

double BHPqueue::update(Univ e, double pri)
{
	assertx(pri>=0);
	int i=find(e);
	if (i<0) return -1;
	double oldpri=ar[i].pri;
	ar[i].pri=pri;
	adjust(i);
	return oldpri;
}

double BHPqueue::enterupdate(Univ e, double pri)
{
	if (m.contains(e)) return update(e,pri);
	return enter(e,pri),-1;
}

void BHPqueue::sort()		// virtual for speed
{
	for (int i=(inum-2)/2;i>=0;i--) padjust(i,0);
}
