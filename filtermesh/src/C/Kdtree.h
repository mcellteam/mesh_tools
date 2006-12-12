// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Kdtree_h
#define Kdtree_h

#include "Array.h"
#include "Hh.h"
#include "Stack.h"
#include "Stat.h"

#include <iostream>

class KdtreeNode;
class Stat;

extern int KDTREE_STATS;
extern Stat SKDsearchnel;

template<class T>
class Kdtree {
  public:
	Kdtree(int pnd, int pmaxlevel=8); // # of dimensions, max # subdiv.
	~Kdtree();
	void clear();
	void allowduplication(double pfsize) { if (!getenv("KDFSIZE")) fsize=pfsize; }

	void enter(T id, const double bb[][2]); // bb copied

	// search is reentrant (for Hlr).
	// Start searching at loc (default root) for objects whose bb
	// intersect the one given.
	// For each object found, call cbfunc with id, bb, and current
	// location in the tree.
	// cbfunc ret: 0=nothing, 1=bb_shrunk, 2=stop_search
	// cbfunc may modify the bb by shrinking it
	typedef int (*CBFUNC)(T id, double bb[][2], KdtreeNode* floc);
	int search(double bb[][2], CBFUNC cbfunc, KdtreeNode* loc=0) const;
				// ret was_stopped

	void print() const { recPrint(root,0); }

  private:
	const int nd;		// number of dimensions of space
	const int maxlevel;	// maximum # of subdivision on each axis
	KdtreeNode* root;
	double fsize;		// ok to duplicate if average edge<fsize

	void recClear(KdtreeNode* n);
	void recDepth(KdtreeNode* n, Stat& stat, int depth) const;
	void recEnter(T id, const double bb[][2],
		      KdtreeNode* n, SArray<double>& aval,
		      int level, double inc, int axis, double avgel);
	int recSearch(KdtreeNode* n, KdtreeNode* lca,
		      double bb[][2], CBFUNC cbfunc, int &nelemvis) const;
	void recPrint(KdtreeNode* n, int l) const;
	DISABLECOPY(Kdtree);
};

class Kdentry {
  public:
	Kdentry() { }
	Univ id;
	double bb[3][2];		// reserve space for 3-dim, allocate more
	// no members can come after bb!
// custom POOLALLOCATION
//	static Pool pool;
	static Kdentry* specialnew(int nd) {
		if (nd <= 3)
			return new Kdentry;
		else
		{
			char *data = new char[sizeof(Kdentry) + (nd-3)*2*sizeof(double)];
			new (data) Kdentry();
			return reinterpret_cast<Kdentry *>(data);
		}
	}
	static void specialdelete(void* p, int nd) {
		Kdentry *k = reinterpret_cast<Kdentry *>(p);
		if (nd <= 3)
			delete k;
		else
		{
			k->Kdentry::~Kdentry();
			delete[] reinterpret_cast<char *>(p);
		}
	}
  private:
	DISABLECOPY(Kdentry);
//	void* operator new(size_t); void operator delete(void*); // disable
};

class KdtreeNode {
  public:
	KdtreeNode() : axis(-1), l(0), h(0) { }
	Stack<Kdentry*> stack;
	int axis;
	double val;
	KdtreeNode* l;	// lower valued subtree
	KdtreeNode* h;	// higher valued subtree
  private:
	DISABLECOPY(KdtreeNode);
};

#include "Kdtree.tcc"

#endif
