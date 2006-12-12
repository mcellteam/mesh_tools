// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef UnionFind_h
#define UnionFind_h

#include "Map.h"

class BUnionFind {
  public:
	BUnionFind();
	~BUnionFind();
	void clear();
	// elements cannot be 0
	int unify(Univ e1, Univ e2); // ret: were_different
	int equal(Univ e1, Univ e2);
  private:
	Map<Univ,Univ> m;
	Univ rep(Univ e, int& present);
	DISABLECOPY(BUnionFind);
};

//----------------------------------------------------------------------------

template<class T>
class UnionFind : public BUnionFind {
  public:
	UnionFind() { }
	~UnionFind() { }
	inline int unify(T e1, T e2)
	{ return BUnionFind::unify(Conv<T>::e(e1),Conv<T>::e(e2)); }
	inline int equal(T e1, T e2)
	{ return BUnionFind::equal(Conv<T>::e(e1),Conv<T>::e(e2)); }
};

#endif
