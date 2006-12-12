// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef HashPoint_h
#define HashPoint_h

#include "HashFloat.h"
#include "HashStruct.h"
#include "Geometry.h"
#include "Pool.h"

class HashPoint {
  public:
	HashPoint(int pnignorebits=8, double psmall=1e-4);
	~HashPoint();
	int enter(const Point& p); // p is copied, ret: index (first is 0)
  private:
	struct Entry {		// moved in for GNUG 2.4.5
		Point p;
		int index;
		POOLALLOCATION(HashPoint::Entry);
	};
// 	typedef HashPointEntry Entry; // for GNUG 2.5.5, no need in 2.5.8
	HashFloat hf;
	HashStruct<HashPoint::Entry> hs;
	int index;		// current index, assigned to next new point
	static Univ hashEntry(const Entry* e);
	static int cmpEntry(const Entry* e1, const Entry* e2);
	DISABLECOPY(HashPoint);
};

#endif
