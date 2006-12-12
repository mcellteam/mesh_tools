// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef HashFloat_h
#define HashFloat_h

#include "Map.h"

/* XXX: What is a good value for pnignorebits? */
class HashFloat {
  public:
	HashFloat(int pnignorebits=8, double psmall=1e-4);
	~HashFloat();
	double enter(double f);	// ret: filtered value
  private:
	Map<Univ,double> m;	// encoded double bucket -> double rep
	int nignorebits;	// # of bits to ignore in FP representation
	double small;		// numbers with abs<small are grouped at 0
	double factor;		// used to access prev and next buckets
	Univ encode(double v) const;
	DISABLECOPY(HashFloat);
};

#endif
