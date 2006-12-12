// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "HashFloat.h"

#include <stdint.h>

/* Idea:
   Each bucket will contain a unique representative floating point number.
   When entering a number into an empty bucket, its neighbors are tested
   first; if one of its neighbors already has a FP number, the new bucket
   inherits this number.
   -> equivalence relation is correct.
   Small numbers (fabs(x)<threshold) are placed in a special bucket that is
   adjacent to numbers that are almost small.
   */

const Univ SMALLKEY=Conv<int>::e(1);
const double SMALLVAL=1e-30;

// sqrt(2)-1, it looks good.
inline double factor(int n) { return 1+pow(.5,23-n)*.414; }

HashFloat::HashFloat(int pnignorebits, double psmall)
: nignorebits(pnignorebits), small(psmall)
{
	factor=::factor(nignorebits);
}

HashFloat::~HashFloat() { }

// Return a key that encodes the bucket in which the value lies.
inline Univ HashFloat::encode(double f) const
{
	if (fabs(f)<=small) {
		return SMALLKEY;
	} else {
		uint64_t v = floattoi64(f);
		return Conv<uint64_t>::e(assertv(v>>nignorebits));
	}
}

double HashFloat::enter(double f)
{
	int foundexact=0;
	Univ bucketn=encode(f);
	double r=m.retrieve(bucketn); // retrieve closest double
	if (r) foundexact=1;
	if (!r) r=m.retrieve(encode(f*factor));
	if (!r) r=m.retrieve(encode(f/factor));
	if (r) {		// found
		// if found in adjacent cell, propagate close value here
		if (!foundexact) m.enter(bucketn,r);
		if (r==SMALLVAL) r=0;
		return r;
	}
	double fe=f;
	if (bucketn==SMALLKEY) fe=SMALLVAL,f=0;
	m.enter(bucketn,fe);
	return f;
}
