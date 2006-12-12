// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Combination.h"
#include "Stack.h"

ALLOCATEPOOL(BCombination);

const BCombination BCombination::EMPTY;

BCombination::BCombination() { }

BCombination::~BCombination() { clear(); }

void BCombination::clear()
{
	m.clear();
}

void BCombination::squeeze() const
{
	BCombination* const vthis=const_cast<BCombination *>(this); // virtual const
	Stack<Univ> stack;
	ForMapKeyValue(m,Univ,e,double,v) {
		if (!v) stack.push(e);
	} EndFor;
	while (!stack.empty()) (void)vthis->m.remove(stack.pop());
}

double BCombination::sum() const
{
	double sum=0;		// could be double
	for (BCombinationIter ci(*this);ci;ci.next()) sum+=ci.value();
	return sum;
}
