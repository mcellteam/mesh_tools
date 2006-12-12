// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Stack.h"

ALLOCATEPOOL(BStack);
ALLOCATEPOOL(BStack::Node);

const BStack BStack::EMPTY;

BStack::~BStack()
{
	clear();
}

void BStack::clear()
{
	Node* last;
	for (Node* n=itop;n;n=last) {
		last=n->next;
		delete n;
	}
	itop=0;
}

int BStack::height() const
{
	int num=0;
	for (BStackIter si(*this);si;si.next()) num++;
	return num;
}

int BStack::contains(Univ e) const
{
	for (BStackIter si(*this);si;si.next())
		if (si()==e) return 1;
	return 0;
}

int BStack::remove(Univ e)
{
	Node* n;
	Node* last=0;
//	for (Node* n=itop;n&&n->data!=e;n=n->next) last=n;
	for (n=itop;n&&n->data!=e;n=n->next) last=n;
	if (n) {
		if (!last) itop=n->next;
		else last->next=n->next;
		delete n;
		return 1;
	} else { return 0; }
}
