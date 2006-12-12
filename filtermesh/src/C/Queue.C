// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Queue.h"

ALLOCATEPOOL(BQueue::Node);

const BQueue BQueue::EMPTY;

BQueue::~BQueue()
{
	clear();
}

void BQueue::clear()
{
	Node* last;
	for (Node* n=ifront;n;n=last) {
		last=n->next;
		delete n;
	}
	ifront=irear=0;
}

int BQueue::length() const
{
	int num=0;
	for (BQueueIter qi(*this);qi;qi.next()) num++;
	return num;
}

int BQueue::contains(Univ e) const
{
	for (BQueueIter qi(*this);qi;qi.next())
		if (qi()==e) return 1;
	return 0;
}

void BQueue::addtoend(BQueue& q)
{
	if (irear) irear->next=q.ifront;
	else ifront=q.ifront;
	irear=q.irear;
	q.ifront=q.irear=0;
}

void BQueue::addtofront(BQueue& q)
{
	q.addtoend(*this);
	addtoend(q);
}
