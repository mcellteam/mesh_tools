// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Pool.h"

// Two reasons to check size_t in new and delete:
// - The class may be derived and may not provide its own new/delete operators
// - GNUG 2.4.5 does not seem to call global operators for array allocation

// New problem in GNUG 2.5.5:
// when a destructor is not explicitly declared for a class,
// delete[] calls pool::delete, and even when !a, and even if using ::delete[]
//  (so then there is no way to deallocate what has been allocated with new[])
// problem is still present in GNUG 2.5.8, so make sure destructors are
// declared for all Pooled objects.

// static void beg_Pool()
// {
// }
// 
// static void end_Pool()
// {
// }
// 
// MODULEDEFINE(Pool);

const int malloc_overhead=16;
const int chunksize=8*1024-malloc_overhead;

static int debug=GetenvValue("DEBUG");

Pool::Pool(unsigned pesize, const char* pname)
: esize(pesize), iname(pname), h(0), chunkh(0), nalloc(0), offset(0)
{
	// make allocated size a multiple of sizeof(Link)!!
	esize=int((esize+sizeof(Link)-1)/sizeof(Link))*sizeof(Link);
	assertx(static_cast<long>(esize) >= static_cast<long>(sizeof(Link)) &&
		(esize%sizeof(Link))==0);
	offset=chunksize-sizeof(Chunk);
}

void* Pool::specialalloc(size_t s)
{
	if (debug>=2) SHOWF("Pool %-20s: alloc size %d\n",iname,s);
	return hmalloc(0,s);
}

void Pool::specialfree(void* p, size_t s)
{
	if (debug>=2) SHOWF("Pool %-20s: free  size %d\n",iname,s);
	hfree(p);
}

Pool::~Pool()
{
	const int found_fix=0;	// use CleanUp(), static Set<> in G3dGL
	int n=0;
	for (Link* p=h;p;p=p->next) n++;
	if (debug>=2 || debug && nalloc || found_fix && n!=nalloc)
		SHOWF("Pool %-20s: %6d/%-6d elements outstanding%s\n",
		      iname,nalloc-n,nalloc,n!=nalloc?" **":"");
	// The Pool may be destroyed before other static structures, so we
	// cannot really report an error here.
	if (n!=nalloc) return;
	for (Chunk* chunk=chunkh;chunk;) {
		char* p=(char*)chunk;
		chunk=chunk->next;
		delete[] (p-offset);
	}
}

void Pool::grow()
{
	assertx(!h);
	const int nelem=offset/esize;
	assertx(nelem>0);
	if (debug>=3) SHOWF("Pool %-20s: allocating new chunk\n",iname);
	char* p=new char[chunksize];
	Chunk* chunk=(Chunk*)(p+offset);
	chunk->next=chunkh;
	chunkh=chunk;
	h=(Link*)p;
	char* l=p+(nelem-1)*esize;
	for (;p<l;p+=esize)
		((Link*)p)->next=(Link*)(p+esize);
	((Link*)p)->next=0;
	nalloc+=nelem;
}
