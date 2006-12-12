// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Spatial.h"
#include "Bbox.h"
#include "Queue.h"
#include "Stack.h"
#include "Stat.h"

const int MAXGN=1290;		// cuberoot(2^31)

#define for3Dindex(ci,bi) \
for (ci[0]=bi[0][0];ci[0]<=bi[1][0];ci[0]++) \
for (ci[1]=bi[0][1];ci[1]<=bi[1][1];ci[1]++) \
for (ci[2]=bi[0][2];ci[2]<=bi[1][2];ci[2]++)


inline int Spatial::indexinbounds(int i) const
{
	return i>=0 && i<gn;
}

inline int Spatial::indicesinbounds(int ci[3]) const
{
	return indexinbounds(ci[0]) && indexinbounds(ci[1]) &&
		indexinbounds(ci[2]);
}

inline int Spatial::float2index(double fd) const
{
	double f=fd;
	if (f<0) assertx(f>-.01),f=0;
	if (f>=.99999) assertx(f<1.01),f=.99999;
	return int(f*gn);
}

inline double Spatial::index2float(int i) const
{
	return i*gni;
}

inline void Spatial::point2indices(const Point& p, int ci[3]) const
{
	for (int c=0;c<3;c++) ci[c]=float2index(p[c]);
}

inline void Spatial::indices2point(int ci[3], Point& p) const
{
	for (int c=0;c<3;c++) p[c]=index2float(ci[c]);
}

inline void Spatial::indices2bbox(int ci[3], Bbox& bb) const
{
	indices2point(ci,bb[0]);
	bb[1]=bb[0]+Vector(gni,gni,gni);
}

inline int Spatial::encode(int ci[3]) const
{
	return (ci[0]*MAXGN+ci[1])*MAXGN+ci[2];
}

inline void Spatial::decode(int en, int ci[3]) const
{
	ci[2]=en%MAXGN; en/=MAXGN;
	ci[1]=en%MAXGN; en/=MAXGN;
	ci[0]=en;
}

Spatial::Spatial(int pgn) : gn(pgn)
{
	assertx(gn<=MAXGN);
	gni=1/double(gn);
}

Spatial::~Spatial() { }

void Spatial::pqrefine(Pqueue<Univ>&, const Point&) const { }

//*** PointSpatial

PointSpatial::PointSpatial(int pgn) : Spatial(pgn) { }

PointSpatial::~PointSpatial()
{
	clear();
}

void PointSpatial::clear() {
	ForMapValue(map,int,Stack<Node*>*,stack) {
		SSTAT(Spspcelln,stack->height());
		while (!stack->empty()) delete stack->pop();
		delete stack;
	} EndFor;
}

void PointSpatial::enter(Univ id, const Point* pp)
{
	int ci[3]; point2indices(*pp,ci); assertx(indicesinbounds(ci));
	int en=encode(ci);
	Stack<PointSpatial::Node*>* cell=map.retrieve(en);
	if (!cell) map.enter(en,cell=new Stack<PointSpatial::Node*>);
	cell->push(new PointSpatial::Node(id,pp));
}

void PointSpatial::remove(Univ id, const Point* pp)
{
	int ci[3]; point2indices(*pp,ci); assertx(indicesinbounds(ci));
	int en=encode(ci);
	Stack<Node*>* cellold=map.get(en);
	Stack<Node*>* cellnew=0;
	int nfound=0;
	while (!cellold->empty()) {
		Node* e=cellold->pop();
		if (e->id==id) {
			nfound++;
			delete e;
		} else {
			if (!cellnew) cellnew=new Stack<Node*>;
			cellnew->push(e);
		}
	}
	assertx(nfound==1);
	delete map.remove(en);
	if (cellnew) map.enter(en,cellnew);
}

void PointSpatial::addcell(int ci[3], Pqueue<Univ>& pq, const Point& pcenter,
			   Set<Univ>&) const
{
	int en=encode(ci);
	Stack<PointSpatial::Node*>* cell=map.retrieve(en);
	if (!cell) return;
	ForStack(*cell,Node*,e) {
		pq.enter(Conv<PointSpatial::Node*>::e(e),dist2(pcenter,*e->p));
	} EndFor;
}

Univ PointSpatial::pqid(Univ id) const
{
	Node* e=Conv<Node*>::d(id);
	return e->id;
}

//*** ObjectSpatial

ObjectSpatial::ObjectSpatial(int pgn, DISTF papproxf2, DISTF pexactf2)
: Spatial(pgn), approxf2(papproxf2), exactf2(pexactf2) { }

ObjectSpatial::~ObjectSpatial()
{
	clear();
}

void ObjectSpatial::clear()
{
	ForMapValue(map,int,Stack<Univ>*,stack) {
		SSTAT(Sospcelln,stack->height());
		delete stack;
	} EndFor;
}

void ObjectSpatial::enter(Univ id, const Point& startp,
			  int (*fcontains)(const Bbox& bb))
{
	Set<int> set;
	Queue<int> queue;
	int ncubes=0;
	int ci[3]; point2indices(startp,ci); assertx(indicesinbounds(ci));
	int enf=encode(ci);
	set.enter(enf);
	queue.enqueue(enf);
	while (!queue.empty()) {
		int en=queue.dequeue();
		decode(en,ci);
		Bbox bb;
		indices2bbox(ci,bb);
		if (!fcontains(bb)) {
			if (en!=enf) continue; // for numerics, en==enf special
		} else {
			Stack<Univ>* cell=map.retrieve(en);
			if (!cell) map.enter(en,cell=new Stack<Univ>);
			cell->push(id);
			ncubes++;
		}
		int bi[2][3];
		for (int c=0;c<3;c++)
			bi[0][c]=max(ci[c]-1,0),bi[1][c]=min(ci[c]+1,gn-1);
		for3Dindex(ci,bi) {
			int en=encode(ci);
			if (set.add(en)) queue.enqueue(en);
		}
	}
	SSTAT(Sospobcells,ncubes);
}

void ObjectSpatial::searchsegment(const Point& p1, const Point& p2,
				  int (*ftest)(Univ id)) const
{
	Set<Univ> set;
	int shouldstop=0;
	double maxe=0;
	for (int c=0;c<3;c++) {
		assertx(p1[c]>=0 && p1[c]<=1);
		assertx(p2[c]>=0 && p2[c]<=1);
		maxe=max(maxe, fabs(p2[c]-p1[c]));
	}
	int ni=float2index(maxe)+2; // 2 there just to be safe
	Vector v=(p2-p1)/ni;
	Point p=p1;
	int pci[3]; point2indices(p,pci);
	int pen=0;
	for (int i=0;;i++) {
		int cci[3]; point2indices(p,cci);
		assertx(indicesinbounds(cci)); // optional
		int bi[2][3];
		bi[0][0]=min(cci[0],pci[0]); bi[1][0]=max(cci[0],pci[0]);
		bi[0][1]=min(cci[1],pci[1]); bi[1][1]=max(cci[1],pci[1]);
		bi[0][2]=min(cci[2],pci[2]); bi[1][2]=max(cci[2],pci[2]);
		int ci[3];
		for3Dindex(ci,bi) {
			int en=encode(ci);
			if (en==pen) continue;
			Stack<Univ>* cell=map.retrieve(en);
			if (!cell) continue;
			ForStack(*cell,Univ,e) {
				if (set.add(e) && ftest(e)) shouldstop=1;
			} EndFor;
		}
		if (i==ni || shouldstop) break;
		pci[0]=cci[0]; pci[1]=cci[1]; pci[2]=cci[2];
		pen=encode(pci);
		p+=v;
	}
	if (!shouldstop) assertw(!compare(p,p2,1e-5));
}

void ObjectSpatial::addcell(int ci[3], Pqueue<Univ>& pq, const Point& pcenter,
			    Set<Univ>& set) const
{
	int en=encode(ci);
	Stack<Univ>* cell=map.retrieve(en);
	if (!cell) return;
	ForStack(*cell,Univ,e) {
		if (set.add(e)) pq.enter(e,approxf2(pcenter,e));
	} EndFor;
}

void ObjectSpatial::pqrefine(Pqueue<Univ>& pq, const Point& pcenter) const
{
	Univ id=pq.min();
	double oldv=pq.minpriority();
	double newv=exactf2(pcenter,id);
	if (newv==oldv) return;
	if (newv<oldv-1e-12 && Warning("newv<oldv"))
		SHOWN(oldv),SHOWN(newv);
	assertx(pq.removemin()==id);
	pq.enter(id,newv);
}

Univ ObjectSpatial::pqid(Univ id) const
{
	return id;
}

//*** SpatialSearch

SpatialSearch::SpatialSearch(const Spatial& psp, const Point& pp,
			     double pmaxdis)
: sp(psp), pcenter(pp), maxdis(pmaxdis), disbv2(0), ncellsv(0), nelemsv(0)
{
	int ci[3];
	sp.point2indices(pcenter,ci);
	assertx(sp.indicesinbounds(ci));
	for (int i=0;i<2;i++) for (int c=0;c<3;c++) ssi[i][c]=ci[c];
	consider(ci);
	getclosestnextcell();
}

SpatialSearch::~SpatialSearch()
{
	SSTAT(Sssncellsv,ncellsv);
	SSTAT(Sssnelemsv,nelemsv);
}

int SpatialSearch::done()
{
	for (;;) {
		if (!pq.empty()) return 0;
		if (disbv2>=square(maxdis)) return 1;
		expandsearchspace();
	}
}

Univ SpatialSearch::next(double* pdis2)
{
	Univ u;
	for (;;) {
		if (pq.empty()) assertx(!done()); // refill pq
		double dis2=pq.minpriority();
		if (dis2>disbv2) {
			expandsearchspace();
			continue;
		}
		u=pq.min();
		sp.pqrefine(pq,pcenter);
		if (pq.min()!=u || pq.minpriority()!=dis2) continue;
		if (pdis2) *pdis2=pq.minpriority();
		u=pq.removemin();
		break;
	}
	return sp.pqid(u);
}

void SpatialSearch::consider(int ci[3])
{
	ncellsv++;
	int n=pq.num();
	sp.addcell(ci,pq,pcenter,setevis);
	nelemsv+=pq.num()-n;
}

void SpatialSearch::getclosestnextcell()
{
	double mindis=1e10;
	for (int c=0;c<3;c++) {
		if (ssi[0][c]>0) {
			double a=pcenter[c]-sp.index2float(ssi[0][c]);
			if (a<mindis) mindis=a,axis=c,dir=0;
		}
		if (ssi[1][c]<sp.gn-1) {
			double a=sp.index2float(ssi[1][c]+1)-pcenter[c];
			if (a<mindis) mindis=a,axis=c,dir=1;
		}
	}
	// mindis may be big if all of space has been searched
	disbv2=square(mindis);
}

void SpatialSearch::expandsearchspace()
{
	assertx(axis>=0 && axis<3 && dir>=0 && dir<=1);
	int bi[2][3];
	for (int i=0;i<2;i++) for (int c=0;c<3;c++) bi[i][c]=ssi[i][c];
	ssi[dir][axis]+=dir?1:-1;
	bi[0][axis]=bi[1][axis]=ssi[dir][axis];
	// consider the layer whose axis's value is ssi[dir][axis]
	int ci[3];
	for3Dindex(ci,bi) consider(ci);
	getclosestnextcell();
}
