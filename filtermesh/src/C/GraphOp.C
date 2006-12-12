// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Geometry.h"
#include "GraphOp.h"
#include "UnionFind.h"
#include "Spatial.h"
#include "Stat.h"
#include "Queue.h"
#include "Array.h"


void GraphSymmetricClosure(BGraph& g)
{
	for (BGraphVertices gv(g);gv;gv.next())
		for (BGraphEdges ge(g,gv());ge;ge.next())
			if (!g.contains(ge(),gv())) g.enter(ge(),gv());
}

//*** Dijkstra

Dijkstra::Dijkstra(const BGraph& pg, Univ pvs,
		   double (*pfdist)(Univ v1, Univ v2))
: g(pg), vs(pvs), fdist(pfdist)
{
	pq.enter(vs,0);
}

Dijkstra::~Dijkstra() { }

int Dijkstra::done() { return pq.empty(); }

Univ Dijkstra::next(double& dis)
{
	assertx(!pq.empty());
	double dmin=pq.minpriority();
	Univ vmin=pq.removemin();
	set.enter(vmin);
	for (BGraphEdges ge(g,vmin);ge;ge.next()) {
		Univ v=ge();
		if (set.contains(v)) continue;
		double pnd=dmin+fdist(vmin,v); // possibly smaller distance
		double cd=pq.retrieve(v);      // current distance
		if (cd<0) {
			pq.enter(v,pnd);
		} else if (pnd<cd) {
			pq.update(v,pnd);
		}
	}
	dis=dmin;
	return vmin;
}


//*** GraphMst KRUSKAL

struct GMSTedgenode {
	GMSTedgenode(Univ a1, Univ a2) : v1(a1), v2(a2) { }
	Univ v1;
	Univ v2;
	POOLALLOCATION(GMSTedgenode);
};

ALLOCATEPOOL(GMSTedgenode);

int GraphMst(const BGraph& undirectedg,
	     double (*fdist)(Univ v1, Univ v2),
	     BGraph& gnew)
{
	int nv=0, nebefore=0, neconsidered=0, neadded=0;
	Pqueue<GMSTedgenode*> pq;
	for (BGraphVertices gv(gnew);gv;gv.next()) {
		nv++;
		for (BGraphEdges ge(undirectedg,gv());ge;ge.next()) {
			Univ v1=gv(), v2=ge();
			if (v1<v2) continue;
			nebefore++;
			pq.enter(new GMSTedgenode(v1,v2),fdist(v1,v2));
		}
	}
	UnionFind<Univ> uf;
	for (;!pq.empty();) {
		GMSTedgenode* e=pq.removemin();
		neconsidered++;
		Univ v1=e->v1, v2=e->v2;
		delete e;
		if (!uf.unify(v1,v2)) continue;
		gnew.enteru(v1,v2);
		neadded++;
		if (neadded==nv-1) break;
	}
	for (;!pq.empty();) delete pq.removemin();
	SHOWF("GraphMst: %d vertices, %d/%d edges considered, %d output\n",
	      nv,neconsidered,nebefore,neadded);
	return neadded==nv-1;
}

BGraph* newGraphMst(const BGraph& undirectedg,
		    double (*fdist)(Univ v1, Univ v2))
{
	BGraph* gnew=new BGraph;
	for (BGraphVertices gv(undirectedg);gv;gv.next()) gnew->enter(gv());
	if (!GraphMst(undirectedg,fdist,*gnew)) delete gnew,gnew=0;
	return gnew;
}


//*** GraphMst PRIM

Graph<int>* newGraphMst(int num, double (*fdist)(int v1, int v2))
{
        int i;
        int j;
	const double INF=1e30f;
	SArray<double> lowcost(num);
	SArray<int> closest(num);
	Graph<int>* gnew=new Graph<int>;
//	for (int i=0;i<num;i++) gnew->enter(i);
	for (i=0;i<num;i++) gnew->enter(i);
	for (i=1;i<num;i++) {
		lowcost[i]=fdist(0,i);
		closest[i]=0;
	}
	for (i=1;i<num;i++) {
		double minf=INF; int minj=-1;
//		for (int j=1;j<num;j++)
		for (j=1;j<num;j++)
			if (lowcost[j]<minf) minf=lowcost[j],minj=j;
		assertx(minj>=0);
		gnew->enteru(minj,closest[minj]);
		lowcost[minj]=INF;
		for (j=1;j<num;j++) {
			if (lowcost[j]==INF) continue;
			double pnd=fdist(minj,j);
			if (pnd<lowcost[j]) lowcost[j]=pnd,closest[j]=minj;
		}
	}
	return gnew;
}


//*** GraphQuickMst

// Try to build the EMST of the num points pa using all edges with length less
// than thresh.
// Uses modified Prim's, where only edges of length<thresh are considered.
// Uses a HPqueue because it can no longer afford to find min in O(n) time.
static Graph<int>* newtryemst(double thresh, const Point pa[], int num,
			      const PointSpatial& sp)
{
        int i;
	Graph<int>* gnew=new Graph<int>;
	SArray<int> inset(num); // vertices already added to mst
	SArray<int> closest(num); // for !inset[i], closest inset[] so far
//	for (int i=0;i<num;i++) gnew->enter(i);
	for (i=0;i<num;i++) gnew->enter(i);
	for (i=0;i<num;i++) inset[i]=0;
	HPqueue<int> pq;
	pq.enter(0,0);
	for (;!pq.empty();) {
		int i=pq.removemin();
		assertx(i>=0 && i<num && !inset[i]); // optional
		if (i) gnew->enteru(i,closest[i]);
		inset[i]=1;
		SpatialSearch ss(sp,pa[i]);
		for (;;) {
			if (ss.done()) break;
			double dis2; int j=Conv<int>::d(ss.next(&dis2));
			if (dis2>square(thresh)) break;
			if (inset[j]) continue;
			double cd2=pq.retrieve(j);
			if (cd2<0) {
				pq.enter(j,dis2);
				closest[j]=i;
			} else if (dis2<cd2) {
				assertx(pq.update(j,dis2));
				closest[j]=i;
			}
		}
	}
	int allfound=1;
	for (i=0;i<num;i++) allfound&=inset[i];
	if (!allfound) delete gnew,gnew=0;
	return gnew;
}

Graph<int>* newGraphQuickEmst(const Point pa[], int num,
			      const PointSpatial& sp)
{
	Graph<int>* gnew;
	const double initf=.02;
	int n=0;
	for (double f=initf;;f*=1.6) {
		n++;
		if ((gnew=newtryemst(f,pa,num,sp)) != 0) break;
	}
	SHOWF("GraphQuickEmst: had to do %d approximate Emst's\n",n);
	return gnew;
}


Stat* newGraphEdgeStats(const BGraph& g, double (*fdist)(Univ v1, Univ v2))
{
	Stat* stat=new Stat;
	for (BGraphVertices gv(g);gv;gv.next())
		for (BGraphEdges ge(g,gv());ge;ge.next())
			*stat+=fdist(gv(),ge());
	return stat;
}


Graph<int>* newGraphEKClosest(const Point pa[], int num, int kcl,
			      const PointSpatial& sp)
{
        int i;
	Graph<int>* gnew=new Graph<int>;
//	for (int i=0;i<num;i++) gnew->enter(i);
	for (i=0;i<num;i++) gnew->enter(i);
	for (i=0;i<num;i++) {
		SpatialSearch ss(sp,pa[i]);
		for (int nn=0;nn<kcl+1;nn++) {
			assertx(!ss.done());
			int j=Conv<int>::d(ss.next());
			if (j==i) continue;
			gnew->enter(i,j);
		}
	}
	return gnew;
}


//*** GraphComponent

void GraphComponent::next()
{
	Queue<Univ> queue;
	set.enter(gv());
	queue.enqueue(gv());
	for (;!queue.empty();) {
		Univ v=queue.dequeue();
		for (BGraphEdges ge(g,v);ge;ge.next())
			if (set.add(ge())) queue.enqueue(ge());
	}
	for (gv.next();gv;gv.next())
		if (!set.contains(gv())) break;
}

int GraphNumComponents(const BGraph& g)
{
	int n=0;
	for (GraphComponent gc(g);gc;gc.next()) n++;
	return n;
}
