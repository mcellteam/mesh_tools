// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Contour.h"
#include "MeshOp.h"
#include "Stat.h"

using std::ostream;

class ContourNode {
    friend class Contour;
	ContourNode(int pen) : en(pen), evaled(0), cubestate(NOTHING) { }
	int en;			// encode vertex index
				//* vertex info
	int evaled;		// vertex value has been evaluated
	double val;		// vertex value
	Point p;		// position of point in grid
				//* cube info
	enum { NOTHING, QUEUED, VISITED } cubestate;
  private:
	POOLALLOCATION(ContourNode);
};

ALLOCATEPOOL(ContourNode);

const int MAXGN=894;		// cuberoot(2**31/3)

#define FORCUBEOFFSET(i,j,k) \
for (i=0;i<2;i++) \
for (j=0;j<2;j++) \
for (k=0;k<2;k++)

const double Contour::UNDEF=1e31f;

Contour::Contour(int pis3D, double (*pfeval)(const Point& p),
		 int pgn, ostream* pos)
: is3D(pis3D), feval(pfeval), gn(pgn), ios(pos), gnf(pgn), gni(1/gnf),
  foutputcontour(0), foutputborder(0), mesh(0), bigmeshfaces(0),
  vertextol(0), ncvisited(0), ncundef(0), ncnothing(0),
  nvevaled(0), nvzero(0), nvundef(0), nedegen(0)
{
	assertx(feval && gn>0);
	assertx(gn+1<=MAXGN);	// must leave room for (0..gn) inclusive
}

Contour::~Contour()
{
	ForMapValue(m,int,ContourNode*,n) { delete n; } EndFor;
	assertx(queue.empty());	// optional
	if (ios) {
		*ios<<hform("# March:\n");
		*ios<<hform("# visited %d cubes",ncvisited);
		*ios<<hform(" (%d were undefined, %d contained nothing)\n",
			    ncundef,ncnothing);
		*ios<<hform("# evaluated %d vertices",nvevaled);
		*ios<<hform(" (%d were zero, %d were undefined)\n",
			    nvzero,nvundef);
		*ios<<hform("# encountered %d tough edges\n",nedegen);
	}
}

void Contour::setOutputContour(void (*pfoutputcontour)(const Polygon& poly))
{
	assertx(!mesh);
	foutputcontour=pfoutputcontour;
}

void Contour::setOutputBorder(void (*pfoutputborder)(const Polygon& poly))
{
	foutputborder=pfoutputborder;
}

void Contour::setOutputMesh(GMesh& pmesh)
{
	assertx(is3D && !foutputcontour);
	mesh=&pmesh;
}

void Contour::bigMeshFaces()
{
	bigmeshfaces=1;
}

void Contour::setVertexTolerance(double tol)
{
	vertextol=tol;
}

inline int Contour::encode(int ci[3]) const
{
	return (ci[0]*MAXGN+ci[1])*MAXGN+ci[2];
}

inline void Contour::decode(int en, int ci[3]) const
{
	ci[2]=en%MAXGN; en/=MAXGN;
	ci[1]=en%MAXGN; en/=MAXGN;
	ci[0]=en;
}

inline int Contour::inbound(int ci[3]) const
{
	return (ci[0]>=0 && ci[1]>=0 && ci[2]>=0 &&
		ci[0]<gn && ci[1]<gn && ci[2]<gn);
}

int Contour::marchFrom(const Point& startp)
{
	assertx(startp[0]>=0 && startp[1]>=0 && startp[2]>=0);
	assertx(startp[0]<=1 && startp[1]<=1 && startp[2]<=1);
	if (!is3D) assertx(!startp[0]);
	int cc[3];
	for (int c=0;c<3;c++) cc[c]=::min(int(startp[c]*gnf),gn-1);
	return marchfrom(cc);
}

int Contour::marchNear(const Point& startp)
{
	assertx(startp[0]>=0 && startp[1]>=0 && startp[2]>=0);
	assertx(startp[0]<=1 && startp[1]<=1 && startp[2]<=1);
	if (!is3D) assertx(!startp[0]);
	int cc[3];
	for (int c=0;c<3;c++) cc[c]=::min(int(startp[c]*gnf),gn-1);
	int ret=0;
	for (int i=-1;i<2;i++)
		for (int j=-1;j<2;j++)
			for (int k=-1;k<2;k++) {
				if (!is3D && i) continue;
				int ci[3];
				ci[0]=cc[0]+i; ci[1]=cc[1]+j; ci[2]=cc[2]+k;
				if (inbound(ci)) ret+=marchfrom(ci);
			}
	return ret;
}

int Contour::marchfrom(int cc[3])
{
	int cncvisited=ncvisited;
	if (!is3D) assertx(!cc[0]);
	int en=encode(cc);
	ContourNode* n=m.retrieve(en);
	if (!n) m.enter(en,n=new ContourNode(en));
	if (n->cubestate==ContourNode::VISITED) return 0;
	assertx(n->cubestate==ContourNode::NOTHING);
	queue.enqueue(en);
	n->cubestate=ContourNode::QUEUED;
	for (;!queue.empty();) {
		en=queue.dequeue();
		considercube(en);
	}
	if (ncvisited-cncvisited==1) ncnothing++;
	return ncvisited-cncvisited;
}

void Contour::considercube(int encube)
{
        int c;
	ncvisited++;
	int cc[3]; decode(encube,cc);
	if (!is3D) assertx(!cc[0]);
	ContourNode* na[2][2][2];
	int cundef=0;
	int cd[3];
	FORCUBEOFFSET(cd[0],cd[1],cd[2]) {
		if (!is3D && cd[0]) continue;
		int ci[3];
//		for (int c=0;c<3;c++) ci[c]=cc[c]+cd[c];
		for (c=0;c<3;c++) ci[c]=cc[c]+cd[c];
		int en=encode(ci);
		ContourNode* n=m.retrieve(en);
		if (!n) m.enter(en,n=new ContourNode(en));
		na[cd[0]][cd[1]][cd[2]]=n;
		if (!n->evaled) {
			n->evaled=1;
			for (c=0;c<3;c++) n->p[c]=::min(ci[c]*gni,1.);
			n->val=feval(n->p);
			nvevaled++;
			if (!n->val) nvzero++;
			if (n->val==UNDEF) nvundef++;
		}
		if (n->val==UNDEF) cundef=1;
	}
	ContourNode* n=na[0][0][0];
	assertx(n->cubestate==ContourNode::QUEUED);
	n->cubestate=ContourNode::VISITED;
	if (cundef) {
		ncundef++;
	} else if (foutputcontour || mesh) {
		is3D?contourcube(na):contoursquare(na);
	}
	for (int d=is3D?0:1;d<3;d++)
		pushneighbors(d,cc,na);
}

void Contour::pushneighbors(int d, const int cc[3], ContourNode* na[2][2][2])
{
        int c;
	int d1=(d+1)%3;
	int d2=(d+2)%3;
	int cd[3];
	for (int i=0;i<2;i++) {
		double fmin=1e30,fmax=-1e30;
		cd[d]=i;
		int fundef=0;
		for (cd[d1]=0;cd[d1]<2;cd[d1]++) {
			for (cd[d2]=0;cd[d2]<2;cd[d2]++) {
				if (!is3D && cd[0]) continue;
				double f=na[cd[0]][cd[1]][cd[2]]->val;
				if (f==UNDEF) fundef=1;
				if (f<fmin) fmin=f;
				if (f>fmax) fmax=f;
			}
		}
		cd[d]=i?1:-1;
		cd[d1]=cd[d2]=0;
		int ci[3]; // indice of node for neighboring cube;
//		for (int c=0;c<3;c++) ci[c]=cc[c]+cd[c];
		for (c=0;c<3;c++) ci[c]=cc[c]+cd[c];
		// note: fmin<0 since 0 is arbitrarily taken to be positive
		if (!fundef && fmin<0 && fmax>=0 && inbound(ci)) {
			int en=encode(ci);
			ContourNode* n=m.retrieve(en);
			if (!n) m.enter(en,n=new ContourNode(en));
			if (n->cubestate==ContourNode::NOTHING) {
				n->cubestate=ContourNode::QUEUED;
				queue.enqueue(en);
			}
		} else if (foutputborder) {
			// output boundary
			cd[d]=i;
			static Polygon poly; // not reentrant
			poly.init(0);
			for (cd[d1]=0;cd[d1]<2;cd[d1]++) {
				int sw=cd[d]^cd[d1]; // 0 or 1
				for (cd[d2]=sw;cd[d2]==0||cd[d2]==1;
				     cd[d2]+=sw?-1:1) {
					if (!is3D && cd[0]) continue;
					Point p;
					for (c=0;c<3;c++) ci[c]=cc[c]+cd[c];
					for (c=0;c<3;c++) p[c]=ci[c]*gni;
					poly+=p;
				}
			}
			foutputborder(poly);
		}
	}
}

void Contour::contourcube(ContourNode* na[2][2][2])
{
	if (mesh) { cubemesh(na); return; }
	// do Kuhn 6-to-1 triangulation of cube
	ContourNode* n4[4]; // tetrahedron
	n4[0]=na[0][0][0]; n4[1]=na[0][0][1];
	n4[2]=na[1][0][1]; n4[3]=na[0][1][0];
	contourtetrahedron(n4);
	n4[0]=na[0][0][0]; n4[1]=na[1][0][1];
	n4[2]=na[1][0][0]; n4[3]=na[0][1][0];
	contourtetrahedron(n4);
	n4[0]=na[1][0][1]; n4[1]=na[1][1][0];
	n4[2]=na[1][0][0]; n4[3]=na[0][1][0];
	contourtetrahedron(n4);
	n4[0]=na[0][1][0]; n4[1]=na[0][1][1];
	n4[2]=na[0][0][1]; n4[3]=na[1][0][1];
	contourtetrahedron(n4);
	n4[0]=na[1][1][1]; n4[1]=na[0][1][1];
	n4[2]=na[0][1][0]; n4[3]=na[1][0][1];
	contourtetrahedron(n4);
	n4[0]=na[1][1][1]; n4[1]=na[0][1][0];
	n4[2]=na[1][1][0]; n4[3]=na[1][0][1];
	contourtetrahedron(n4);
}

void Contour::contourtetrahedron(ContourNode* n4[4])
{
        int i;
	int nposi=0;
//	for (int i=0;i<4;i++) if (n4[i]->val>=0) nposi++;
	for (i=0;i<4;i++) if (n4[i]->val>=0) nposi++;
	if (nposi==0 || nposi==4) return;
	int j;
	for (i=0,j=3;i<j;) {
		if (n4[i]->val>=0) { i++; continue; }
		if (n4[j]->val<0) { j--; continue; }
		swap(&n4[i],&n4[j]);
		i++; j--;
	}
	ContourNode* n3[3][2];
	switch (nposi) {
	  case 1:
		n3[0][0]=n4[0]; n3[0][1]=n4[1];
		n3[1][0]=n4[0]; n3[1][1]=n4[2];
		n3[2][0]=n4[0]; n3[2][1]=n4[3];
		outputtriangle(n3);
	  bcase 2:
		n3[0][0]=n4[0]; n3[0][1]=n4[2];
		n3[1][0]=n4[0]; n3[1][1]=n4[3];
		n3[2][0]=n4[1]; n3[2][1]=n4[3];
		outputtriangle(n3);
		n3[0][0]=n4[0]; n3[0][1]=n4[2];
		n3[1][0]=n4[1]; n3[1][1]=n4[3];
		n3[2][0]=n4[1]; n3[2][1]=n4[2];
		outputtriangle(n3);
	  bcase 3:
		n3[0][0]=n4[0]; n3[0][1]=n4[3];
		n3[1][0]=n4[1]; n3[1][1]=n4[3];
		n3[2][0]=n4[2]; n3[2][1]=n4[3];
		outputtriangle(n3);
	  bdefault: assertnever("internal");
	}
}

void Contour::outputtriangle(ContourNode* n3[3][2])
{
	static Polygon poly(3);
	for (int i=0;i<3;i++) {
		ContourNode* np=n3[i][0]; ContourNode* nn=n3[i][1];
		poly[i]=computepoint(np->p,nn->p,np->val,nn->val);
	}
	Vector normal=cross(poly[0],poly[1],poly[2]);
	if (dot(normal,n3[0][0]->p-n3[0][1]->p)<0) swap(&poly[0],&poly[1]);
	foutputcontour(poly);
}

// Find the intersection point on the segment between pp and pn whose values
//  are vp and vn.
// Return the point as p.
// If maxdiff is nonzero, return a positive pointing normal if the gradient
//  along this edge is greater than previously seen (as determined by maxdiff).
// To handle degeneracies for non-mesh output, require intersections to be at
//  least a fixed percentage distance away from the vertices (fdegen).
// If vertextol, do binary search until uncertainty is less than vertextol.
Point Contour::computepoint(const Point& pp, const Point& pn,
			    double vp, double vn)
{
        int neval;
	Point pm; double fm;
	NEST {
		double v0=vp, v1=vn;
		Point p0=pp, p1=pn;
		double f0=0, f1=1;
//		for (int neval=0;;) {
		for (neval=0;;) {
			assertx(v0>=0 && v1<0 && f0<f1);
			double b1=v0/(v0-v1);
			// guarantee quick convergence
			if (vertextol) b1=::min(::max(b1,.05),.95);
			fm=f0*(1-b1)+f1*b1;
			pm=interp(p0,p1,1-b1);
			if (!vertextol) break;
			double vm=feval(pm);
			neval++;
			if (vm<0) f1=fm,p1=pm,v1=vm;
			else f0=fm,p0=pm,v0=vm;
			if (dist2(p0,p1)<=square(vertextol)) break;
		}
		SSTAT(SContneval,neval);
	}
	if (!mesh) {
		const double fs=.001;
		if (fm<fs) fm=fs,nedegen++,pm=interp(pp,pn,1-fm);
		else if (fm>1-fs) fm=1-fs,nedegen++,pm=interp(pp,pn,1-fm);
	}
	return pm;
}


void Contour::contoursquare(ContourNode* na[2][2][2])
{
	ContourNode* n3[3];
	n3[0]=na[0][0][0]; n3[1]=na[0][1][1]; n3[2]=na[0][0][1];
	contourtriangle(n3);
	n3[0]=na[0][0][0]; n3[1]=na[0][1][0]; n3[2]=na[0][1][1];
	contourtriangle(n3);
}

void Contour::contourtriangle(ContourNode* n3[3])
{
        int i;
	int nposi=0;
//	for (int i=0;i<3;i++) if (n3[i]->val>=0) nposi++;
	for (i=0;i<3;i++) if (n3[i]->val>=0) nposi++;
	if (nposi==0 || nposi==3) return;
	int j;
	for (i=0,j=2;i<j;) {
		if (n3[i]->val>=0) { i++; continue; }
		if (n3[j]->val<0) { j--; continue; }
		swap(&n3[i],&n3[j]);
		i++; j--;
	}
	ContourNode* n2[2][2];
	switch (nposi) {
	  case 1:
		n2[0][0]=n3[0]; n2[0][1]=n3[1]; n2[1][0]=n3[0]; n2[1][1]=n3[2];
		outputline(n2);
	  bcase 2:
		n2[0][0]=n3[0]; n2[0][1]=n3[2]; n2[1][0]=n3[1]; n2[1][1]=n3[2];
		outputline(n2);
	  bdefault: assertnever("internal");
	}
}

void Contour::outputline(ContourNode* n2[2][2])
{
	static Polygon poly(2);
	for (int i=0;i<2;i++) {
		ContourNode* np=n2[i][0]; ContourNode* nn=n2[i][1];
		poly[i]=computepoint(np->p,nn->p,np->val,nn->val);
	}
	Vector v=poly[1]-poly[0];
	assertx(!v[0]);		// optional
	Vector normal(0,-v[2],v[1]); // 90 degree rotation
	if (dot(normal,n2[0][0]->p-n2[0][1]->p)<0) swap(&poly[0],&poly[1]);
	foutputcontour(poly);
}

void Contour::cubemesh(ContourNode* na[2][2][2])
{
	Map<Vertex,Vertex> mapsucc;
	ForIndex(d,3) {
		ForIndex(v,2) {
			examineface(d,v,na,mapsucc);
		} EndFor;
	} EndFor;
	for (;!mapsucc.empty();) {
		Vertex vf=0; int minvi=1<<30; // find min to be portable
		ForMapKey(mapsucc,Vertex,Vertex,v) {
			int vi=mesh->vertexid(v);
			if (vi<minvi) minvi=vi,vf=v;
		} EndFor;
		static Array<Vertex> va; va.init(0);
		for (Vertex v=vf;;) {
			va+=v;
			v=assertv(mapsucc.remove(v));
			if (v==vf) break;
		}
		Face f=assertv(mesh->createFace(va,va.num()));
		if (va.num()>3 && !bigmeshfaces) {
			// If 6 or more edges, may have 2 edges on same cube
			//  face, then must introduce new vertex to be safe.
			if (va.num()>=6) (void)CenterSplitFace(*mesh,f);
			else assertx(TriangulateFace(*mesh,f));
		}
	}
}

void Contour::examineface(int d, int v, ContourNode* na[2][2][2],
			  Map<Vertex,Vertex>& mapsucc)
{
	int d1=(d+1)%3, d2=(d+2)%3, cd[3];
	cd[d]=v;
	ContourNode* naf[4];
	int i=0;
	// Gather 4 cube vertices in a consistent order
	for (cd[d1]=0;cd[d1]<2;cd[d1]++) {
		int sw=cd[d]^cd[d1]; // 0 or 1
		for (cd[d2]=sw;cd[d2]==0||cd[d2]==1;
		     cd[d2]+=sw?-1:1) {
			naf[i++]=na[cd[0]][cd[1]][cd[2]];
		}
	}
	int nneg=0;
	double sumval=0;
	for (i=0;i<4;i++) {
		double val=naf[i]->val;
		if (val<0) nneg++;
		sumval+=val;
	}
	// If pedantic, could sort the vals before summing.
	for (i=0;i<4;i++) {
		int i1=(i+1)%4, i2=(i+2)%4, i3=(i+3)%4;
		if (!(naf[i]->val<0 && naf[i1]->val>=0)) continue;
		// have start of edge
		assertx(nneg>=1 && nneg<=3); // optional
		int ie;	// end of edge
		if (nneg==1) {
			ie=i3;
		} else if (nneg==3) {
			ie=i1;
		} else if (naf[i2]->val>=0) {
			ie=i2;
		} else if (sumval<0) {
			ie=i1;
		} else {
			ie=i3;
		}
		Vertex v1=getvertexonedge(naf[i1],naf[i]);
		Vertex v2=getvertexonedge(naf[ie],naf[(ie+1)%4]);
		mapsucc.enter(v2,v1);		 // to get face order correct
	}
}

// First node is positive, second node is negative
Vertex Contour::getvertexonedge(ContourNode* n1, ContourNode* n2)
{
	int ene;
	NEST {
		int cc1[3],cc2[3];
		decode(n1->en,cc1);
		decode(n2->en,cc2);
		int d=-1;
		for (int c=0;c<3;c++) if (cc1[c]!=cc2[c]) assertx(d<0),d=c;
		assertx(d>=0);
		assertx(abs(cc1[d]-cc2[d])==1);
		ContourNode* n=(cc1[d]<cc2[d])?n1:n2;
		ene=n->en*3+d;
	}
	Vertex v=mapenev.retrieve(ene);
	if (!v) {
		mapenev.enter(ene,v=mesh->createVertex());
		mesh->setPoint(v,computepoint(n1->p,n2->p,n1->val,n2->val));
	}
	return v;
}
