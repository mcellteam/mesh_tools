// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Mesh.h"
#include "Set.h"
#include "Array.h"
#include "Random.h"

ALLOCATEPOOL(MVertex);
ALLOCATEPOOL(MFace);
ALLOCATEPOOL(MEdge);

// Note: all changes to mesh go through enterHEdge() and removeHEdge() !!
//  How elegant that is!

//            v2
//       /    ^    \        *
//   vo1  f1 e|  f2  vo2
//       \    |    /
//            v1            if (isBoundary(e)), f2==0, vo2==0

//  Vertex v1 -> stack of HEdges (v1,?)

//*** Mesh

int Mesh::sdebug=GetenvValue("DEBUG");

Mesh::Mesh() : vertexnum(1), facenum(1), nedges(0) { }

Mesh::~Mesh()
{
	clear();
}

void Mesh::clear()
{
	for (;; ) {
		MapIter<int,Face> mi(id2face);
		if (!mi) break;
		destroyFace(mi.value());
	}
	for (;; ) {
		MapIter<int,Vertex> mi(id2vertex);
		if (!mi) break;
		destroyVertex(mi.value());
	}
	vertexnum=1;
	facenum=1;
	nedges=0;
}

void Mesh::copy(const Mesh& m)
{
	int oldv=setDebug(-1);
	copyVertices(m);
	Array<Vertex> va;
	ForMeshFace(m,f) {
		m.vertices(f,va);
		for (int i=0;i<va.num();i++)
			va[i]=idvertex(m.vertexid(va[i]));
		assertx(createFaceI(m.faceid(f),va,va.num()));
	} EndFor;
	setDebug(oldv);
}

void Mesh::copyVertices(const Mesh& m)
{
	assertx(!numVertices());
	ForMeshVertex(m,v) {
		(void)createVertexI(m.vertexid(v));
	} EndFor;
}

//** Raw manipulation

Vertex Mesh::createVertex()
{
	return createVertexI(vertexnum);
}

Vertex Mesh::createVertexI(int id)
{
	assertx(id>=1);
	Vertex v=new MVertex;
	// v->ste set to empty
	v->id=id;
	v->flags=0;
	v->point=Point(0,0,0);
	if (sdebug>=0) assertx(!id2vertex.contains(id));
	id2vertex.enter(id,v);
	vertexnum=max(vertexnum,id+1);
	return v;
}

void Mesh::destroyVertex(Vertex v)
{
	assertx(v->ste.empty());
	assertx(id2vertex.remove(v->id));
	if (vertexnum-1==v->id) vertexnum--;
	delete v;
}

Face Mesh::createFace(Vertex va[], int nv)
{
	return createFaceI(facenum,va,nv);
}

Face Mesh::createFace(Vertex v1, Vertex v2, Vertex v3)
{
	Vertex va[3];
	va[0]=v1; va[1]=v2; va[2]=v3;
	return createFace(va,3);
}

Face Mesh::createFaceI(int id, Vertex va[], int nv)
{
	assertx(id>=1);
	assertx(nv>=3);
	if (sdebug<0) {
		// trust everything is ok
	} else if (nv==3) {	// cheap check
		if (va[0]==va[1] || va[1]==va[2] || va[0]==va[2]) return 0;
	} else {
		Set<Vertex> setv;
		for (int i=0;i<nv;i++) {
			Vertex v=va[i];
			valid(v);
			if (!setv.add(v)) return 0;
		}
	}
	if (sdebug>=0) {
		Vertex vo=va[nv-1];
		for (int i=0;i<nv;i++) {
			if (lookupEdge(vo,va[i])) return 0;
			vo=va[i];
		}
	}
	// all systems go.
	Face f=new MFace;
	// f->erep defined below
	f->id=id;
	f->flags=0;
	if (sdebug>=1) assertx(!id2face.contains(id));
	id2face.enter(id,f);
	MEdge edummy;
	Edge ep=&edummy;
	edummy.vert=va[0];
	for (int i=0;i<nv;i++) {
		Vertex v2=va[i+1==nv?0:i+1];
		Edge e=new MEdge;
		e->prev=ep;
		// e->next is set below
		// e->sym is set in enterHEdge
		e->vert=v2;
		e->face=f;
		e->flags=0;
		// assertx(e->prev->vert==v1);
		enterHEdge(e);
		ep=e;
	}
	Edge elast=ep;
	for (;;) {
		if (ep->prev==&edummy) {
			ep->prev=elast;
			elast->next=ep;
			break;
		}
		ep->prev->next=ep;
		ep=ep->prev;
	}
	f->erep=elast;		// any edge is ok here
	facenum=max(facenum,id+1);
	if (sdebug>=3) OK();
	return f;
}

void Mesh::destroyFace(Face f)
{
	Stack<Edge> stack;
	ForFaceHEdge(*this,f,e) {
		removeHEdge(e);
		stack.push(e);
	} EndFor;
	while (!stack.empty()) delete stack.pop();
	assertx(id2face.remove(f->id));
	if (facenum-1==f->id) facenum--;
	delete f;
	if (sdebug>=3) OK();
}

int Mesh::substituteFaceVertex(Face f, Vertex vold, Vertex vnew)
{
	assertx(vold!=vnew);
	if (sdebug>=1) valid(vold),valid(vnew);
	ForFaceVertex(*this,f,v) {
		if (v==vnew) return 0;
	} EndFor;
	Edge e=unrepv2(clwEdge(vold,f),vold);
	assertx(e->vert==vold);	// optional
	if (lookupEdge(e->prev->vert,vnew)) return 0;
	if (lookupEdge(vnew,e->next->vert)) return 0;
	// all systems go
	removeHEdge(e);
	removeHEdge(e->next);
	e->vert=vnew;
	enterHEdge(e);
	enterHEdge(e->next);
	return 1;
}

//** Vertex

int Mesh::isNice(Vertex v) const
{
        Edge e;
	Edge er=erep(v);
	if (!er) return 1;
	int ne=0;
//	for (Edge e=er;;) {
	for (e=er;;) {
		ne++;
		e=clwHEdge(e);
		if (!e || e==er) break;
	}
//	if (e!=er) for (e=er;(e=ccwHEdge(e)) != 0;) ne++;
	if (ne!=v->ste.height())
            printf("Mesh::isNice - ne!=v->ste.height(): vertex [%.15lg %.15lg %.15lg]\n", v->point[0], v->point[1], v->point[2]);
	return ne==v->ste.height();
}

int Mesh::degree(Vertex v) const
{
	int i=0;
	ForStack(v->ste,Edge,e) i+=isBoundary(e)?2:1; EndFor;
	return i;
}

int Mesh::numBoundaries(Vertex v) const
{
	int i=0;
	ForStack(v->ste,Edge,e) {
		if (isBoundary(e)) i++;
	} EndFor;
	return i;
}

int Mesh::isBoundary(Vertex v) const
{
	assertx(erep(v) && isNice(v));
	return numBoundaries(v);
}

Edge Mesh::oppEdge(Vertex v, Face f) const
{
        int i;
	assertx(isTriangle(f));
	Edge e=erep(f);
//	for (int i=0;i<3;i++,e=e->next)
	for (i=0;i<3;i++,e=e->next)
		if (e->vert==v) break;
	assertx(i<3);		// optional
	return rep(e->prev);
}

Vertex Mesh::oppVertex(Vertex v, Edge e) const
{
	if (vertex1(e)==v) return vertex2(e);
	if (vertex2(e)==v) return vertex1(e);
	assertnever("Vertex not on Edge"); return 0;
}

int Mesh::vertices(Vertex v, Array<Vertex>& va) const
{
        Vertex ve;
	va.init(0);
	Vertex v0=mostClwVertex(v);
//	for (Vertex ve=v0;ve;) {
	for (ve=v0;ve;) {
		va+=ve;
		if ((ve=ccwVertex(v,ve))==v0) break;
	}
	return ve?1:0;
}

int Mesh::faces(Vertex v, Array<Face>& fa) const
{
	fa.init(0);
	Edge e=mostClwHEdge(v);
	for (;e;) {
		Face f=e->face;
		if (fa.num() && f==fa[0]) break;
		fa+=f;
		e=ccwHEdge(e);
	}
	return e?1:0;
}

Vertex Mesh::mostClwVertex(Vertex v) const
{
	Edge e=mostClwHEdge(v);
	return e?e->next->vert:0;
}

Vertex Mesh::mostCcwVertex(Vertex v) const
{
	Edge e=mostCcwHEdge(v);
	return e?e->prev->vert:0;
}

Vertex Mesh::clwVertex(Vertex v, Vertex vext) const
{
	Edge e=lookupEdge(vext,v);
	return e?e->next->vert:0;
}

Vertex Mesh::ccwVertex(Vertex v, Vertex vext) const
{
	Edge e=lookupEdge(v,vext);
	return e?e->prev->prev->vert:0;
}

Face Mesh::mostClwFace(Vertex v) const
{
	Edge e=mostClwHEdge(v);
	return e?e->face:0;
}

Face Mesh::mostCcwFace(Vertex v) const
{
	Edge e=mostCcwHEdge(v);
	return e?e->face:0;
}

Face Mesh::clwFace(Vertex v, Face f) const
{
	return oppFace(f,ccwEdge(v,f));
}

Face Mesh::ccwFace(Vertex v, Face f) const
{
	return oppFace(f,clwEdge(v,f));
}

Edge Mesh::mostClwEdge(Vertex v) const
{
	Edge e=mostClwHEdge(v);
	return e?rep(e->next):0;
}

Edge Mesh::mostCcwEdge(Vertex v) const
{
	Edge e=mostCcwHEdge(v);
	return e?rep(e):0;
}

Edge Mesh::clwEdge(Vertex v, Edge e) const
{
	Edge he=unrepv2(e,v);
	return he?rep(he->next):0;
}

Edge Mesh::ccwEdge(Vertex v, Edge e) const
{
	Edge he=unrepv1(e,v);
	return he?rep(he->prev):0;
}

//** Face

int Mesh::isNice(Face f) const
{
	// non-nice iff all its neighbors are the same
	Face fp=0; int i=0;
	ForFaceFace(*this,f,f2) {
		if (++i==1) { fp=f2; }
		else if (i==2) return f2!=fp;
	} EndFor;
	return 1;
}

int Mesh::numVertices(Face f) const
{
	int i=0;
	ForFaceHEdge(*this,f,e) {
		(void)e;
		i++;
	} EndFor;
	return i;
}

int Mesh::numBoundaryEdges(Face f) const
{
	int i=0;
	ForFaceHEdge(*this,f,e) {
		if (isBoundary(e)) i++;
	} EndFor;
	return i;
}

int Mesh::isTriangle(Face f) const
{
	Edge e=erep(f);
	return e->next->next->next==e;
}

int Mesh::isBoundary(Face f) const
{
	ForFaceHEdge(*this,f,e) {
		if (isBoundary(e)) return 1;
	} EndFor;
	ForFaceVertex(*this,f,v) {
		if (isBoundary(v)) return 1;
	} EndFor;
	return 0;
}

Face Mesh::oppFace(Face f, Edge e) const
{
	if (face1(e)==f) return face2(e);
	if (face2(e)==f) return face1(e);
	assertnever("Face not adjacent to Edge"); return 0;
}

void Mesh::vertices(Face f, Array<Vertex>& va) const
{
	va.init(0);
	ForFaceVertex(*this,f,v) {
		va+=v;
	} EndFor;
}

void Mesh::vertices(Face f, Vertex va[3]) const
{
	Edge e=erep(f), ef=e;
	va[0]=e->vert; e=e->next;
	va[1]=e->vert; e=e->next;
	va[2]=e->vert; e=e->next;
	assertx(e==ef);
}

Edge Mesh::clwEdge(Edge e, Face f) const
{
	if (face1(e)==f) return rep(e->prev);
	if (face2(e)==f) return rep(e->sym->prev);
	assertnever("Face not adjacent to Edge"); return 0;
}

Edge Mesh::ccwEdge(Edge e, Face f) const
{
	if (face1(e)==f) return rep(e->next);
	if (face2(e)==f) return rep(e->sym->next);
	assertnever("Face not adjacent to Edge"); return 0;
}

//** Edge

int Mesh::isBoundary(Edge e) const
{
	return e->sym?0:1;
}

Vertex Mesh::vertex1(Edge e) const
{
	return e->prev->vert;
}

Vertex Mesh::vertex2(Edge e) const
{
	return e->vert;
}

Face Mesh::face1(Edge e) const
{
	return e->face;
}

Face Mesh::face2(Edge e) const
{
	return e->sym?e->sym->face:0;
}

Vertex Mesh::sideVertex1(Edge e) const
{
	assertx(e->next->next->next==e);
	return e->next->vert;
}

Vertex Mesh::sideVertex2(Edge e) const
{
	if (!e->sym) return 0;
	return sideVertex1(e->sym);
}

Vertex Mesh::sideVertex(Edge e, Face f) const
{
	if (face1(e)==f) return sideVertex1(e);
	if (face2(e)==f) return sideVertex2(e);
	assertnever("Face not adjacent to Edge"); return 0;
}

Edge Mesh::oppBoundary(Edge e, Vertex v) const
{
	assertx(isBoundary(e));
	Edge enext;
	if (vertex1(e)==v) {
		e=e->prev;
		for (;(enext=ccwHEdge(e)) != 0;) e=enext;
	} else if (vertex2(e)==v) {
		for (;(enext=clwHEdge(e)) != 0;) e=enext;
		e=e->next;
	} else { assertnever("Vertex not on Edge"); }
	return rep(e);
}

Edge Mesh::clwBoundary(Edge e) const
{
	return oppBoundary(e,vertex2(e));
}

Edge Mesh::ccwBoundary(Edge e) const
{
	return oppBoundary(e,vertex1(e));
}

//** Other associations

Edge Mesh::queryEdge(Vertex v, Vertex w) const
{
	Edge e=lookupEdge(v,w);
	if (e) return rep(e);
	e=lookupEdge(w,v);
	if (e) assertx(e==rep(e)); // optional
	return e;
}

Edge Mesh::edge(Vertex v, Vertex w) const
{
	return assertv(queryEdge(v,w));
}

Edge Mesh::orderedEdge(Vertex v1, Vertex v2) const
{
	Edge e=assertv(lookupEdge(v1,v2));
	assertx(e==rep(e));
	return e;
}

Face Mesh::face(Vertex v1, Vertex v2, Vertex v3) const
{
	Edge e=edge(v1,v2);
	Face f;
	f=face1(e); if (f && isTriangle(f) && sideVertex1(e)==v3) return f;
	f=face2(e); if (f && isTriangle(f) && sideVertex2(e)==v3) return f;
	assertnever("Face not adjacent to Vertices"); return 0;
}

Edge Mesh::clwEdge(Vertex v, Face f) const
{
	Edge ef=0;
	ForVertexHEdge(*this,v,e) {
		if (e->face==f) { ef=e; break; }
	} EndFor;
	if (!ef) assertnever("Face not adjacent to Vertex");
	return rep(ef);
}

Edge Mesh::ccwEdge(Vertex v, Face f) const
{
	Edge ef=0;
	ForVertexHEdge(*this,v,e) {
		if (e->face==f) { ef=e; break; }
	} EndFor;
	if (!ef) assertnever("Face not adjacent to Vertex");
	return rep(ef->next);
}

//** Counting routines

int Mesh::numVertices() const
{
	return id2vertex.num();
}

int Mesh::numFaces() const
{
	return id2face.num();
}

int Mesh::numEdges() const
{
	return nedges;
}

Vertex Mesh::randomVertex(Random& r) const
{
	MapIter<int,Vertex> mi(id2vertex,r);
	assertx(mi); return mi.value();
}

Face Mesh::randomFace(Random& r) const
{
	MapIter<int,Face> mi(id2face,r);
	assertx(mi); return mi.value();
}

Edge Mesh::randomEdge(Random& r) const
{
	Face f=randomFace(r);
	int vi=r.getint()%numVertices(f);
	Edge ef=0;
	ForFaceEdge(*this,f,e) {
		if (!vi--) { ef=e; break; }
	} EndFor;
	return assertv(ef);
}

//** Mesh operations

int Mesh::niceEdgeCollapse(Edge e) const
{
	if (sdebug>=1) valid(e);
	Vertex v1=vertex1(e), v2=vertex2(e);
	Face f1=face1(e), f2=face2(e);
	assertx(isTriangle(f1));
	if (f2) assertx(isTriangle(f2));
	if (sdebug>=1) assertx(isNice(v1) && isNice(v2));
	// Requirements:
	//* 1 - If v1 and v2 are both boundary, (v1,v2) is a boundary edge
	if (!isBoundary(e) && isBoundary(v1) && isBoundary(v2)) return 0;
	//* 2 - For all vertices adjacent to both v1 and v2, exists a face
	Vertex vo1=sideVertex1(e), vo2=sideVertex2(e);
	Set<Vertex> set;
	ForVertexVertex(*this,v1,v) {
		if (v!=vo1 && v!=vo2) set.enter(v);
	} EndFor;
	ForVertexVertex(*this,v2,v) {
		if (v!=vo1 && v!=vo2 && !set.add(v)) return 0;
	} EndFor;
	//* 3 - two small base cases: single face and tetrahedron
	if (set.num()==2 && isBoundary(e))
		return 0;	// single face
	if (set.num()==2 && !isBoundary(v1) && !isBoundary(v2))
		return 0;	// tetrahedron
	return 1;
}

int Mesh::legalEdgeCollapse(Edge e) const
{
	if (sdebug>=1) valid(e);
	Vertex v1=vertex1(e), v2=vertex2(e);
	Vertex vo1=sideVertex1(e), vo2=sideVertex2(e); // vo2 may be 0
	// Check that substituting v2 to v1 will not duplicate an edge in any
	// faces adjacent to v2 (besides f1 and f2).
	// (case of Face vertices being duplicated cannot happen here
	//  since only f1 and f2 can have both v1 and v2)
	ForVertexVertex(*this,v2,v) {
		if (v==v1 || v==vo1 || v==vo2) continue;
		if (queryEdge(v,v1)) return 0;
	} EndFor;
	return 1;
}

int Mesh::legalEdgeSwap(Edge e) const
{
	if (isBoundary(e)) return 0;
	// illegal if cross edge already exists (as in tetrahedron)
	if (queryEdge(sideVertex1(e),sideVertex2(e))) return 0;
	return 1;
}

void Mesh::collapseEdge(Edge e)
{
	if (sdebug>=1) valid(e);
	Vertex v1=vertex1(e), v2=vertex2(e);
	Face f1=face1(e), f2=face2(e);
	assertx(isTriangle(f1));
	if (f2) assertx(isTriangle(f2));
	// Destroy f1+f2, non-nice vertices may result
	destroyFace(f1);
	if (f2) destroyFace(f2);
	// Change remaining faces around v2 to have v1 instead
	ForVertexFace(*this,v2,f) {
		assertx(substituteFaceVertex(f,v2,v1));	// may die if illegal
	} EndFor;
	// Destroy vertex v2
	destroyVertex(v2);
}

Vertex Mesh::splitEdge(Edge e, int id)
{
	if (sdebug>=1) valid(e);
	Vertex v2=vertex2(e);
	Vertex vo1=sideVertex1(e), vo2=sideVertex2(e); // implies triangles
	// Create new vertex
	Vertex v=id?createVertexI(id):createVertex();
	// In faces adjacent to e, change references to v2 to v.
	Stack<Face> stackf;
	ForEdgeFace(*this,e,f) {
		stackf.push(f);
	} EndFor;
	while (!stackf.empty())
		assertx(substituteFaceVertex(stackf.pop(),v2,v));
	// Fill in the 3 or 4 sided hole (vo1,v,vo2,v2) with two new faces
	assertx(createFace(v,v2,vo1));
	if (vo2) assertx(createFace(v,vo2,v2));
	return v;
}

Edge Mesh::swapEdge(Edge e)
{
	if (sdebug>=1) valid(e);
	assertx(!isBoundary(e));
	Vertex v1=vertex1(e), v2=vertex2(e);
	Face f1=face1(e), f2=face2(e);
	Vertex vo1=sideVertex1(e), vo2=sideVertex2(e);
	// Destroying f1+f2 may create non-nice vertices, but no problem
	destroyFace(f1);
	destroyFace(f2);
	// createFace()'s may die if illegal
	assertx(createFace(v1,vo2,vo1));
	assertx(createFace(v2,vo1,vo2));
	return edge(vo1,vo2);
}

void Mesh::OK() const
{
	// Check that sets and maps have unique entries
	id2vertex.OK();
	id2face.OK();
	// Check consistency of id2x (one way)
	ForMapKeyValue(id2vertex,int,id,Vertex,v) {
		valid(v);
		assertx(v->id==id);
	} EndFor;
	ForMapKeyValue(id2face,int,id,Face,f) {
		valid(f);
		assertx(f->id==id);
	} EndFor;
	// Look over Vertices
	Set<Edge> sete;
	ForMeshVertex(*this,v1) {
		assertx(id2vertex.get(v1->id)==v1);
		Set<Vertex> set;
		ForStack(v1->ste,Edge,e) {
			// Check that Edges are valid
			assertx(e->next && e->prev && e->face && e->vert);
			assertx(e->prev->vert==v1);
			Vertex v2=e->vert;
			// Check that v2's are not duplicated
			assertx(set.add(v2));
			// Check that sym matches that in ste
			assertx(lookupEdge(v2,v1)==e->sym);
			// Check that Edge sym is valid
			if (e->sym) {
				assertx(e->sym->sym==e);
				assertx(e->sym->vert==e->prev->vert);
				assertx(e->sym->prev->vert==e->vert);
			}
			// Check that edges are unique
			assertx(sete.add(e));
			// Check that Faces reachable from Edges are valid
			valid(e->face);
			// Check that each Edge appears in its face
			int found=0;
			ForFaceHEdge(*this,e->face,ee) {
				if (ee==e) { found=1; break; }
			} EndFor;
			assertx(found);
		} EndFor;
	} EndFor;
	// Look over Faces
	ForMeshFace(*this,f) {
		assertx(id2face.get(f->id)==f);
		// Check that Face has valid erep
		assertx(erep(f)->face==f);
		valid(rep(erep(f)));
		Set<Vertex> set;
		ForFaceVertex(*this,f,v) {
			// Check Face contains no duplicate Vertices
			assertx(set.add(v));
			// Check Vertices in face are valid
			valid(v);
		} EndFor;
		// Check Edges in Faces are one-to-one with ste (1)
		ForFaceHEdge(*this,f,e) {
			assertx(e->face==f);
			assertx(sete.remove(e));
		} EndFor;
	} EndFor;
	// Check Edges in Faces are one-to-one with ste (2)
	assertx(sete.empty());
	// Check nedges
	NEST {
		int i=0;
		ForMeshEdge(*this,e) {
			(void)e; i++;
		} EndFor;
		assertx(nedges==i);
	}
}

int Mesh::isNice() const
{
	if (sdebug>=1) OK();
	ForMeshVertex(*this,v) {
		if (!isNice(v)) return 0;
	} EndFor;
	ForMeshFace(*this,f) {
		if (!isNice(f)) return 0;
	} EndFor;
	return 1;
}

void Mesh::valid(Vertex v) const
{
	assertx(v==id2vertex.get(vertexid(v)));
}

void Mesh::valid(Face f) const
{
	assertx(f==id2face.get(faceid(f)));
}

void Mesh::valid(Edge e) const
{
	assertx(e && e->next && e->prev);
	valid(vertex1(e));
	valid(vertex2(e));
	valid(face1(e));
	assertx(queryEdge(vertex1(e),vertex2(e))==e); // e==rep(e)
}

const char* Mesh::edgestring(Edge e) const
{
	valid(e);
	return hform("Edge(%d,%d)",vertexid(vertex1(e)),vertexid(vertex2(e)));
}

int Mesh::setDebug(int psdebug)
{
	int oldv=sdebug;
	sdebug=psdebug;
	return oldv;
}

Vertex Mesh::idvertex(int i) const
{
	return id2vertex.get(i);
}

int Mesh::vertexid(Vertex v) const
{
	return v->id;
}

Face Mesh::idface(int i) const
{
	return id2face.get(i);
}

int Mesh::faceid(Face f) const
{
	return f->id;
}

void Mesh::renumber()
{
	int id=1;
	ForMeshOrderedVertex(*this,v) {
		assertx(id2vertex.remove(v->id)==v);
		v->id=id;
		id2vertex.enter(id,v);
		id++;
	} EndFor;
	id=1;
	ForMeshOrderedFace(*this,f) {
		assertx(id2face.remove(f->id)==f);
		f->id=id;
		id2face.enter(id,f);
		id++;
	} EndFor;
}

//** Mesh protected

Edge Mesh::erep(Vertex v) const
{
//	if(v->ste.empty())
//            printf("\nMesh::erep - v->ste.empty()==true: vertex(%.15lg, %.15lg, %.15lg)\n\n", v->point[0], v->point[1], v->point[2]);
	return v->ste.empty()?0:v->ste.top()->prev;
}

Edge Mesh::erep(Face f) const
{
	return f->erep;
}

Edge Mesh::mostClwHEdge(Vertex v) const
{
	assertx(isNice(v));
	Edge e=erep(v);
	if (!e) return 0;
	for (Edge ef=e;;) {
		Edge en=clwHEdge(e);
		if (!en || en==ef) break;
		e=en;
	}
	return e;
}

Edge Mesh::mostCcwHEdge(Vertex v) const
{
	assertx(isNice(v));
	Edge e=erep(v);
	if (!e) return 0;
	for (Edge ef=e;;) {
		Edge en=ccwHEdge(e);
		if (!en || en==ef) break;
		e=en;
	}
	return e;
}

Edge Mesh::clwHEdge(Edge e) const
{
	return e->next->sym;
}

Edge Mesh::ccwHEdge(Edge e) const
{
	return e->sym?e->sym->prev:0;
}

Edge Mesh::unrepv1(Edge e, Vertex v) const
{
	if (vertex1(e)==v) return e;
	if (vertex2(e)==v) return e->sym;
	assertnever("Vertex not on Edge"); return 0;
}

Edge Mesh::unrepv2(Edge e, Vertex v) const
{
	if (vertex2(e)==v) return e;
	if (vertex1(e)==v) return e->sym;
	assertnever("Vertex not on Edge"); return 0;
}

Edge Mesh::lookupEdge(Vertex v1, Vertex v2) const
{
	assertx(v1 && v2);
	if (sdebug>=1) valid(v1),valid(v2);
	ForStack(v1->ste,Edge,e) {
		if (e->vert==v2) return e;
	} EndFor;
	return 0;
}

// Must have:		vert, prev->vert, face
// Defines:		sym
// Must define later:	next
void Mesh::enterHEdge(Edge e)
{
	Vertex v1=vertex1(e), v2=vertex2(e);
	Stack<Edge>& stack=v1->ste;
	if (sdebug>=1) assertx(!stack.contains(e));
	stack.push(e);
	e->sym=lookupEdge(v2,v1);
	if (e->sym) {
		assertx(!e->sym->sym);
		e->sym->sym=e;
		// insertion invalidates flag and info fields of sym
		e->sym->flags=0;
		delete e->sym->info,e->sym->info=0;
	} else {
		nedges++;
	}
}

// Must have:		vert, prev->vert, sym
// Must do later:	delete e
// Note: must be careful: could have e==erep(e->face) ! :
//   destroyFace() : deletes face -> ok.
//   substituteFaceVertex() : reintroduces edge -> ok.
void Mesh::removeHEdge(Edge e)
{
	Vertex v1=vertex1(e), v2=vertex2(e);
	if (sdebug>=1) assertx(e->sym==lookupEdge(v2,v1));
	if (e->sym) {
		e->sym->sym=0;
		// removal invalidates flag and info fields of sym
		e->sym->flags=0;
		delete e->sym->info,e->sym->info=0;
	} else {
		nedges--;
	}
	assertx(v1->ste.remove(e)); // slow, shucks
	e->flags=0;
	delete e->info,e->info=0;
}

//*** Info

const char* MeshInfo::getString() const { return 0; }
// virtual const char* MeshInfo::getString() const { return 0; }

void Mesh::setInfo(Vertex v, MeshInfo* newi)
{
	if (sdebug>=1) valid(v);
	delete v->info,v->info=newi;
}

void Mesh::setInfo(Face f, MeshInfo* newi)
{
	if (sdebug>=1) valid(f);
	delete f->info,f->info=newi;
}

void Mesh::setInfo(Edge e, MeshInfo* newi)
{
	if (sdebug>=1) valid(e);
	delete e->info,e->info=newi;
}

//*** ITERATORS
//** Mesh
//* MeshHEdgeIter

// JW, 2006-10-26: YUCK!  The below construction relies upon the fact that
//                 Stack<Edge> doesn't add anything to the binary
//                 representation of the object.  If it did, the iterator could
//                 be pointing off to something invalid...
MeshHEdgeIter::MeshHEdgeIter(const Mesh& m)
: mesh(m), mi(m.id2vertex),
  si(*(mi?&mi.value()->ste:(Stack<Edge> *)&Stack<Edge>::EMPTY)) { }

//* MeshVertexIter

MeshVertexIter::MeshVertexIter(const Mesh& m) : mi(m.id2vertex) { }

//* MeshOrderedVertexIter

MeshOrderedVertexIter::MeshOrderedVertexIter(const Mesh& m)
{
	ForMeshVertex(m,v) {
		int vi=m.vertexid(v);
		assertw1(vi<1e8); // precision available in double
		pq.enterUnsorted(v,vi);
	} EndFor;
	pq.sort();
}

//* MeshFaceIter

MeshFaceIter::MeshFaceIter(const Mesh& m): mi(m.id2face) { }

//* MeshOrderedFaceIter

MeshOrderedFaceIter::MeshOrderedFaceIter(const Mesh& m)
{
	ForMeshFace(m,f) {
		int fi=m.faceid(f);
		assertw1(fi<1e8); // precision available in double
		pq.enterUnsorted(f,fi);
	} EndFor;
	pq.sort();
}

//* MeshEdgeIter

MeshEdgeIter::MeshEdgeIter(const Mesh& m) : mesh(m), it(m) { }

//** Vertex
//* VertexHEdgeIter

VertexHEdgeIter::VertexHEdgeIter(const Mesh& m, Vertex v) : si(v->ste) { }

Edge VertexHEdgeIter::next()
{
	return si?si.next()->prev:0;
}

//* VertexVertexIter

VertexVertexIter::VertexVertexIter(const Mesh& m, Vertex v)
: it(m,v), extrav(0) { }

Vertex VertexVertexIter::next()
{
	if (extrav) { Vertex vr=extrav; extrav=0; return vr; }
	Edge e=it.next();
	if (!e) return 0;
	if (!e->next->sym) extrav=e->next->vert;
	return e->prev->vert;
}

//* VertexFaceIter

VertexFaceIter::VertexFaceIter(const Mesh& m, Vertex v) : it(m,v) { }

Face VertexFaceIter::next()
{
	Edge e=it.next();
	return e?e->face:0;
}

//* VertexEdgeIter

VertexEdgeIter::VertexEdgeIter(const Mesh& m, Vertex v)
: mesh(m), it(m,v), extrae(0) { }

Edge VertexEdgeIter::next()
{
	if (extrae) { Edge er=extrae; extrae=0; return er; }
	Edge e=it.next();
	if (!e) return 0;
	if (!e->next->sym) extrae=mesh.rep(e->next);
	return mesh.rep(e);
}

//** Face
//* FaceHEdge is inline

//* FaceVertexIter is inline

//* FaceFaceIter

FaceFaceIter::FaceFaceIter(const Mesh& m, Face f) : it(m,f) { }

Face FaceFaceIter::next()
{
	for (;;) {
		Edge e=it.next();
		if (!e) return 0;
		if (e->sym) return e->sym->face;
	}
}

//* FaceEdgeIter

FaceEdgeIter::FaceEdgeIter(const Mesh& m, Face f) : mesh(m), it(m,f) { }

Edge FaceEdgeIter::next()
{
	Edge e=it.next();
	return e?mesh.rep(e):0;
}

//** Edge
//* EdgeVertexIter

Vertex EdgeVertexIter::next()
{
	return i>=2?0:!i++?e->vert:e->prev->vert;
}

//* EdgeFaceIter

Face EdgeFaceIter::next()
{
	return i>=2?0:!i++?e->face:e->sym?e->sym->face:0;
}

//** VertexCcw
//* VertexCcwVertex

VertexCcwVertexIter::VertexCcwVertexIter(const Mesh& m, Vertex vp)
: mesh(m), v(vp), vf(mesh.mostClwVertex(v)), vc(vf) { }

Vertex VertexCcwVertexIter::next()
{
	Vertex vr=vc;
	if (vc) vc=mesh.ccwVertex(v,vc);
	if (vc==vf) vc=0;
	return vr;
}
