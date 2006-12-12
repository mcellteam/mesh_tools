// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#ifndef Mesh_h
#define Mesh_h

#include "Map.h"
#include "Set.h"
#include "Stack.h"
#include "Pqueue.h"
#include "Geometry.h"		// because of Point, too bad.
#include "Array.h"
#include "Pool.h"

// template<class T> class Array; // GNUG had problems in Recon.C and G3dX.C
class Random;

class MVertex;	typedef MVertex* Vertex;
class MFace;	typedef MFace* Face;
class MEdge;	typedef MEdge* Edge;
class MeshInfo;

// Mesh: a set of Vertices, Faces, and Edges and their topological relations.
// Properties:
//   (1a) - faces do not contain duplicate vertices
//   (1b) - mesh is orientable and edges are shared by 1 or 2 faces
//   (2a) - vertices are nice (contain at most 1 (possibly partial) face ring)
//   (2b) - faces are nice: (a,b,c) implies no (a,c,b)
//   (3a) - all faces are triangular
// a Mesh must always satisfy (1a) and (1b); such a mesh is called "legal"
// a Mesh is "nice" if it also satisfies (2a) and (2b)
// a Mesh is a "nice triangular mesh" if in addition it satisfies (3a)
//
// MVertex allocates space for Point, but Mesh does not use them; GMesh does.

class Mesh {
  public:
	Mesh();
	virtual ~Mesh();
	void clear();
	// copy does not carry flags or info
	void copy(const Mesh& m); // must be empty, no virtual constructors
	void copyVertices(const Mesh& m); // must be empty
// Raw manipulation functions, may lead to non-nice Meshes.
	// always legal
	Vertex createVertex();
	// die if degree(v)>0
	virtual void destroyVertex(Vertex v);
	// ret 0: illegal if has duplicate vertices or if existing edge
	Face createFace(Vertex va[], int nv);
	Face createFace(Vertex v1, Vertex v2, Vertex v3);
	// always legal
	virtual void destroyFace(Face f);
	// ret success (illegal if introduces duplicate vertex or edge)
	virtual int substituteFaceVertex(Face f, Vertex vold, Vertex vnew);
// Vertex
	int isNice(Vertex v) const;
	int degree(Vertex v) const; // == number of adjacent vertices/edges
	int numBoundaries(Vertex v) const; // 0/1 for nice vertex
	int isBoundary(Vertex v) const; // isNice(v), degree(v)>0
	Edge oppEdge(Vertex v, Face f) const; // isTriangle(f)
	Vertex oppVertex(Vertex v, Edge e) const;
	// isNice(v), ccw order
	int vertices(Vertex v, Array<Vertex>& va) const; // ret:is_ring
	int faces(Vertex v, Array<Face>& fa) const;	// ret:is_ring
	// all mostClw and mostCcw assert isNice(v)
	// move about vertices adjacent to a vertex
	Vertex mostClwVertex(Vertex v) const; // if !bnd,ret any; may ret 0
	Vertex mostCcwVertex(Vertex v) const; // if !bnd,ret any; may ret 0
	Vertex clwVertex(Vertex v, Vertex vext) const; // may ret 0
	Vertex ccwVertex(Vertex v, Vertex vext) const; // may ret 0
	// move about faces adjacent to a vertex
	Face mostClwFace(Vertex v) const; // if !bnd,ret any; may ret 0
	Face mostCcwFace(Vertex v) const; // if !bnd,ret any; may ret 0
	Face clwFace(Vertex v, Face f) const; // may ret 0
	Face ccwFace(Vertex v, Face f) const; // may ret 0
	// move about edges adjacent to a vertex
	Edge mostClwEdge(Vertex v) const; // if !bnd,ret any; may ret 0
	Edge mostCcwEdge(Vertex v) const; // if !bnd,ret any; may ret 0
	Edge clwEdge(Vertex v, Edge e) const; // may ret 0
	Edge ccwEdge(Vertex v, Edge e) const; // may ret 0
// Face
	int isNice(Face f) const;
	int numVertices(Face f) const;
	int numBoundaryEdges(Face f) const;
	int isTriangle(Face f) const;
	int isBoundary(Face f) const; // == has a boundary vertex
	Face oppFace(Face f, Edge e) const; // ret 0 if isBoundary(e)
	// ccw order
	void vertices(Face f, Array<Vertex>& va) const;
	void vertices(Face f, Vertex va[3]) const; // isTriangle(f)
	// move about edges adjacent to a face
	Edge clwEdge(Edge e, Face f) const;
	Edge ccwEdge(Edge e, Face f) const;
// Edge
	int isBoundary(Edge e) const;
	Vertex vertex1(Edge e) const;
	Vertex vertex2(Edge e) const;
	Face face1(Edge e) const;
	Face face2(Edge e) const; // ret 0 if isBoundary(e)
	Vertex sideVertex1(Edge e) const; // isTriangle(face1())
	Vertex sideVertex2(Edge e) const; // isTriangle(face2()); may ret 0
	Vertex sideVertex(Edge e, Face f) const; // isTriangle(f)
	Edge oppBoundary(Edge e, Vertex v) const; // isBoundary(e)
	Edge clwBoundary(Edge e) const;		  // isBoundary(e)
	Edge ccwBoundary(Edge e) const;		  // isBoundary(e)
// Other associations
	// obtain edge from vertices
	Edge queryEdge(Vertex v, Vertex w) const;
	Edge edge(Vertex v, Vertex w) const; // asserts it exists
	Edge orderedEdge(Vertex v1, Vertex v2) const; // asserts orientation!
	// obtain face from vertices
	Face face(Vertex v1, Vertex v2, Vertex v3) const; // ordered ccw
	// given vertex and adjacent face, get 2 related edges
	Edge clwEdge(Vertex v, Face f) const;
	Edge ccwEdge(Vertex v, Face f) const;
// Counting routines (fast)
	int numVertices() const;
	int numFaces() const;
	int numEdges() const;
// Random access (fast), assert there exist at least one
	Vertex randomVertex(Random& r) const;
	Face randomFace(Random& r) const;
	Edge randomEdge(Random& r) const; // unbiased for closed triangle mesh
// Flags
	// get selected flag(s)
	int gflag(int flagmask) const; // flag for whole mesh
	int flag(Vertex v, int flagmask) const;
	int flag(Face f, int flagmask) const;
	int flag(Edge e, int flagmask) const;
	// set or reset selected flag(s), return previous state of flag(s)
	int gmodFlag(int flagmask, int bool_flag);
	int modFlag(Vertex v, int flagmask, int bool_flag);
	int modFlag(Face f, int flagmask, int bool_flag);
	int modFlag(Edge e, int flagmask, int bool_flag);
// Info
	MeshInfo* info(Vertex v) const;
	void setInfo(Vertex v, MeshInfo* newi);
	MeshInfo* info(Face f) const;
	void setInfo(Face f, MeshInfo* newi);
	MeshInfo* info(Edge e) const;
	void setInfo(Edge e, MeshInfo* newi);
// Mesh operations (die unless triangular meshes!)
	// would collapse preserve a nice mesh?
	int niceEdgeCollapse(Edge e) const;
	// would collapse be legal?
	int legalEdgeCollapse(Edge e) const;
	// would edge swap be legal?  (legal implies nice here)
	int legalEdgeSwap(Edge e) const;
	// die if !legalEdgeCollapse(e)
	// remove (v1,?), (v2,?), v2, f1,[f2]; add edges (v1,?)
	virtual void collapseEdge(Edge e);
	// splitEdge(e) always legal
	// remove (v1,v2), (vo1,v2),[(vo2,v2)]
	// add vnew, (vnew,v2,vo1),[(vnew,vo2,v2)]
	// add (vnew,?), (vo1,v2), [(vo2,v2)]
	virtual Vertex splitEdge(Edge e, int vid=0);
	// ret 0 if !legalEdgeSwap(e)
	// remove (v1,v2),(v1,vo1),(v2,vo1),(v1,vo2),(v2,vo2),
	// add (vo1,vo2), (v1,vo1),(v2,vo1),(v1,vo2),(v2,vo2),
	// remove f1,f2, add (vo1,vo2,v2), (v1,vo2,vo1)
	virtual Edge swapEdge(Edge e);
// Misc
	void OK() const;
	int isNice() const;
	void valid(Vertex v) const;
	void valid(Face f) const;
	void valid(Edge e) const;
	const char* edgestring(Edge e) const;
	static int setDebug(int psdebug);
// Ids
	Vertex idvertex(int i) const;
	int vertexid(Vertex v) const;
	Face idface(int i) const;
	int faceid(Face f) const;
	void renumber();	// renumber vertices and faces
// Discouraged:
	// next 2 die if index is already used
	virtual Vertex createVertexI(int id);
	virtual Face createFaceI(int id, Vertex va[], int nv); // may ret 0
  protected:
	static int sdebug;	// 0=no, 1=min, 2=max
    friend class MeshHEdgeIter;
    friend class MeshVertexIter;
    friend class MeshFaceIter;
    friend class MeshEdgeIter;
    friend class VertexHEdgeIter;
    friend class VertexEdgeIter;
    friend class FaceHEdgeIter;
    friend class FaceEdgeIter;
  private:
	int iflags;
	Map<int,Vertex> id2vertex; // also act as set of vertices and faces
	Map<int,Face> id2face;
	int vertexnum;		// id to assign to next new vertex
	int facenum;		// id to assign to next new face
	int nedges;
	Edge erep(Vertex v) const; // not necessarily a rep edge
	Edge erep(Face v) const;   // not necessarily a rep edge
	// all HEdge point toward v
	Edge mostClwHEdge(Vertex v) const; // isNice(v), may ret 0
	Edge mostCcwHEdge(Vertex v) const; // isNice(v), may ret 0
	Edge clwHEdge(Edge e) const;	   // may ret 0
	Edge ccwHEdge(Edge e) const;	   // may ret 0
	Edge rep(Edge e) const;
	Edge unrepv1(Edge e, Vertex v) const;
	Edge unrepv2(Edge e, Vertex v) const;
	Edge lookupEdge(Vertex v1, Vertex v2) const; // not always a rep()
	void enterHEdge(Edge e); // must have: vert, prev->vert, face
	void removeHEdge(Edge e); // call 'delete e' later
	DISABLECOPY(Mesh);
};

//*** Iterators.
// Iterators can crash if continued after any change in the Mesh.
// HEdge iterators should not be used by the general public.

// Mesh iterators do not have order
#define ForMeshHEdge(m,e) \
{ Edge e; for (MeshEdgeIter zz(m); (e=zz.next()) != 0;) {
#define ForMeshVertex(m,v) \
{ Vertex v; for (MeshVertexIter zz(m); (v=zz.next()) != 0;) {
#define ForMeshOrderedVertex(m,v) \
{ Vertex v; for (MeshOrderedVertexIter zz(m); (v=zz.next()) != 0;) {
#define ForMeshFace(m,f) \
{ Face f; for (MeshFaceIter zz(m); (f=zz.next()) != 0;) {
#define ForMeshOrderedFace(m,f) \
{ Face f; for (MeshOrderedFaceIter zz(m); (f=zz.next()) != 0;) {
#define ForMeshEdge(m,e) \
{ Edge e; for (MeshEdgeIter zz(m); (e=zz.next()) != 0;) {

// Vertex iterators do not specify order! Work correctly on non-nice vertices.
#define ForVertexHEdge(m,vv,e) \
{ Edge e; for (VertexHEdgeIter zz(m,vv); (e=zz.next()) != 0;) {
#define ForVertexVertex(m,vv,v) \
{ Vertex v; for (VertexVertexIter zz(m,vv); (v=zz.next()) != 0;) {
#define ForVertexFace(m,vv,f) \
{ Face f; for (VertexFaceIter zz(m,vv); (f=zz.next()) != 0;) {
#define ForVertexEdge(m,vv,e) \
{ Edge e; for (VertexEdgeIter zz(m,vv); (e=zz.next()) != 0;) {

// Face iterators all go CCW
#define ForFaceHEdge(m,ff,e) \
{ Edge e; for (FaceHEdgeIter zz(m,ff); (e=zz.next()) != 0;) {
#define ForFaceVertex(m,ff,v) \
{ Vertex v; for (FaceVertexIter zz(m,ff); (v=zz.next()) != 0;) {
#define ForFaceFace(m,ff,f) \
{ Face f; for (FaceFaceIter zz(m,ff); (f=zz.next()) != 0;) {
#define ForFaceEdge(m,ff,e) \
{ Edge e; for (FaceEdgeIter zz(m,ff); (e=zz.next()) != 0;) {
	
#define ForEdgeVertex(m,ee,v) \
{ Vertex v; for (EdgeVertexIter zz(m,ee); (v=zz.next()) != 0;) {
#define ForEdgeFace(m,ee,f) \
{ Face f; for (EdgeFaceIter zz(m,ee); (f=zz.next()) != 0;) {

#define ForVertexCcwVertex(m,vv,v) \
{ Vertex v; for (VertexCcwVertexIter zz(m,vv); (v=zz.next()) != 0;) {

//-----------------------------------------------------------------------------

class MeshInfo {		// base class from which to derive
  public:
	MeshInfo() { }
	virtual ~MeshInfo() { }
	// string representation to output upon mesh write (0=none)
	virtual const char* getString() const;
  protected:
	// disable here to be safe, may be defined for derived classes
	DISABLECOPY(MeshInfo);
};

class MVertex {
    friend class Mesh;
    friend class GMesh;
    friend class MeshHEdgeIter;
    friend class VertexHEdgeIter;
  private:
	Stack<Edge> ste;
	int id;
	int flags;
	MeshInfo* info;
	Point point;
	MVertex() : ste(), id(0), flags(0), info(0), point() { }
	~MVertex() { delete info; }
	POOLALLOCATION(MVertex);
	DISABLECOPY(MVertex);
};

class MFace {
    friend class Mesh;
    friend class GMesh;
  private:
	Edge erep;
	int id;
	int flags;
	MeshInfo* info;
	MFace() : erep(), id(0), flags(0), info(0) { }
	~MFace() { delete info; }
	POOLALLOCATION(MFace);
	DISABLECOPY(MFace);
};

class MEdge {
    friend class Mesh;
    friend class GMesh;
    friend class VertexHEdgeIter;
    friend class VertexVertexIter;
    friend class VertexFaceIter;
    friend class VertexEdgeIter;
    friend class FaceHEdgeIter;
    friend class FaceVertexIter;
    friend class FaceFaceIter;
    friend class EdgeVertexIter;
    friend class EdgeFaceIter;
  private:
	Edge prev;		// previous edge in ring around face
	Edge next;		// next edge in ring around face
	Edge sym;		// pointer to symmetric half-edge
	Vertex vert;		// vertex to which edge is pointing
	Face face;		// face to which edge belongs
	int flags;
	MeshInfo* info;
	MEdge() : prev(), next(), sym(), vert(), face(), flags(0), info(0) { }
	~MEdge() { delete info; }
	POOLALLOCATION(MEdge);
	DISABLECOPY(MEdge);
};

//* Mesh

class MeshHEdgeIter {
  public:
	MeshHEdgeIter(const Mesh& m);
	Edge next();
  private:
	const Mesh& mesh;
	MapIter<int,Vertex> mi;
	StackIter<Edge> si;
};

class MeshVertexIter {
  public:
	MeshVertexIter(const Mesh& m);
	Vertex next();
  private:
	MapIter<int,Vertex> mi;
};

class MeshOrderedVertexIter {
  public:
	MeshOrderedVertexIter(const Mesh& m);
	Vertex next();
  private:
	Pqueue<Vertex> pq;
};

class MeshFaceIter {
  public:
	MeshFaceIter(const Mesh& m);
	Face next();
  private:
	MapIter<int,Face> mi;
};

class MeshOrderedFaceIter {
  public:
	MeshOrderedFaceIter(const Mesh& m);
	Face next();
  private:
	Pqueue<Face> pq;
};

class MeshEdgeIter {
  public:
	MeshEdgeIter(const Mesh& m);
	Edge next();
  private:
	const Mesh& mesh;
	MeshHEdgeIter it;
};

//* Vertex

class VertexHEdgeIter {
  public:
	VertexHEdgeIter(const Mesh& m, Vertex v);
	Edge next();
  private:
	StackIter<Edge> si;
};

class VertexVertexIter {
  public:
	VertexVertexIter(const Mesh& m, Vertex v);
	Vertex next();
  private:
	VertexHEdgeIter it;
	Vertex extrav;
};

class VertexFaceIter {
  public:
	VertexFaceIter(const Mesh& m, Vertex v);
	Face next();
  private:
	VertexHEdgeIter it;
};

class VertexEdgeIter {
  public:
	VertexEdgeIter(const Mesh& m, Vertex v);
	Edge next();
  private:
	const Mesh& mesh;
	VertexHEdgeIter it;
	Edge extrae;
};

//* Face

class FaceHEdgeIter {
  public:
	FaceHEdgeIter(const Mesh& m, Face f);
	Edge next();
  private:
	Edge ef;
	Edge e;
};

class FaceVertexIter {
  public:
	FaceVertexIter(const Mesh& m, Face f);
	Vertex next();
  private:
	FaceHEdgeIter it;
};

class FaceFaceIter {
  public:
	FaceFaceIter(const Mesh& m, Face f);
	Face next();
  private:
	FaceHEdgeIter it;
};

class FaceEdgeIter {
  public:
	FaceEdgeIter(const Mesh& m, Face f);
	Edge next();
  private:
	const Mesh& mesh;
	FaceHEdgeIter it;
};

//* Edge

class EdgeVertexIter {
  public:
	EdgeVertexIter(const Mesh&, Edge ee) : e(ee), i(0) { }
	Vertex next();
  private:
	Edge e;
	int i;
};

class EdgeFaceIter {
  public:
	EdgeFaceIter(const Mesh&, Edge ee) : e(ee), i(0) { }
	Face next();
  private:
	Edge e;
	int i;
};

class VertexCcwVertexIter {
  public:
	VertexCcwVertexIter(const Mesh& m, Vertex vp);
	Vertex next();
  private:
	const Mesh& mesh;
	Vertex v;
	Vertex vf;
	Vertex vc;
};

//*** inlines

inline Edge Mesh::rep(Edge e) const
{
	return e->sym && e->vert->id<e->sym->vert->id ? e->sym : e;
}

//** Mesh iter

inline Edge MeshHEdgeIter::next()
{
	if (si) return si.next();
	if (!mi) return 0;
	for (;;) {
		mi.next();
		if (!mi) return 0;
		si.reinit(mi.value()->ste);
		if (si) return si.next();
	}
}

inline Vertex MeshVertexIter::next()
{
	Vertex vr=0;
	if (mi) vr=mi.value(),mi.next();
	return vr;
}

inline Vertex MeshOrderedVertexIter::next()
{
	return pq.empty()?0:pq.removemin();
}

inline Face MeshFaceIter::next()
{
	Face fr=0;
	if (mi) fr=mi.value(),mi.next();
	return fr;
}

inline Face MeshOrderedFaceIter::next()
{
	return pq.empty()?0:pq.removemin();
}

inline Edge MeshEdgeIter::next()
{
	Edge er;
	for (;;) {
		er=it.next();
		if (!er || er==mesh.rep(er)) break;
	}
	return er;
}

//** Face iter

inline FaceHEdgeIter::FaceHEdgeIter(const Mesh& m, Face f)
{
	e=ef=m.erep(f);
}

inline Edge FaceHEdgeIter::next()
{
	if (!e) return 0;
	Edge er=e;
	e=e->next;
	if (e==ef) e=0;
	return er;
}

inline FaceVertexIter::FaceVertexIter(const Mesh& m, Face f) : it(m,f) { }

inline Vertex FaceVertexIter::next()
{
	Edge e=it.next();
	return e?e->vert:0;
}

//*** Info

inline MeshInfo* Mesh::info(Vertex v) const { return v->info; }
inline MeshInfo* Mesh::info(Face f) const { return f->info; }
inline MeshInfo* Mesh::info(Edge e) const { return e->info; }

//*** Flags

inline int Mesh::gflag(int flm) const { return iflags&flm; }
inline int Mesh::flag(Vertex v, int flm) const { return v->flags&flm; }
inline int Mesh::flag(Face f, int flm) const { return f->flags&flm; }
inline int Mesh::flag(Edge e, int flm) const { return e->flags&flm; }

inline int Mesh::gmodFlag(int flm, int bool_flag) {
	int st=gflag(flm);
	if (bool_flag) { if (st!=flm) iflags|=flm; }
	else { if (st) iflags&=~flm; }
	return st;
}
inline int Mesh::modFlag(Vertex v, int flm, int bool_flag) {
	int st=flag(v,flm);
	if (bool_flag) { if (st!=flm) v->flags|=flm; }
	else { if (st) v->flags&=~flm; }
	return st;
}
inline int Mesh::modFlag(Face f, int flm, int bool_flag) {
	int st=flag(f,flm);
	if (bool_flag) { if (st!=flm) f->flags|=flm; }
	else { if (st) f->flags&=~flm; }
	return st;
}
inline int Mesh::modFlag(Edge e, int flm, int bool_flag) {
	int st=flag(e,flm);
	if (bool_flag) { if (st!=flm) e->flags|=flm; }
	else { if (st) e->flags&=~flm; }
	return st;
}

#endif
