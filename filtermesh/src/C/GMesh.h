// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Gmesh_h
#define Gmesh_h

#include "Mesh.h"
#include "Geometry.h"
#include <iostream>

class Polygon; class WA3dStream; class A3dVertexColor;

class GMesh : public Mesh {
  public:
	GMesh();
	~GMesh();
	// copy carries flags but not info
	void copy(const GMesh& m); // must be empty, copies Point
	void copyVertices(const GMesh& m); // must be empty
// Geometry
	const Point& point(Vertex v) const;
	void setPoint(Vertex v, const Point& p);
	void polygon(Face f, Polygon& polygon) const;
	void points(Face f, Point pa[3]) const;
	double length2(Edge e) const;
	double length(Edge e) const;
// standard I/O for my meshes (see format below)
	void read(std::istream& is);	// read a whole mesh, discard comments
	void readline(char* s); // no '\n' required
	static int recognizeLine(const char* s);
	void write(std::ostream& os) const;
	void write(WA3dStream& oa3d, const A3dVertexColor& col) const;
	std::ostream* recordChanges(std::ostream* pos); // pos may be 0, ret old
// Flag bits
	enum { ALL=~0 };	// all flags
	// Predefined Vertex,Face,Edge flag bits
	// first 2 are parsed when reading
	enum { CUSPV=1<<30 };	// "cusp" on Vertex
	enum { SHARPE=1<<30 };	// "sharp" on Edge
	enum { ELEMOK=1<<29 };	// used to tag Vertex,Face,Edge
	// Predefined Mesh flag
	enum { MESHOK=1<<30 };
// Override Mesh virtuals
	void destroyVertex(Vertex v);
	void destroyFace(Face f);
	int substituteFaceVertex(Face f, Vertex vold, Vertex vnew);
	// do appropriate actions with geometry and flags
	void collapseEdge(Edge e);
	Vertex splitEdge(Edge e, int id=0);
	Edge swapEdge(Edge e);
// Discouraged:
	Vertex createVertexI(int id);
	Face createFaceI(int id, Vertex va[], int nv);
  private:
	std::ostream* ios;		// for recordChanges
};

// I/O Mesh Format (Vertices and Faces must fit on one line)
//   (vertex numbers begin with 1)
//   Vertex vi  x y z [{other_info}]
//   Face fi  vi1 vi2 ... vin [{other_info}]
//   MVertex vi newx newy newz
//   Ecol v1 v2
//   Eswa v1 v2
//   Espl v1 v2 vnew
// Example:
//   Vertex 1  1.5e2 0 1.5 {normal=(0,1,0)}
//   Vertex 2  0 1.5 0
//   Face 1  1 2 3
//   Face 2  2 3 4 5 {color=red, phong=2}
//  fi may be zero, in which case a number is assigned

class MeshInfoString : public MeshInfo {
  public:
	MeshInfoString(const char* ps);	// uses copy of string
	~MeshInfoString();
	const char* getString() const;
	void setString(const char* ps);	// uses copy of string
  private:
	char* s;
};

//----------------------------------------------------------------------

inline const Point& GMesh::point(Vertex v) const { return v->point; }

#endif
