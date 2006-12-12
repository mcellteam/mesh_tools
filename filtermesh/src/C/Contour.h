// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Contour_h
#define Contour_h

#include "Geometry.h"
#include "GMesh.h"
#include "Polygon.h"
#include "Map.h"
#include "Queue.h"

#include <iostream>

class ContourNode;

// Extract the zeroset of a scalar function defined in the unit cube.
// Can do either 3D or 2D contours.
// If the zero-curve of a 2D function is specified, all x coordinates must be
//  zero everywhere.
class Contour {
  public:
	Contour(int pis3D, double (*pfeval)(const Point& p),
		int pgn, std::ostream* pos=&std::cerr);
	~Contour();
	void setOutputContour(void (*pfoutputcontour)(const Polygon& poly));
	void setOutputBorder(void (*pfoutputborder)(const Polygon& poly));
	void setOutputMesh(GMesh& pmesh); // only for 3D
	void bigMeshFaces();
	void setVertexTolerance(double tol);
	int marchFrom(const Point& startp); // ret number of new cubes visited:
					    // 0=revisit, 1=no_surf, >1=new
	// call marchFrom() on all cells near startp
	int marchNear(const Point& startp); // ret number of new cubes visited:
	static const double UNDEF;	    // undefined scalar VALUE
  private:
	int is3D;
	double (*feval)(const Point& p);
	int gn;
	std::ostream* ios;
	double gnf;		// double(gn)
	double gni;		// double(1/gn);
	void (*foutputcontour)(const Polygon& poly);
	void (*foutputborder)(const Polygon& poly);
	GMesh* mesh;
	int bigmeshfaces;
	double vertextol;	// note: 0 is special: infinite tolerance
	Map<int,ContourNode*> m; // vertex index -> ContourNode*
	Queue<int> queue;	 // cubes queued to be visited
	int ncvisited;
	int ncundef;
	int ncnothing;
	int nvevaled;
	int nvzero;
	int nvundef;
	int nedegen;
	Map<int,Vertex> mapenev; // for mesh, encoded edge -> Vertex
	int encode(int ci[3]) const;
	void decode(int en, int ci[3]) const;
	int inbound(int ci[3]) const;
	int marchfrom(int cc[3]);
	void considercube(int encube);
 	void pushneighbors(int d, const int cc[3], ContourNode* na[2][2][2]);
	void contourcube(ContourNode* na[2][2][2]);
	void contourtetrahedron(ContourNode* n4[4]);
	void outputtriangle(ContourNode* n3[3][2]);
	Point computepoint(const Point& pp, const Point& pn,
			  double vp, double vn);
	void contoursquare(ContourNode* na[2][2][2]);
	void contourtriangle(ContourNode* n3[3]);
	void outputline(ContourNode* n2[2][2]);
	void cubemesh(ContourNode* na[2][2][2]);
	void examineface(int d, int v, ContourNode* na[2][2][2],
			 Map<Vertex,Vertex>& mapsucc);
	Vertex getvertexonedge(ContourNode* n1, ContourNode* n2);
	DISABLECOPY(Contour);
};

#endif
