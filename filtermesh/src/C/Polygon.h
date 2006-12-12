// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Polygon_h
#define Polygon_h

#include "Geometry.h"
#include "Array.h"
#include "Pool.h"

#include <iostream>

class A3dElem;
class Bbox;

class Polygon : public Array<Point> {
  public:
	Polygon(int nv=0);
	Polygon(const A3dElem& el);
	~Polygon();
	void copy(const Polygon& poly);
	void getbbox(Bbox& bb) const;
	Vector getnormaldir() const; // non-normalized normal
	Vector getnormal() const; // user should check !iszero()
	double getplanec(const Vector& pnor) const;
	double gettolerance(const Vector& pnor, double d) const;
	double getarea() const;
	// can call centroid() on a Polygon!
	// Finds intersection of polygon with halfspace in +hn direction
	int intersectHyperplane(const Point& hp,
				const Vector& hn); // ret: is_modified
	int intersectBbox(const Bbox& bb);	   // ret: is_modified
	int intersectSegment(const Point& p1, const Point& p2,
			     Point& pint) const; // ret: is_intersection
	// Intersect with plane defined (planenor,planed,planetol).
	// Report intersection as array of points pa.
	void intersectPlane(const Vector& polynor, const Vector& planenor,
			    double planed, double planetol,
			    Array<Point>& pa) const;
	// Intersect 2 polygons.  np is number of points allocated and returned
	// np is always even.  pa[2i]..pa[2i+1] are segments
	friend void IntersectPolyPoly(const Polygon& p1, const Polygon& p2,
				      Array<Point>& pa);
	int pointinside(const Vector& pnor,
			const Point& point) const; // ret: point_is_inside
	int isconvex() const;
	friend std::ostream& operator<<(std::ostream& s, const Polygon& poly);
	POOLALLOCATION(Polygon);
};

extern int IntersectPlaneSegment(const Vector& normal, double d,
				 const Point& p1, const Point& p2,
				 Point& pint); // ret is_intersection

// Find a vector perpendicular to v (any such vector)
extern Vector OrthogonalVector(const Vector& v);

// Assign a predictable direction to v
extern void VectorStandardDirection(Vector& v);

//----------------------------------------------------------------------------

inline Polygon::Polygon(int nv) : Array<Point>(nv) { }

#endif
