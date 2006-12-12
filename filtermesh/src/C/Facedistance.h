// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Facedistance_h
#define Facedistance_h

#include "Geometry.h"

// Find a lower bound on the distance between p and triangle (p1,p2,p3)
inline double LBDistPTri(const Point& p, const Point& p1,
			const Point& p2, const Point& p3);

extern double DistPointTriangle2(const Point& p, const Point& p1,
				const Point& p2, const Point& p3);

// Given point p and triangular face with vertices p1,p2,p3,
//  compute distance squared dis2, projection barycentric coordinates ba
//   (some may be negative if p projects outside the triangle),
//  convex barycentric coordinates of the closest point cba within
//  the triangle, and clp (the closest point in the triangle).
extern void ProjectPTri(const Point& p, const Point& p1,
			const Point& p2, const Point& p3,
			double* dis2, Bary* ba, Bary* cba, Point* clp);

//----------------------------------------------------------------------------

inline double LBDistPTri(const Point& p, const Point& p1,
			const Point& p2, const Point& p3)
{
	register double d=0;
	for (int c=0;c<3;c++) {
		register double v1=p1[c], mi=v1, ma=v1, v2=p2[c], v3=p3[c];
		if (v2<mi) mi=v2; else if (v2>ma) ma=v2;
		if (v3<mi) mi=v3; else if (v3>ma) ma=v3;
		register double v=p[c], a;
		if ((a=v-ma)>0) { if (a>d) d=a; }
		else if ((a=mi-v)>0) { if (a>d) d=a; }
	}
	return d;
}

#endif
