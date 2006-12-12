// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef GeomOp_h
#define GeomOp_h

#if 0
{
	double ang[3];
	FrameToEuclideanAngles(f,ang);
	EuclideanAnglesToFrame(ang,f);
}
#endif

#include "Geometry.h"

//*** Radii

// Compute the circumscribed radius of the 3 points p0, p1, p2.
extern double CircumRadius(const Point& p0, const Point& p1, const Point& p2);

// Compute the inscribed radius of the 3 points p0, p1, p2.
extern double InRadius(const Point& p0, const Point& p1, const Point& p2);

//*** Misc

extern double DihedralAngleCos(const Point& p1, const Point& p2,
			      const Point& po1, const Point& po2);

extern double SolidAngle(const Point& p, const Point pa[], int np);

//*** Frames and Euclidean angles

// Compute Euclidean angles of f.
extern void FrameToEuclideanAngles(const Frame& f, double ang[3]);

// Modify f by setting v[0],v[1],v[2] according to Euclidean angles.
// f.p is unchanged.
extern void EuclideanAnglesToFrame(const double ang[3], Frame& f);

// Modify f so that its x axis points towards p and its y axis is vertical;
// f.p is ignored and unchanged.
extern void FrameAimAt(Frame& f, const Vector& v);

// Modify f so that its y axis lies in the xy plane
extern void FrameMakeLevel(Frame& f);

// Modify f so that its x and y axes lie in the xy plane.
extern void FrameMakeHoriz(Frame& f);

#endif
