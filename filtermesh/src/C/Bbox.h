// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Bbox_h
#define Bbox_h

#include "Geometry.h"

class Bbox {
  public:
	Bbox() { }
	Bbox(const Point& pmin, const Point& pmax);
	~Bbox() { }
	void clear();
	void infinite();
	Point& operator[](int i) { return p[i]; }
	const Point& operator[](int i) const { return p[i]; }
	void takeunion(const Bbox& bb);
	void takeunion(const Point& pp);
	void intersect(const Bbox& bb);
	int inside(const Bbox& bb) const;
	// uniform scaling into unit cube, centered on x & y, rest at z=0
	Frame getFrameToCube() const;
	Frame getFrameToSmallCube() const;
  private:
	Point p[2];
	// shallow copy ok
};

#endif
