// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef PolygonSpatial_h
#define PolygonSpatial_h

#include "Spatial.h"
#include "Polygon.h"

class PolygonSpatial : public ObjectSpatial {
  public:
	PolygonSpatial(int pgn);
	~PolygonSpatial();
	// clear() inherited from ObjectSpatial, does not delete Polygons!
	// deleteclear() must be implemented by user
	void enter(const Polygon* poly); // poly is not copied!
	int firstAlongSegment(const Point& p1, const Point& p2,
			      const Polygon*& poly, Point& pint) const;
};

#endif
