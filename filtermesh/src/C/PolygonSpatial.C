// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "PolygonSpatial.h"
#include "Facedistance.h"
#include "Bbox.h"

static double polygonapproxdistance2(const Point& p, Univ id)
{
	Polygon& poly=*Conv<Polygon*>::d(id);
	assertx(poly.num()==3);
	return square(LBDistPTri(p,poly[0],poly[1],poly[2]));
}

static double polygondistance2(const Point& p, Univ id)
{
	Polygon& poly=*Conv<Polygon*>::d(id);
	assertx(poly.num()==3);
	return DistPointTriangle2(p,poly[0],poly[1],poly[2]);
}

PolygonSpatial::PolygonSpatial(int pgn)
: ObjectSpatial(pgn,polygonapproxdistance2,polygondistance2) { }

PolygonSpatial::~PolygonSpatial() { }


static struct {
	const Polygon* polyorigp;
	Polygon poly;
	Bbox bbox;
} senter;

static int polygoninbbox(const Bbox& bb)
{
	for (int c=0;c<3;c++)
		if (senter.bbox[0][c]>bb[1][c] ||
		    senter.bbox[1][c]<bb[0][c]) return 0;
	int modif=senter.poly.intersectBbox(bb);
	int ret=senter.poly.num()>0;
	if (modif) senter.poly.copy(*senter.polyorigp);
	return ret;
}

void PolygonSpatial::enter(const Polygon* poly)
{
	assertx(poly->num()==3);
	senter.polyorigp=poly;
	senter.poly.copy(*senter.polyorigp);
	senter.poly.getbbox(senter.bbox);
	ObjectSpatial::enter(Conv<const Polygon*>::e(poly),
			     (*poly)[0],polygoninbbox);
}


static struct Ssinter {
	const Point* p1;
	const Point* p2;
	Vector vray;
	Point pint;
	int foundint;
	const Polygon* poly;
	double fmin;
} sinter;

static int testpolygonwithray(Univ id)
{
	const Polygon* poly=Conv<Polygon*>::d(id);
	Point pint;
	if (!poly->intersectSegment(*sinter.p1,*sinter.p2,pint)) return 0;
	double f=dot(pint-*sinter.p1,sinter.vray);
	if (!sinter.foundint || f<sinter.fmin)
		sinter.fmin=f,sinter.pint=pint,sinter.poly=poly;
	sinter.foundint=1;
	return 1;
}

int PolygonSpatial::firstAlongSegment(const Point& p1, const Point& p2,
				      const Polygon*& poly, Point& pint) const
{
	sinter.p1=&p1;
	sinter.p2=&p2;
	sinter.vray=p2-p1;
	sinter.foundint=0;
	searchsegment(p1,p2,testpolygonwithray);
	poly=sinter.poly;
	pint=sinter.pint;
	return sinter.foundint;
}
