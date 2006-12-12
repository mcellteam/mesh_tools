// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Bbox.h"

Bbox::Bbox(const Point& pmin, const Point& pmax)
{
	p[0]=pmin;
	p[1]=pmax;
}

void Bbox::clear()
{
	const double big=1e30;
	p[0]=Point(big,big,big);
	p[1]=Point(-big,-big,-big);
}

void Bbox::infinite()
{
	const double big=1e30;
	p[0]=Point(-big,-big,-big);
	p[1]=Point(big,big,big);
}

void Bbox::takeunion(const Bbox& bb)
{
	for (int c=0;c<3;c++) {
		p[0][c]=min(p[0][c],bb[0][c]);
		p[1][c]=max(p[1][c],bb[1][c]);
	}
}

void Bbox::takeunion(const Point& pp)
{
	for (int c=0;c<3;c++) {
		p[0][c]=min(p[0][c],pp[c]);
		p[1][c]=max(p[1][c],pp[c]);
	}
}

void Bbox::intersect(const Bbox& bb)
{
	for (int c=0;c<3;c++) {
		p[0][c]=max(p[0][c],bb[0][c]);
		p[1][c]=min(p[1][c],bb[1][c]);
	}
}

int Bbox::inside(const Bbox& bb) const
{
	for (int c=0;c<3;c++) {
		if (p[0][c]<bb[0][c]) return 0;
		if (p[1][c]>bb[1][c]) return 0;
	}
	return 1;
}

Frame Bbox::getFrameToCube() const
{
        int c;
	Vector di=p[1]-p[0];
	double maxdi=0;
//	for (int c=0;c<3;c++) if (di[c]>maxdi) maxdi=di[c];
	for (c=0;c<3;c++) if (di[c]>maxdi) maxdi=di[c];
	assertx(maxdi);
	Vector center;
	for (c=0;c<3;c++) center[c]=(1-di[c]/maxdi)*.5;
	center[2]=0;	// objects lie at bottom of cube
	return (Frame::translation(-toVector(p[0]))*
		Frame::scaling(1/maxdi)*
		Frame::translation(center));
}

Frame Bbox::getFrameToSmallCube() const
{
	Frame f=getFrameToCube();
	f*=Frame::scaling(.8)*Frame::translation(.1,.1,.1);
	ForIndex(i,3) {
		if (fabs(f[i][i]-1.)<.05) f[i][i]=1.;
		if (fabs(f.p[i])<.05) f.p[i]=0.;
	} EndFor;
	return f;
}

