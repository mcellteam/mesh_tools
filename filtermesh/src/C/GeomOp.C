// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "GeomOp.h"

//*** Radii

// Deduced from book: Coxeter "Geometry".
double CircumRadius(const Point& p0, const Point& p1, const Point& p2)
{
	double a=dist(p0,p1),b=dist(p1,p2),c=dist(p2,p0);
	double s=(a+b+c)*.5;
	double d2=s*(s-a)*(s-b)*(s-c);
	return assertw1(d2>0)? 1e10 : a*b*c*.25/mysqrt(d2);
}

// From Coxeter "Intro to Geometry, Second ed.", page 12, equation  1.531
double InRadius(const Point& p0, const Point& p1, const Point& p2)
{
	// r=d/s
	// d=sqrt(s(s-a)(s-b)(s-c))
	// s=(a+b+c)/2
	double a=dist(p0,p1),b=dist(p1,p2),c=dist(p2,p0);
	double s=(a+b+c)*.5;
	double d2=s*(s-a)*(s-b)*(s-c);
	return assertw1(d2>0)? 0 : mysqrt(d2)/s;
}

//*** Misc

double DihedralAngleCos(const Point& p1, const Point& p2,
		       const Point& po1, const Point& po2)
{
	Vector ves1=cross(p1,p2,po1);
	if (!ves1.normalize()) return -2;
	Vector ves2=cross(p1,po2,p2);
	if (!ves2.normalize()) return -2;
	double d=dot(ves1,ves2);
	if (d<-1) d=-1;
	if (d> 1) d= 1;
	return d;
}

double SolidAngle(const Point& p, const Point pa[], int np)
{
	// solid angle: fraction area covered on sphere centered about p
	//  maximum=4*PI
	// idea: Gauss-Bonnet theorem:
	// integral of curvature + line integral + exterior angles = 2*PI
	// integral of curvature on unit sphere is equal to area
	// line integral along geodesics (great circles) is zero
	// So, solid angle=2*PI- sum of exterior angles on unit sphere
	double sumang=0;
	for (int i=0;i<np;i++) {
		int ip=(i-1+np)%np, in=(i+1)%np;
		Vector top=pa[i]-p;
		if (assertw1(top.normalize())) continue;
		Vector v1=pa[i]-pa[ip];
		v1-=top*dot(v1,top);
		if (assertw1(v1.normalize())) continue;
		Vector v2=pa[in]-pa[i];
		v2-=top*dot(v2,top);
		if (assertw1(v2.normalize())) continue;
		double vcos=dot(v1,v2);
		double vsin=dot(cross(v1,v2),top);
		double ang=atan2(vsin,vcos);
		sumang+=ang;
	}
	double solidang=2*PI-sumang;
	return solidang;
}

//*** Frames and Euclidean angles

inline double myatan2(double y, double x)
{
	return !y&&!x? 0 : atan2(y,x);
}

// angle 0 is z axis (alpha)
// angle 1 is y axis (beta)
// angle 2 is x axis (phi)
// To extract angles, look at tostdf:
//	cos(b)cos(a)	cos(b)sin(a)	-sin(b)
//	?		?		sin(p)cos(b)
//	?		?		cos(p)cos(b)

void FrameToEuclideanAngles(const Frame& f, double ang[3])
{
	ang[0]=myatan2(f[0][1],f[0][0]);
	ang[1]=myatan2(-f[0][2],hypot(f[0][0],f[0][1]));
	ang[2]=myatan2(f[1][2]/hypot(hypot(f[1][0],f[1][1]),f[1][2]),
		       f[2][2]/hypot(hypot(f[2][0],f[2][1]),f[2][2]));
}

void EuclideanAnglesToFrame(const double ang[3], Frame& f)
{
        int c;
	Frame fr=Frame::identity();
//	for (int c=0;c<3;c++)
	for (c=0;c<3;c++)
		fr[c][c]=mag(f.v[c]);
	for (c=0;c<3;c++)
		fr=fr*Frame::rotation(c,ang[2-c]);
	for (c=0;c<3;c++)
		f[c]=fr[c];
}

void FrameAimAt(Frame& f, const Vector& v)
{
	double ang[3];
	ang[0]=myatan2(v[1],v[0]);
	ang[1]=myatan2(-v[2],hypot(v[0],v[1]));
	ang[2]=0;
	EuclideanAnglesToFrame(ang,f);
}

void FrameMakeLevel(Frame& f)
{
	double ang[3];
	FrameToEuclideanAngles(f,ang);
	ang[2]=0;
	EuclideanAnglesToFrame(ang,f);
}

void FrameMakeHoriz(Frame& f)
{
	double ang[3];
	FrameToEuclideanAngles(f,ang);
	ang[1]=0;
	ang[2]=0;
	EuclideanAnglesToFrame(ang,f);
}
