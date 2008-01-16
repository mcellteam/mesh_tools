// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

// 06/12/92     started
// 06/13/92     GNU g++2.0 has trouble with:
//   class Vector { ... static const Vector zero; ... };

#include "Hh.h"
#include "Geometry.h"
#include "Homogeneous.h"	// for Frame <> Matrix4 conversions, centroid
#include "Array.h"

using std::ostream;
int invert(const Frame& fi, Frame& fo);
Frame transpose(const Frame& f);
void transpose(const Frame& fi, Frame& fo);

//*** Affine

ALLOCATEPOOL(Affine);

ostream& operator<<(ostream& s, const Affine& a)
{
	return s << "Affine(" << a[0] << "," << a[1] << "," << a[2] << ")\n";
}

//*** Vector

ostream& operator<<(ostream& s, const Vector& v)
{
	return s << "Vector(" << v[0] << "," << v[1] << "," << v[2] << ")\n";
}

Vector operator*(const Vector& v, const Frame& f)
{
	register double xi=v[0], yi=v[1], zi=v[2];
	return Vector(xi*f[0][0]+yi*f[1][0]+zi*f[2][0],
		      xi*f[0][1]+yi*f[1][1]+zi*f[2][1],
		      xi*f[0][2]+yi*f[1][2]+zi*f[2][2]);
}

Vector& Vector::operator*=(const Frame& f)
{
	return *this=*this*f;
}


Vector operator*(const Frame& f, const Vector& n)
{
	register double xi=n[0], yi=n[1], zi=n[2];
	return Vector(f[0][0]*xi+f[0][1]*yi+f[0][2]*zi,
		      f[1][0]*xi+f[1][1]*yi+f[1][2]*zi,
		      f[2][0]*xi+f[2][1]*yi+f[2][2]*zi);
}

int normalize(const Vector& vi, Vector& vo)
{
	register double sum2=vi[0]*vi[0]+vi[1]*vi[1]+vi[2]*vi[2];
	if (!sum2) return 0;
	register double fac=1/mysqrt(sum2);
	vo[0]=vi[0]*fac; vo[1]=vi[1]*fac; vo[2]=vi[2]*fac;
	return 1;
}

int Vector::normalize()
{
	return ::normalize(*this,*this);
}

Vector oknormalize(const Vector& v)
{
	Vector vr;
	if (normalize(v,vr)) return vr;
	return Vector(0,0,0);
}

int compare(const Vector& v1, const Vector& v2)
{
	double d;
	if ((d=v1[0]-v2[0]) != 0) return d<0?-1:1;
	if ((d=v1[1]-v2[1]) != 0) return d<0?-1:1;
	if ((d=v1[2]-v2[2]) != 0) return d<0?-1:1;
	return 0;
}

int compare(const Vector& v1, const Vector& v2, double tol)
{
	for (int i=0;i<3;i++) {
		double d=v1[i]-v2[i];
		if (d<-tol) return -1;
		if (d>tol) return 1;
	}
	return 0;
}

//*** Point

ostream& operator<<(ostream& s, const Point& p)
{
	return s << "Point(" << p[0] << "," << p[1] << "," << p[2] << ")\n";
}

Point operator*(const Point& p, const Frame& f)
{
	register double xi=p[0], yi=p[1], zi=p[2];
	return Point(xi*f[0][0]+yi*f[1][0]+zi*f[2][0]+f[3][0],
		     xi*f[0][1]+yi*f[1][1]+zi*f[2][1]+f[3][1],
		     xi*f[0][2]+yi*f[1][2]+zi*f[2][2]+f[3][2]);
}

Point& Point::operator*=(const Frame& f)
{
	return *this=*this*f;
}

int compare(const Point& p1, const Point& p2)
{
	double d;
	if ((d=p1[0]-p2[0]) != 0) return d<0?-1:1;
	if ((d=p1[1]-p2[1]) != 0) return d<0?-1:1;
	if ((d=p1[2]-p2[2]) != 0) return d<0?-1:1;
	return 0;
}

int compare(const Point& p1, const Point& p2, double tol)
{
	for (int i=0;i<3;i++) {
		double d=p1[i]-p2[i];
		if (d<-tol) return -1;
		if (d>tol) return 1;
	}
	return 0;
}

Vector cross(const Point& p1, const Point& p2, const Point& p3)
{
	return cross(p2-p1,p3-p1);
}

double area2(const Point& p1, const Point& p2, const Point& p3)
{
	return .25*mag2(cross(p1,p2,p3));
}

Point centroid(const Array<Point>& pa)
{
	assertx(pa.num());
	Homogeneous h;
	for (int i=0;i<pa.num();i++) h+=pa[i];
	return toPoint(h/pa.num());
}

//*** Frame

// destructor used to be non-existent
// introduced here for GNUG 2.5.8; inline {} would complain
//  "function may return with or without a value" when returning a Frame
//  with -O -Wall (not a problem without -O).
// problem linked to: subclass Vector v[3] + inline destructor {}
Frame::~Frame() { }

void Frame::zero()
{
	p=Point(0,0,0);
	v[0]=v[1]=v[2]=Vector(0,0,0);
}

void Frame::ident()
{
	p=Point(0,0,0);
	v[0]=Vector(1,0,0); v[1]=Vector(0,1,0); v[2]=Vector(0,0,1);
}

Frame::Frame(const Vector& v0, const Vector& v1,
		    const Vector& v2, const Point& pp)
: p(pp)
{
	v[0]=v0; v[1]=v1; v[2]=v2;
}

Frame Frame::translation(double x, double y, double z)
{
	return Frame(Vector(1,0,0),Vector(0,1,0),Vector(0,0,1),Point(x,y,z));
}

Frame Frame::translation(const Affine& vv)
{
	return translation(vv[0],vv[1],vv[2]);
}

Frame Frame::scaling(double x, double y, double z)
{
	return Frame(Vector(x,0,0),Vector(0,y,0),Vector(0,0,z),Point(0,0,0));
}

Frame Frame::scaling(double s)
{
	return scaling(s,s,s);
}

Frame Frame::identity() 
{
	return Frame(Vector(1,0,0),Vector(0,1,0),Vector(0,0,1),Point(0,0,0));
}

Frame operator*(const Frame& f1, const Frame& f2)
{
	return Frame(f1.v[0][0]*f2.v[0]+f1.v[0][1]*f2.v[1]+f1.v[0][2]*f2.v[2],
		     f1.v[1][0]*f2.v[0]+f1.v[1][1]*f2.v[1]+f1.v[1][2]*f2.v[2],
		     f1.v[2][0]*f2.v[0]+f1.v[2][1]*f2.v[1]+f1.v[2][2]*f2.v[2],
		     f1.p[0]*f2.v[0]+f1.p[1]*f2.v[1]+f1.p[2]*f2.v[2]+f2.p);
}

Frame& Frame::operator*=(const Frame& f)
{
	return *this=*this*f;
}

int Frame::invert()
{
	return ::invert(*this,*this);
}

Frame inverse(const Frame& f)
{
	Frame fr; assertx(invert(f,fr)); return fr;
}

Frame Frame::operator~() const
{
	return inverse(*this);
}

Frame& Frame::transpose()
{
	::transpose(*this,*this); return *this;
}

Frame transpose(const Frame& f)
{
	Frame fr; transpose(f,fr); return fr;
}

int Frame::isident() const
{
	if (compare(v[0],Vector(1,0,0))) return 0;
	if (compare(v[1],Vector(0,1,0))) return 0;
	if (compare(v[2],Vector(0,0,1))) return 0;
	return !compare(p,Point(0,0,0));
}

int invert(const Frame& fi, Frame& fo)
{
	// &fi==&fo is ok
	Matrix4 m=toMatrix4(fi);
	if (!invert(m,m)) return 0;
	fo=toFrame(m);
	return 1;
}

Frame Frame::rotation(int axis, double angle)
{
	assertx(axis>=0 && axis<=2);
	Frame f=Frame::identity();
	double c=cos(angle),s=sin(angle);
	switch (axis) {
	  case 0:
		f[0][0]=1; f[1][1]=c; f[1][2]=s; f[2][1]=-s; f[2][2]=c;
	  bcase 1:
		f[1][1]=1; f[2][2]=c; f[2][0]=s; f[0][2]=-s; f[0][0]=c;
	  bcase 2:
		f[2][2]=1; f[0][0]=c; f[0][1]=s; f[1][0]=-s; f[1][1]=c;
	}
	return f;
}

void transpose(const Frame& fi, Frame& fo)
{
	if (&fi!=&fo) fo=fi;
	assertx(!compare(fo.p,Point(0,0,0)));
	swap(&fo[1][0],&fo[0][1]);
	swap(&fo[2][0],&fo[0][2]);
	swap(&fo[2][1],&fo[1][2]);
}

void Frame::makerighthand()
{
	if (dot(cross(v[0],v[1]),v[2])<0) v[0]=-v[0];
}

ostream& operator<<(ostream& s, const Frame& f)
{
	return s << "Frame {\n  " << f.v[0] << "  " << f.v[1] << "  " <<
		f.v[2] << "  " << f.p << "}\n";
}
