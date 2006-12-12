// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

// 06/16/92	extracted from Geometry.C
// 02/13/93	moved all inlines to end

#ifndef Geometry_h
#define Geometry_h

#include "Array.h"
#include "Pool.h"
#include <iostream>

class Frame;
class Point;

class Affine {
  public:
	double& operator[](int i) { return c[i]; }
	double operator[](int i) const { return c[i]; }
	int iszero() const;

	operator double const *() const { return c; }

	friend std::ostream& operator<<(std::ostream& s, const Affine& a);
	POOLALLOCATION(Affine);
#if __GNUG__==2 && __GNUC_MINOR__==5
	// causes problems without in 2.5.8: delete[]
	// with in 2.4.5: not implemented: aggregate value in COND_EXPR
	~Affine() {}
#endif
  protected:
	double c[3];
	Affine() {};
	Affine(double x, double y, double z);
};

class Vector : public Affine {
  public:
	Vector() {}
#if __GNUG__==2 && __GNUC_MINOR__==5
	~Vector() {}
#endif
	Vector(double x, double y, double z);
	friend inline Vector toVector(const Point& p);
	Vector operator-() const;
	friend inline Vector operator+(const Vector& v1, const Vector& v2);
	friend inline Vector operator-(const Vector& v1, const Vector& v2);
	Vector& operator+=(const Vector& v);
	Vector& operator-=(const Vector& v);
	friend inline Vector operator*(const Vector& v, double f);
	friend inline Vector operator*(double f, const Vector& v);
	friend inline Vector operator/(const Vector& v, double f);
	Vector& operator*=(double f);
	Vector& operator/=(double f);
	friend Vector operator*(const Vector& v, const Frame& f);
	Vector& operator*=(const Frame& f);
	friend Vector operator*(const Frame& f, const Vector& normal);
	friend inline Vector cross(const Vector& v1, const Vector& v2);
	friend int normalize(const Vector& vi, Vector& vo);
	int normalize();
	friend Vector oknormalize(const Vector& v);
	friend inline double mag2(const Vector& v);
	friend inline double mag(const Vector& v);
	friend inline double dot(const Vector& v1, const Vector& v2);
	friend int compare(const Vector& v1, const Vector& v2);
	friend int compare(const Vector& v1, const Vector& v2, double tol);
	friend std::ostream& operator<<(std::ostream& s, const Vector& v);
	// for stupid cxx
	double& operator[](int i) { return Affine::operator[](i); }
	double operator[](int i) const { return Affine::operator[](i); }
};

class Point : public Affine {
  public:
	Point() {}
#if __GNUG__==2 && __GNUC_MINOR__==5
	~Point() {}
#endif
	Point(double x, double y, double z);
	friend inline Point toPoint(const Vector& v);
	friend inline Vector operator-(const Point& p1, const Point& p2);
	friend inline Point operator+(const Point& p, const Vector& v);
	friend inline Point operator+(const Vector& v, const Point& p);
	friend inline Point operator-(const Point& p, const Vector& v);
	Point& operator+=(const Vector& v);
	Point& operator-=(const Vector& v);
	friend Point operator*(const Point& p, const Frame& f);
	Point& operator*=(const Frame& f);
	friend inline double pvdot(const Point& p, const Vector& v);
	friend inline double dist2(const Point& p1, const Point& p2);
	friend inline double dist(const Point& p1, const Point& p2);
	friend int compare(const Point& p1, const Point& p2);
	friend int compare(const Point& p1, const Point& p2, double tol);
	friend inline Point interp(const Point& p1, const Point& p2,
				   double f1=.5);
	friend inline Point interp(const Point& p1, const Point& p2,
				   const Point& p3, double f1=1./3.,
				   double f2=1./3.);
	friend Vector cross(const Point& p1, const Point& p2, const Point& p3);
	friend double area2(const Point& p1, const Point& p2, const Point& p3);
	friend std::ostream& operator<<(std::ostream& s, const Point& p);
	// for stupid cxx
	double& operator[](int i) { return Affine::operator[](i); }
	double operator[](int i) const { return Affine::operator[](i); }
};

class Frame {
  public:
	Frame() : v(), p() {}
	Frame(const Vector& v0, const Vector& v1,
	      const Vector& v2, const Point& pp);
	~Frame();
	Vector v[3];
	Point p;		// comes immediately after v for operator[]!
	
	void zero();
	void ident();
	Affine& operator[](int i);
	const Affine& operator[](int i) const;
	int isident() const;
	friend Frame operator*(const Frame& f1, const Frame& f2);
	Frame& operator*=(const Frame& f);
	friend int invert(const Frame& fi, Frame& fo);
	int invert();
	friend Frame inverse(const Frame& f);
	Frame operator~() const;
	friend void transpose(const Frame& fi, Frame& fo);
	Frame& transpose();
	friend Frame transpose(const Frame& f);
	void makerighthand();	// flip an axis if necessary
	static Frame translation(const Affine& vv);
	static Frame translation(double x, double y, double z);
	static Frame rotation(int axis, double angle);
	static Frame scaling(double x, double y, double z);
	static Frame scaling(double s);
	static Frame identity();
	friend std::ostream& operator<<(std::ostream& s, const Frame& f);
};

class Bary {			// Barycentric coordinates
  public:
	Bary() { }
	Bary(double x, double y, double z);
	double& operator[](int i) { return p[i]; }
	double operator[](int i) const { return p[i]; }
  private:
	double p[3];
};


extern Point centroid(const Array<Point>& pa);


//----------------------------------------------------------------------------

//*** Affine

inline Affine& Frame::operator[](int i) { return v[i]; }

inline const Affine& Frame::operator[](int i) const { return v[i]; }

inline Affine::Affine(double x, double y, double z) { c[0]=x; c[1]=y; c[2]=z; }

inline int Affine::iszero() const
{
	return !c[0] && !c[1] && !c[2];
}

//*** Vector

inline Vector::Vector(double x, double y, double z) : Affine(x,y,z) { }

inline Vector toVector(const Point& p)
{
	return Vector(p[0],p[1],p[2]);
}

inline Vector Vector::operator-() const
{
	return Vector(-c[0],-c[1],-c[2]);
}

inline Vector operator+(const Vector& v1, const Vector& v2)
{
	return Vector(v1.c[0]+v2.c[0],v1.c[1]+v2.c[1],v1.c[2]+v2.c[2]);
}

inline Vector operator-(const Vector& v1, const Vector& v2)
{
	return Vector(v1.c[0]-v2.c[0],v1.c[1]-v2.c[1],v1.c[2]-v2.c[2]);
}

inline Vector& Vector::operator+=(const Vector& v)
{
	c[0]+=v.c[0]; c[1]+=v.c[1]; c[2]+=v.c[2]; return *this;
}

inline Vector& Vector::operator-=(const Vector& v)
{
	c[0]-=v.c[0]; c[1]-=v.c[1]; c[2]-=v.c[2]; return *this;
}

inline Vector operator*(const Vector& v, double f)
{
	return Vector(v.c[0]*f,v.c[1]*f,v.c[2]*f);
}

inline Vector operator*(double f, const Vector& v)
{
	return Vector(v.c[0]*f,v.c[1]*f,v.c[2]*f);
}

inline Vector operator/(const Vector& v, double f)
{
	return Vector(v.c[0]/f,v.c[1]/f,v.c[2]/f);
}

inline Vector& Vector::operator*=(double f)
{
	c[0]*=f; c[1]*=f; c[2]*=f; return *this;
}

inline Vector& Vector::operator/=(double f)
{
	c[0]/=f; c[1]/=f; c[2]/=f; return *this;
}

inline Vector cross(const Vector& v1, const Vector& v2)
{
	register double v1x=v1[0], v1y=v1[1], v1z=v1[2];
	register double v2x=v2[0], v2y=v2[1], v2z=v2[2];
	return Vector(v1y*v2z-v1z*v2y, v1z*v2x-v1x*v2z, v1x*v2y-v1y*v2x);
}

inline double mag2(const Vector& v)
{
	return v.c[0]*v.c[0]+v.c[1]*v.c[1]+v.c[2]*v.c[2];
}

inline double mag(const Vector& v)
{
	return mysqrt(mag2(v));
}

inline double dot(const Vector& v1, const Vector& v2)
{
	return v1.c[0]*v2.c[0]+v1.c[1]*v2.c[1]+v1.c[2]*v2.c[2];
}

//*** Point

inline Point::Point(double x, double y, double z) : Affine(x,y,z) { }

inline Point toPoint(const Vector& v) { return Point(v[0],v[1],v[2]); }

inline Vector operator-(const Point& p1, const Point& p2)
{
	return Vector(p1.c[0]-p2.c[0],p1.c[1]-p2.c[1],p1.c[2]-p2.c[2]);
}

inline Point operator+(const Point& p, const Vector& v)
{
	return Point(p.c[0]+v[0],p.c[1]+v[1],p.c[2]+v[2]);
}

inline Point operator+(const Vector& v, const Point& p)
{
	return Point(v[0]+p.c[0],v[1]+p.c[1],v[2]+p.c[2]);
}

inline Point operator-(const Point& p, const Vector& v)
{
	return Point(p.c[0]-v[0],p.c[1]-v[1],p.c[2]-v[2]);
}

inline Point& Point::operator+=(const Vector& v)
{
	c[0]+=v[0]; c[1]+=v[1]; c[2]+=v[2]; return *this;
}

inline Point& Point::operator-=(const Vector& v)
{
	c[0]-=v[0]; c[1]-=v[1]; c[2]-=v[2]; return *this;
}

inline double pvdot(const Point& p, const Vector& v)
{
	return p.c[0]*v[0]+p.c[1]*v[1]+p.c[2]*v[2];
}

inline double dist2(const Point& p1, const Point& p2)
{
	register double x=p1.c[0]-p2.c[0],y=p1.c[1]-p2.c[1],z=p1.c[2]-p2.c[2];
	return x*x+y*y+z*z;
}

inline double dist(const Point& p1, const Point& p2)
{
	return mysqrt(dist2(p1,p2));
}

inline Point interp(const Point& p1, const Point& p2, double f1)
{
	double f2=1-f1;
	return Point(f1*p1.c[0]+f2*p2.c[0],
		     f1*p1.c[1]+f2*p2.c[1],
		     f1*p1.c[2]+f2*p2.c[2]);
}

inline Point interp(const Point& p1, const Point& p2, const Point& p3,
		    double f1, double f2) {
	double f3=1-f1-f2;
	return Point(f1*p1.c[0]+f2*p2.c[0]+f3*p3.c[0],
		     f1*p1.c[1]+f2*p2.c[1]+f3*p3.c[1],
		     f1*p1.c[2]+f2*p2.c[2]+f3*p3.c[2]);
}

//*** Bary

inline Bary::Bary(double x, double y, double z) { p[0]=x; p[1]=y; p[2]=z; }

#endif
