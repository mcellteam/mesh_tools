// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Homogeneous_h
#define Homogeneous_h

#include "Geometry.h"
#include "Matrix4.h"

#include <iostream>

class HFrame;

class Homogeneous : public Vector4 {
  public:
	Homogeneous();
	Homogeneous(double x, double y, double z, double w);
	Homogeneous(const Vector4& v);
	Homogeneous(const Point& p);
	Homogeneous(const Vector& p);
	int iszero() const;
	friend inline Homogeneous normalize(const Homogeneous& h);
	Homogeneous& normalize();
	friend inline Homogeneous operator*(const Point& p, double f);
	friend inline Homogeneous operator*(double f, const Point& p);
	// next 3 are intermediates to allow Point/Vector->Homogeneous conv.
	friend inline Homogeneous operator*(const Homogeneous& h,
					    const HFrame& f);
	friend inline Homogeneous operator+(const Homogeneous& h1,
					    const Homogeneous& h2);
	friend inline Homogeneous operator-(const Homogeneous& h1,
					    const Homogeneous& h2);
	Homogeneous& operator+=(const Homogeneous& h);
	Homogeneous& operator-=(const Homogeneous& h);
	friend std::ostream& operator<<(std::ostream& s, const Homogeneous& h);
	// for stupid cxx
	double& operator[](int i) { return Vector4::operator[](i); }
	double operator[](int i) const { return Vector4::operator[](i); }
};

class HFrame : public Matrix4 {
  public:
	HFrame() { }
	HFrame(const Homogeneous& h0, const Homogeneous& h1,
	       const Homogeneous& h2, const Homogeneous& h3);
	HFrame(const Matrix4& m);
	HFrame(const Frame& f);
	friend std::ostream& operator<<(std::ostream& s, const HFrame& f);
	// for stupid cxx
	Vector4& operator[](int i) { return Matrix4::operator[](i); }
	const Vector4& operator[](int i) const {return Matrix4::operator[](i);}
};

inline Point toPoint(const Vector4& h);
inline Vector toVector(const Vector4& h);
extern Matrix4 toMatrix4(const Frame& f);
extern Frame toFrame(const Matrix4& m);

//----------------------------------------------------------------------------

//*** Homogeneous

inline Homogeneous::Homogeneous() : Vector4(0,0,0,0) { }

inline Homogeneous::Homogeneous(double x, double y, double z, double w)
: Vector4(x,y,z,w) { }

inline Homogeneous::Homogeneous(const Vector4& v) : Vector4(v) { }

inline Homogeneous::Homogeneous(const Point& p) : Vector4(p[0],p[1],p[2],1) { }

inline Homogeneous::Homogeneous(const Vector& v): Vector4(v[0],v[1],v[2],0) { }

inline int Homogeneous::iszero() const
{
	return !c[0] && !c[1] && !c[2] && !c[3];
}

inline Homogeneous normalize(const Homogeneous& h)
{
	register double a=1/assertv(h[3]);
	return Homogeneous(h[0]*a,h[1]*a,h[2]*a,h[3]*a);
}

inline Homogeneous& Homogeneous::normalize()
{
	return *this=::normalize(*this);
}

inline Homogeneous operator*(const Point& p, double f)
{
	return Homogeneous(p[0]*f,p[1]*f,p[2]*f,f);
}

inline Homogeneous operator*(double f, const Point& p)
{
	return p*f;
}

inline Homogeneous operator*(const Homogeneous& h, const HFrame& f)
{
	return Homogeneous(reinterpret_cast<Vector4 const &>(h) * reinterpret_cast<Matrix4 const &>(f));
}

inline Homogeneous operator+(const Homogeneous& h1, const Homogeneous& h2)
{
	return Homogeneous(reinterpret_cast<Vector4 const &>(h1) + reinterpret_cast<Vector4 const &>(h2));
}

inline Homogeneous operator-(const Homogeneous& h1, const Homogeneous& h2)
{
	return Homogeneous(reinterpret_cast<Vector4 const &>(h1) - reinterpret_cast<Vector4 const &>(h2));
}

inline Homogeneous& Homogeneous::operator+=(const Homogeneous& h)
{
	return *this = *this + h;
}

inline Homogeneous& Homogeneous::operator-=(const Homogeneous& h)
{
	return *this = *this - h;
}

//*** HFrame

inline HFrame::HFrame(const Homogeneous& h0, const Homogeneous& h1,
		      const Homogeneous& h2, const Homogeneous& h3)
: Matrix4(h0,h1,h2,h3) { }

inline HFrame::HFrame(const Matrix4& m) : Matrix4(m) { }

inline HFrame::HFrame(const Frame& f) : Matrix4(toMatrix4(f)) { }

//*** other

inline Point toPoint(const Vector4& h)
{
	const double TOLERANCE=1e-5;
	assertx(fabs(h[3]-1)<TOLERANCE);
	return Point(h[0],h[1],h[2]);
}

inline Vector toVector(const Vector4& h)
{
	const double TOLERANCE=1e-5;
	assertx(fabs(h[3])<TOLERANCE);
	return Vector(h[0],h[1],h[2]);
}

#endif
