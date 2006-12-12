// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Matrix4_h
#define Matrix4_h

#include <iostream>

class Matrix4;

class Vector4 {
  public:
	Vector4() { }
	Vector4(double x, double y, double z, double w);
	double& operator[](int i) { return c[i]; }
	double operator[](int i) const { return c[i]; }
	friend inline Vector4 operator+(const Vector4& v1, const Vector4& v2);
	friend inline Vector4 operator-(const Vector4& v1, const Vector4& v2);
	Vector4& operator+=(const Vector4& v);
	Vector4& operator-=(const Vector4& v);
	friend inline Vector4 operator*(const Vector4& v, double f);
	friend inline Vector4 operator*(double f, const Vector4& v);
	Vector4& operator*=(double f);
	friend inline Vector4 operator/(const Vector4& v, double f);
	Vector4& operator/=(double f);
	friend Vector4 operator*(const Vector4& v, const Matrix4& m);
	friend Vector4 operator*(const Matrix4& m, const Vector4& v);
	Vector4& operator*=(const Matrix4& m);
	friend std::ostream& operator<<(std::ostream& s, const Vector4& v);

	operator double const *() const { return c; }
  protected:
	double c[4];
};

class Matrix4 {
  public:
	Matrix4() { }
	Matrix4(const Vector4& v0, const Vector4& v1,
		const Vector4& v2, const Vector4& v3);
	Vector4& operator[](int i) { return v[i]; }
	const Vector4& operator[](int i) const { return v[i]; }
	void zero();
	void ident();
	friend Matrix4 operator*(const Matrix4& m1, const Matrix4& m2);
	Matrix4& operator*=(const Matrix4& m);
	friend int invert(const Matrix4& mi, Matrix4& mo);
	friend Matrix4 inverse(const Matrix4& m);
	Matrix4 operator~() const;
	Matrix4& transpose();
	friend Matrix4 transpose(const Matrix4& m);
	friend std::ostream& operator<<(std::ostream& s, const Matrix4& m);
	Matrix4& fix();		// make affine if close to affine
  protected:
	Vector4 v[4];
};

//----------------------------------------------------------------------------

//*** Vector4

inline Vector4::Vector4(double x, double y, double z, double w)
{
	c[0]=x; c[1]=y; c[2]=z; c[3]=w;
}

inline Vector4 operator+(const Vector4& v1, const Vector4& v2)
{
	return Vector4(v1.c[0]+v2.c[0],v1.c[1]+v2.c[1],
		       v1.c[2]+v2.c[2],v1.c[3]+v2.c[3]);
}

inline Vector4 operator-(const Vector4& v1, const Vector4& v2)
{
	return Vector4(v1.c[0]-v2.c[0],v1.c[1]-v2.c[1],
		       v1.c[2]-v2.c[2],v1.c[3]-v2.c[3]);
}

inline Vector4& Vector4::operator+=(const Vector4& v)
{
	c[0]+=v.c[0]; c[1]+=v.c[1]; c[2]+=v.c[2]; c[3]+=v.c[3]; return *this;
}

inline Vector4& Vector4::operator-=(const Vector4& v)
{
	c[0]-=v.c[0]; c[1]-=v.c[1]; c[2]-=v.c[2]; c[3]-=v.c[3]; return *this;
}

inline Vector4 operator*(const Vector4& v, double f)
{
	return Vector4(v.c[0]*f,v.c[1]*f,v.c[2]*f,v.c[3]*f);
}

inline Vector4 operator*(double f, const Vector4& v)
{
	return Vector4(v.c[0]*f,v.c[1]*f,v.c[2]*f,v.c[3]*f);
}

inline Vector4& Vector4::operator*=(double f)
{
	c[0]*=f; c[1]*=f; c[2]*=f; c[3]*=f; return *this;
}

inline Vector4 operator/(const Vector4& v, double f)
{
	return v*(1/f);
}

inline Vector4& Vector4::operator/=(double f)
{
	return *this*=(1/f);
}

//*** Matrix4

inline Matrix4::Matrix4(const Vector4& v0, const Vector4& v1,
			const Vector4& v2, const Vector4& v3)
{
	v[0]=v0; v[1]=v1; v[2]=v2; v[3]=v3;
}

#endif
