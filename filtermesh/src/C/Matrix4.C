// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Matrix4.h"

using std::ostream;

const double PERSPTHRESH=1e-3;

//*** Vector4

Vector4 operator*(const Vector4& v, const Matrix4& m)
{
	Vector4 vr;
	for (int i=0;i<4;i++) {
		double a=0;
		for (int j=0;j<4;j++)
			a+=v[j]*m[j][i];
		vr[i]=a;
	}
	return vr;
}

Vector4 operator*(const Matrix4& m, const Vector4& v)
{
	Vector4 vr;
	for (int i=0;i<4;i++) {
		double a=0;
		for (int j=0;j<4;j++)
			a+=m[i][j]*v[j];
		vr[i]=a;
	}
	return vr;
}

Vector4& Vector4::operator*=(const Matrix4& m)
{
	return *this=*this*m;
}

//*** Matrix4

void Matrix4::zero()
{
	for (int i=0;i<4;i++) for (int j=0;j<4;j++) v[i][j]=0;
}

void Matrix4::ident()
{
	zero();
	v[0][0]=v[1][1]=v[2][2]=v[3][3]=1;
}

Matrix4 operator*(const Matrix4& m1, const Matrix4& m2)
{
	Matrix4 mr;
	for (int i=0;i<4;i++)
		for (int j=0;j<4;j++) {
			double a=0;
			for (int k=0;k<4;k++)
				a+=m1[i][k]*m2[k][j];
			mr[i][j]=a;
		}
	return mr;
}

Matrix4& Matrix4::operator*=(const Matrix4& m)
{
	return *this=*this*m;
}

int invert(const Matrix4& mi, Matrix4& mo)
{
	double t[4][8];
	int i,j;
	for (i=0;i<4;i++)
		for (j=0;j<4;j++) {
			t[i][j]=mi[i][j];
			t[i][j+4]=0;
		}
	for (i=0;i<4;i++)
		t[i][i+4]=1;
	for (i=0;i<4;i++) {
		if (i<3) {	// swap row with largest front coefficient
			double a=fabs(t[i][i]),ab;
			int m=i;
			for (int l=i+1;l<4;l++)
				if ((ab=fabs(t[l][i]))>a)
					a=ab,m=l;
			if (m!=i)
				for (j=0;j<8;j++)
					swap(&t[i][j],&t[m][j]);
		}
		if (!t[i][i])
			return 0;
		for (j=0;j<4;j++) {
			if (i==j) {
				double a=1/t[i][i];
				for (int k=0;k<8;k++)
					t[j][k]*=a;
			} else {
				double a=-t[j][i]/t[i][i];
				for (int k=0;k<8;k++)
					t[j][k]+=a*t[i][k];
			}
		}
	}
	for (i=0;i<4;i++)
		for (j=0;j<4;j++)
			mo[i][j]=t[i][j+4];
	return 1;
}

Matrix4 inverse(const Matrix4& m)
{
	Matrix4 mr; assertx(invert(m,mr)); return mr;
}

Matrix4 Matrix4::operator~() const
{
	return inverse(*this);
}

Matrix4& Matrix4::transpose()
{
	swap(&v[0][1],&v[1][0]); swap(&v[0][2],&v[2][0]);
	swap(&v[0][3],&v[3][0]); swap(&v[1][2],&v[2][1]);
	swap(&v[1][3],&v[3][1]); swap(&v[2][3],&v[3][2]);
	return *this;
}

Matrix4 transpose(const Matrix4& m)
{
	Matrix4 mr=m; return mr.transpose();
}

Matrix4& Matrix4::fix()
{
	if (fabs(v[0][3])+fabs(v[1][3])+fabs(v[2][3])+
	    fabs(v[3][3]-1)>PERSPTHRESH) return *this;
	v[0][3]=v[1][3]=v[2][3]=0;
	v[3][3]=1;
	return *this;
}

ostream& operator<<(ostream& s, const Vector4& v)
{
	return s << "Vector4(" << v[0] << "," << v[1] <<
		"," << v[2] << "," << v[3] << ")\n";
}

ostream& operator<<(ostream& s, const Matrix4& m)
{
	s << "Matrix4 {\n";
	for (int i=0;i<4;i++)
		s << "  " << m[i];
	s << "}\n";
	return s;
}
