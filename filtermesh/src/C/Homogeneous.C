// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Homogeneous.h"

using std::ostream;

//*** Homogeneous

ostream& operator<<(ostream& s, const Homogeneous& h)
{
	return s << "Homogeneous(" << h.c[0] << "," << h.c[1] <<
		"," << h.c[2] << "," << h.c[3] << ")\n";
}

//*** HFrame

ostream& operator<<(ostream& s, const HFrame& f)
{
	s << "HFrame {\n";
	for (int i=0;i<4;i++)
		s << "  " << f[i];
	s << "}\n";
	return s;
}

Matrix4 toMatrix4(const Frame& f)
{
	Matrix4 m;
	m[0][3]=m[1][3]=m[2][3]=0; m[3][3]=1;
	m[0][0]=f.v[0][0]; m[0][1]=f.v[0][1]; m[0][2]=f.v[0][2];
	m[1][0]=f.v[1][0]; m[1][1]=f.v[1][1]; m[1][2]=f.v[1][2];
	m[2][0]=f.v[2][0]; m[2][1]=f.v[2][1]; m[2][2]=f.v[2][2];
	m[3][0]=f.p[0]; m[3][1]=f.p[1]; m[3][2]=f.p[2];
	return m;
}

Frame toFrame(const Matrix4& m)
{
	Frame f;
	const double TOLERANCE=1e-5;
	if (fabs(m[0][3])>TOLERANCE ||
	    fabs(m[1][3])>TOLERANCE ||
	    fabs(m[2][3])>TOLERANCE ||
	    fabs(m[3][3]-1)>TOLERANCE)
		if (Warning("Frame matrix strange")) SHOW(m);
	f.v[0][0]=m[0][0]; f.v[0][1]=m[0][1]; f.v[0][2]=m[0][2];
	f.v[1][0]=m[1][0]; f.v[1][1]=m[1][1]; f.v[1][2]=m[1][2];
	f.v[2][0]=m[2][0]; f.v[2][1]=m[2][1]; f.v[2][2]=m[2][2];
	f.p[0]=m[3][0]; f.p[1]=m[3][1]; f.p[2]=m[3][2];
	return f;
}
