// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Principal.h"
#include "Homogeneous.h"
#include "Recipes.h"		// jacobi(), eigsrt()

void PrincipalComponents(const Point pa[], int n, Frame& f, double eimag[3])
{
        int i;
	double** mcovar=matrix(1,3,1,3);
	double* d=vector(1,3);
	double** v=matrix(1,3,1,3);
	assertx(n>0);
	Homogeneous hp;
//	for (int i=0;i<n;i++) hp+=pa[i];
	for (i=0;i<n;i++) hp+=pa[i];
	Point avgp=toPoint(hp/n);
	for (int c0=0;c0<3;c0++)
		for (int c1=0;c1<=c0;c1++) {
			double a=0;
			for (i=0;i<n;i++)
				a+=(pa[i][c0]-avgp[c0])*(pa[i][c1]-avgp[c1]);
			a/=n;
			mcovar[1+c0][1+c1]=a;
			if (c1<c0) mcovar[1+c1][1+c0]=a;
		}
	int nrot;
	jacobi(mcovar,3,d,v,&nrot);
	eigsrt(d,v,3);
	for (i=0;i<3;i++) {
		double a=d[1+i];
		if (a<0) a=0;	// for numerics
		a=mysqrt(a);
		eimag[i]=a;
		if (!a) a=1e-15; // very small but non-zero vector
		f.v[i]=Vector(v[1][1+i],v[2][1+i],v[3][1+i])*a;
	}
	f.p=avgp;
	f.makerighthand();
	free_matrix(mcovar,1,3,1,3);
	free_vector(d,1,3);
	free_matrix(v,1,3,1,3);
}
