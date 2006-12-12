// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Recipes_h
#define Recipes_h

#ifndef ANSI
#define ANSI
#define ANSIWASNOTDEFINED
#endif

extern "C" {
#include "nr.h"
#include "nrutil.h"
}

#ifdef ANSIWASNOTDEFINED
#undef ANSIWASNOTDEFINED
#undef ANSI
#endif

extern void
LineOptimization(double (*func)(double), double t1, double t2,
		 double tol, double& tmin, double& emin, int& neval);
void
LineOptimization(double (*func)(double), double (*dfunc)(double),
		 double t1, double t2,
		 double tol, double& tmin, double& emin, int& neval);

extern void
Mnbrak(double* ax, double* bx, double* cx,
       double* fa, double* fb, double* fc,
       double (*func)(double), int& neval);

extern double
Golden(double ax, double bx, double cx, double (*f)(double),
       double tol, double* xmin, int& neval);

extern
double Dbrent(double ax, double bx, double cx,
	     double (*func)(double), double (*dfunc)(double),
	     double tol, double* xmin, int& neval);

#endif
