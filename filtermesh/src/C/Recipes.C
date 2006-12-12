// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Recipes.h"		// for new declarations only

static int debug=GetenvValue("DEBUG");

void
LineOptimization(double (*func)(double), double t1, double t2,
		 double tol, double& tmin, double& emin, int& neval)
{
	double t3, e1, e2, e3;
	int neval1, neval2;
	Mnbrak(&t1,&t2,&t3,&e1,&e2,&e3,func,neval1);
	if (debug) {
		SHOWF("mnbrak t1=%12g e1=%g\n",t1,e1);
		SHOWF("mnbrak t2=%12g e2=%g\n",t2,e2);
		SHOWF("mnbrak t3=%12g e3=%g\n",t3,e3);
	}
	emin=Golden(t1,t2,t3,func,tol,&tmin,neval2);
	neval=neval1+neval2;
}

void
LineOptimization(double (*func)(double), double (*dfunc)(double),
		 double t1, double t2,
		 double tol, double& tmin, double& emin, int& neval)
{
	double t3, e1, e2, e3;
	int neval1, neval2;
	Mnbrak(&t1,&t2,&t3,&e1,&e2,&e3,func,neval1);
	if (debug) {
		SHOWF("mnbrak t1=%12g e1=%g\n",t1,e1);
		SHOWF("mnbrak t2=%12g e2=%g\n",t2,e2);
		SHOWF("mnbrak t3=%12g e3=%g\n",t3,e3);
	}
	// Dbrent will reevaluate func and dfunc at t2
	// it is unfortunate but does not cost all that much.
	// It is necessary since func may not be last evaluated at t2,
	// and dfunc may assumed that func has just been called there.
	emin=Dbrent(t1,t2,t3,func,dfunc,tol,&tmin,neval2);
	neval=neval1+neval2;
}


#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)

void
Mnbrak(double* ax, double* bx, double* cx,
       double* fa, double* fb, double* fc,
       double (*func)(double), int& neval)
{
	neval=0;
	double ulim,u,r,q,fu,dum;
	
	*fa=(*func)(*ax); neval++;
	*fb=(*func)(*bx); neval++;
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum);
		SHFT(dum,*fb,*fa,dum);
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx); neval++;
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u); neval++;
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u); neval++;
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u); neval++;
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx));
				SHFT(*fb,*fc,fu,(*func)(u)); neval++;
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u); neval++;
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u); neval++;
		}
		SHFT(*ax,*bx,*cx,u);
		SHFT(*fa,*fb,*fc,fu);
	}
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SIGN
#undef SHFT 

#define R 0.61803399
#define C (1.0-R)
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double
Golden(double ax, double bx, double cx, double (*func)(double),
       double tol, double* xmin, int& neval)
{
	neval=0;
	double f0,f1,f2,f3,x0,x1,x2,x3;

	x0=ax;
	x3=cx;
	if (fabs(cx-bx) > fabs(bx-ax)) {
		x1=bx;
		x2=bx+C*(cx-bx);
	} else {
		x2=bx;
		x1=bx-C*(bx-ax);
	}
	f1=(*func)(x1); neval++;
	f2=(*func)(x2); neval++;
	while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
		if (f2 < f1) {
			SHFT(x0,x1,x2,R*x1+C*x3);
			SHFT(f0,f1,f2,(*func)(x2)); neval++;
		} else {
			SHFT(x3,x2,x1,R*x2+C*x0);
			SHFT(f3,f2,f1,(*func)(x1)); neval++;
		}
	}
	if (f1 < f2) {
		*xmin=x1;
		return f1;
	} else {
		*xmin=x2;
		return f2;
	}
}

#undef C
#undef R 


#define ITMAX 100
#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

double
Dbrent(double ax, double bx, double cx,
       double (*func)(double), double (*dfunc)(double),
       double tol, double* xmin, int& neval)
{
	neval=0;
	int iter,ok1,ok2;
	double a,b,d=99e30,d1,d2,du,dv,dw,dx,e=0.0;
	double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*func)(x); neval++;
	dw=dv=dx=(*dfunc)(x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol1=tol*fabs(x)+ZEPS;
		tol2=2.0*tol1;
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			d1=2.0*(b-a);
			d2=d1;
			if (dw != dx)  d1=(w-x)*dx/(dx-dw);
			if (dv != dx)  d2=(v-x)*dx/(dx-dv);
			u1=x+d1;
			u2=x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			olde=e;
			e=d;
			if (ok1 || ok2) {
				if (ok1 && ok2)
					d=(fabs(d1) < fabs(d2) ? d1 : d2);
				else if (ok1)
					d=d1;
				else
					d=d2;
				if (fabs(d) <= fabs(0.5*olde)) {
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				} else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
				}
			} else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		} else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
		if (fabs(d) >= tol1) {
			u=x+d;
			fu=(*func)(u); neval++;
		} else {
			u=x+SIGN(tol1,d);
			fu=(*func)(u); neval++;
			if (fu > fx) {
				*xmin=x;
				return fx;
			}
		}
		du=(*dfunc)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			MOV3(v,fv,dv, w,fw,dw)
			MOV3(w,fw,dw, x,fx,dx)
			MOV3(x,fx,dx, u,fu,du)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				MOV3(v,fv,dv, w,fw,dw)
				MOV3(w,fw,dw, u,fu,du)
			} else if (fu < fv || v == x || v == w) {
				MOV3(v,fv,dv, u,fu,du)
			}
		}
	}
	nrerror("Too many iterations in routine DBRENT");
	return 0;
}

#undef ITMAX
#undef ZEPS
#undef SIGN
#undef MOV3 
