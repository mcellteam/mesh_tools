// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Facedistance.h"


double DistPointTriangle2(const Point& p, const Point& p1,
			 const Point& p2, const Point& p3)
{
	double dis2;
	ProjectPTri(p,p1,p2,p3,&dis2,0,0,0);
	return dis2;
}

static void projecth(const Point& p, const Point& p1,
		     const Point& p2, const Point& p3, const Bary& ba,
		     double* pdis2, Bary* pcba, Point* pclp)
{
	// projection lies outside triangle, so more work is needed.
	Point proj=interp(p1,p2,p3,ba[0],ba[1]);
	Point pf[3]; pf[0]=p1,pf[1]=p2,pf[2]=p3;
	Bary cba;
	double mind2=1e30;
	for (int i=0;i<3;i++) {
		if (ba[(i+2)%3]>=0) continue;
		// project proj onto segment pf[(i+0)%3]--pf[(i+1)%3]
		Vector vvi=pf[(i+1)%3]-pf[i];
		Vector vppi=proj-pf[i];
		double d12sq=mag2(vvi);
		double don12=dot(vvi,vppi);
		if (don12<=0) {
			double d2=dist2(pf[i],proj);
			if (d2>=mind2) continue;
			mind2=d2; cba[i]=1; cba[(i+1)%3]=0; cba[(i+2)%3]=0;
		} else if (don12>=d12sq) {
			double d2=dist2(pf[(i+1)%3],proj);
			if (d2>=mind2) continue;
			mind2=d2; cba[i]=0; cba[(i+1)%3]=1; cba[(i+2)%3]=0;
		} else {
			double a=don12/d12sq;
			cba[i]=1-a; cba[(i+1)%3]=a; cba[(i+2)%3]=0;
			break;
		}
	}
	if (pcba) (*pcba)[0]=cba[0],(*pcba)[1]=cba[1],(*pcba)[2]=cba[2];
	if (pclp || pdis2) {
		Point clp=interp(p1,p2,p3,cba[0],cba[1]);
		if (pclp) *pclp=clp;
		if (pdis2) *pdis2=dist2(p,clp);
	}
}

// two bad cases:
// - v2==0 or v3==0 (two points of triangle are same) -> ok
// - v23*v23==v22*v33 (area zero but neither v2 nor v3 zero) -> ouch
void ProjectPTri(const Point& p, const Point& p1,
		 const Point& p2, const Point& p3,
		 double* pdis2, Bary* pba, Bary* pcba, Point* pclp)
{
	// Vector v2=p2-p1, v3=p3-p1, vpp1=p-p1;
	// register double v22=mag2(v2), v33=mag2(v3), v23=dot(v2,v3);
	// register double v2pp1=dot(v2,vpp1), v3pp1=dot(v3,vpp1);
	register double p1x=p1[0], p1y=p1[1], p1z=p1[2];
	register double v2x=p2[0]-p1x, v2y=p2[1]-p1y, v2z=p2[2]-p1z;
	register double v3x=p3[0]-p1x, v3y=p3[1]-p1y, v3z=p3[2]-p1z;
	register double vpp1x=p[0]-p1x, vpp1y=p[1]-p1y, vpp1z=p[2]-p1z;
	register double v22=v2x*v2x+v2y*v2y+v2z*v2z;
	register double v33=v3x*v3x+v3y*v3y+v3z*v3z;
	register double v23=v2x*v3x+v2y*v3y+v2z*v3z;
	register double v2pp1=v2x*vpp1x+v2y*vpp1y+v2z*vpp1z;
	register double v3pp1=v3x*vpp1x+v3y*vpp1y+v3z*vpp1z;
	if (!v22) v22=1;	// recover if v2==0
	if (!v33) v33=1;	// recover if v3==0
	register double a2,a3;
	register double denom=(v33-v23*v23/v22);
	if (!denom) {
		a2=a3=1./3.;	// recover if v23*v23==v22*v33
	} else {
		a3=(v3pp1-v23/v22*v2pp1)/denom;
		a2=(v2pp1-a3*v23)/v22;
	}
	register double a1=1-a2-a3;
	if (pba) (*pba)[0]=a1,(*pba)[1]=a2,(*pba)[2]=a3;
	if (a1<0 || a2<0 || a3<0) {
		Bary ba; ba[0]=a1,ba[1]=a2,ba[2]=a3;
		projecth(p,p1,p2,p3,ba,pdis2,pcba,pclp);
		return;
	}
	// fast common case
	if (pcba) (*pcba)[0]=a1,(*pcba)[1]=a2,(*pcba)[2]=a3;
	if (!pclp && !pdis2)
		return;
	// Point clp=interp(p1,p2,p3,a1,a2);
	register double clpx=p1x+a2*v2x+a3*v3x;
	register double clpy=p1y+a2*v2y+a3*v3y;
	register double clpz=p1z+a2*v2z+a3*v3z;
	if (pclp) {
		double* pclpf=&(*pclp)[0];
		pclpf[0]=clpx; pclpf[1]=clpy; pclpf[2]=clpz;
	}
	if (pdis2) {
		// *pdis2=dist2(p,clp);
		register double x=p[0]-clpx, y=p[1]-clpy, z=p[2]-clpz;
		*pdis2=x*x+y*y+z*z;
	}
}
