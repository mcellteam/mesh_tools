// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Mklib.h"

const double p90=torad(90);
const double n90=torad(-90);
const double p18=torad(180);

static double ft0,ft1,ft2,ft3,ft4,ft5;

static void wedge0(Mklib& mklib, int);
static void poly0(Mklib& mklib, int);

Mklib::Mklib(Mk3d& mk3d) : mk(mk3d), ismooth(1) { }

Mklib::~Mklib() { }

int Mklib::smooth() { return ismooth; }

void Mklib::beginSmooth(int psmooth)
{
	Ssmooth.push(ismooth);
	ismooth=psmooth;
}

void Mklib::endSmooth()
{
	assertx(!Ssmooth.empty());
	ismooth=Ssmooth.pop();
}

void Mklib::OtoU(FUNC func, int n)
{
	mk.push(); {
		mk.rotate(1,p90);
		mk.rotate(0,n90);
		mk.translate(.5,0,0);
		func(*this,n);
	} mk.pop();
}

void Mklib::squareO()
{
	mk.point(0,-.5,-.5);
	mk.point(0,+.5,-.5);
	mk.point(0,+.5,+.5);
	mk.point(0,-.5,+.5);
	mk.endpolygon();
}

void Mklib::squareXY()
{
	MKNEW( mk.translate(.5,.5,0); mk.rotate(1,p90); squareO(); );
}

void Mklib::squareU()
{
	MKNEW( mk.rotate(1,n90); squareO(); );
}

void Mklib::cubeO()
{
	MKNEW( mk.translate(+.0,+.0,-.5); mk.rotate(1,p90); squareO(); );
	MKNEW( mk.translate(+.0,+.0,+.5); mk.rotate(1,n90); squareO(); );
	MKNEW( mk.translate(-.5,+.0,+.0); mk.rotate(2,p18); squareO(); );
	MKNEW( mk.translate(+.5,+.0,+.0); mk.rotate(2,0  ); squareO(); );
	MKNEW( mk.translate(+.0,-.5,+.0); mk.rotate(2,n90); squareO(); );
	MKNEW( mk.translate(+.0,+.5,+.0); mk.rotate(2,p90); squareO(); );
}

void Mklib::cubeXYZ()
{
	MKNEW( mk.translate(.5,.5,.5); cubeO(); );
}

void Mklib::cubeU()
{
	MKNEW( mk.translate(0,0,.5); cubeO(); );
}

void Mklib::polygonO(int n)
{
	for (int i=0;i<n;i++) {
		double a=2*PI*i/n;
		mk.point(0,cos(a),sin(a));
	}
	mk.endpolygon();
}

void Mklib::polygonU(int n)
{
	MKNEW( mk.rotate(1,n90); mk.rotate(0,n90); polygonO(n); );
}

void Mklib::circleOf(FUNC func, int n)
{
	double a=2*PI/n,h=cos(a*.5);
	mk.push(); {
		mk.rotate(1,p90);
		mk.rotate(2,p90);
		mk.rotate(2,a*.5);
		for (int i=0;i<n;i++) {
			MKNEW( mk.translate(h,0,0); func(*this,i); );
			mk.rotate(2,a);
		}
	} mk.pop();
}

void Mklib::circleOfU(FUNC func, int n)
{
	double a=2*PI/n,h=cos(a*.5);
	mk.push(); {
		mk.rotate(2,a*.5);
		for (int i=0;i<n;i++) {
			MKNEW( mk.translate(h,0,0); func(*this,i); );
			mk.rotate(2,a);
		}
	} mk.pop();
}

void Mklib::radiusOfU(FUNC func, int n)
{
	double a=2*PI/n;
	mk.push(); {
		for (int i=0;i<n;i++) {
			MKNEW( mk.translate(1,0,0); func(*this,i); );
			mk.rotate(2,a);
		}
	} mk.pop();
}

void Mklib::ringU(int n, double h, double r0, double r1, double a0, double a1)
{
	assertx((r0>0 || r1>0) && r0>=0 && r1>=0);
	mk.push(); {
		if (r0<=0) {
			mk.translate(0,0,h);
			mk.rotate(0,p18);
			ringU(n,h,r1,r0,-a1,-a0);
		} else {
			ft0=r1/r0; ft3=h/r0;
			ft4=sin(PI/n); ft5=cos(PI/n);
			ft1=tan(a0); ft2=tan(a1);
			mk.scale(r0);
			circleOfU(wedge0,n);
		}
	} mk.pop();
}

static void wedge0(Mklib& mklib, int)
{
	mklib.mk.point(0,ft4,0);
	if (mklib.smooth()) mklib.mk.normal(ft5,ft4,ft1);
	mklib.mk.point((ft0-1)*ft5,ft4*ft0,ft3);
	if (ft0>1e-6) {
		if (mklib.smooth()) mklib.mk.normal(ft5,ft4,ft2);
		mklib.mk.point((ft0-1)*ft5,-ft4*ft0,ft3);
		if (mklib.smooth()) mklib.mk.normal(ft5,-ft4,ft2);
	} else {
		if (mklib.smooth()) mklib.mk.normal(1,0,ft2);
	}
	mklib.mk.point(0,-ft4,0);
	if (mklib.smooth()) mklib.mk.normal(ft5,-ft4,ft1);
	mklib.mk.endpolygon();
}

void Mklib::flatringU(int n, double h, double r0, double r1)
{
	assertx(r0*h || r0*r0-r1);
	double a0=atan2(r0*h,r0*r0-r1);
	ringU(n,h,r0,r1,a0,a0);
}

void Mklib::polyhole(int n, double r1)
{
	ft1=r1; ft2=PI/n;
	circleOfU(poly0,n);
}

static void poly0(Mklib& mklib, int)
{
	double s=sin(ft2),h=cos(ft2);
	mklib.mk.point(0,s,0);
	mklib.mk.point(-h*(1-ft1),+s*ft1,0);
	mklib.mk.point(-h*(1-ft1),-s*ft1,0);
	mklib.mk.point(0,-s,0);
	mklib.mk.endpolygon();
}

void Mklib::volumeringU(int n, double r1)
{
	MKNEW( mk.rotate(0,p18); polyhole(n,r1); );
	MKNEW( mk.translate(0,0,1); polyhole(n,r1); );
	MKNEW( tubeU(n); );
	MKNEW( mk.scale(r1,r1,1);
	     mk.beginForceFlip(1); tubeU(n); mk.endForceFlip(); );
}

void Mklib::tubeU(int n)
{
	ringU(n,1,1,1,0,0);
}

void Mklib::cylinderU(int n)
{
	tubeU(n);
	MKNEW( mk.rotate(0,p18); polygonU(n); );
	MKNEW( mk.translate(0,0,1); polygonU(n); );
}

void Mklib::capU(int n)
{
	ringU(n,1,1,0,PI/4,PI/4);
}

void Mklib::coneU(int n)
{
	capU(n);
	MKNEW( mk.rotate(0,p18); polygonU(n); );
}

void Mklib::sphere(int nlat, int nlong)
{
	gsphere(nlat,nlong,0);
}

void Mklib::hemisphere(int nlat, int nlong)
{
	gsphere(nlat*2,nlong,1);
}

void Mklib::gsphere(int nlat, int nlong, int hemi)
{
	assertx(nlat>1 && nlong>2);
	for (int i=0;i<nlat;i++) {
		if (hemi && i<nlat<2) continue;
		double a1=(-.5+double(i)/nlat)*PI;
		double a2=(-.5+double(i+1)/nlat)*PI;
		double c1=fabs(cos(a1)),c2=fabs(cos(a2));
		double s1=sin(a1),s2=sin(a2);
		if (i==0) a1=-.499999*PI,c1=0;
		else if (i==nlat-1) a2=.499999*PI,c2=0;
		MKNEW( mk.translate(0,0,s1); ringU(nlong,s2-s1,c1,c2,a1,a2); );
	}
}

void Mklib::tetra()
{
	double xp=1/mysqrt(3),xn=xp/-2,yp=.5,yn=-.5,zp=1.5/mysqrt(6),zn=zp/-3;
	mk.point(xp,0,zn);
	mk.point(xn,yn,zn);
	mk.point(xn,yp,zn);
	mk.endpolygon();
	mk.point(xp,0,zn);
	mk.point(0,0,zp);
	mk.point(xn,yn,zn);
	mk.endpolygon();
	mk.point(xn,yn,zn);
	mk.point(0,0,zp);
	mk.point(xn,yp,zn);
	mk.endpolygon();
	mk.point(xn,yp,zn);
	mk.point(0,0,zp);
	mk.point(xp,0,zn);
	mk.endpolygon();
}

void Mklib::tetraU()
{
	MKNEW( mk.translate(.5/mysqrt(6),0,0); tetra(); );
}
