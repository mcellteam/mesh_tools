// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Mklib_h
#define Mklib_h

#include "Mk3d.h"
#include "Stack.h"

class Mklib {
  public:
	Mklib(Mk3d& mk3d);
	~Mklib();
	typedef void (*FUNC)(Mklib& mklib, int i);
	int smooth();
	void beginSmooth(int psmooth);
	void endSmooth();
	// transformation on object 
	// (-.5,-.5,-.5)<>(.5,.5,.5) with primary axis +x, secondary axis +y
	// -->  (-.5,-.5,0)<>(.5,.5,1) with primary axis +z, secondary axis +x
	void OtoU(FUNC func, int n);
	// unit square centered at origin, facing +x axis
	void squareO();
	// unit square between (0,0,0) and (1,1,0), facing +z
	void squareXY();
	// square above origin, in xy plane, facing +z axis
	void squareU();
	// unit cube centered at origin
	void cubeO();
	// unit cube between (0,0,0) and (1,1,1)
	void cubeXYZ();
	// cube with center of bottom face at origin
	void cubeU();
	// regular polygon, radius 1, normal to x axis, vertex on y axis
	void polygonO(int n);
	// polygon facing +z axis, vertex on +x axis
	void polygonU(int n);
	// radius 1 along +x axis, calls func with +x axis normal to circle
	// s=sin(PI/n)   h=cos(PI/n)
	// scaled to touch @(0,-s,0)&(0,+s,0), and center of circle @(-h,0,0)
	void circleOf(FUNC func, int n);
	void circleOfU(FUNC func, int n);
	// radius 1 along +z axis, calls func with +x axis normal to circle
	// not scaled -> center of circle @(-1,0,0)
	void radiusOfU(FUNC func, int n);
	// height h, 2 radii and normal angles w/respect XY plane
	void ringU(int n, double h, double r0, double r1, double a0, double a1);
	// ring with angles such that it is flat (not smooth)
	void flatringU(int n, double h, double r0, double r1);
	// circle with a hole in it
	void polyhole(int n, double r1);
	// a "discrete torus" with rectangular cross section
	void volumeringU(int n, double r1);
	// height 1, radius 1 open in +z axis, vertex on +x axis
	void tubeU(int n);
	// cylinder==tube with closed ends
	void cylinderU(int n);
	// height 1, radius 1, bottom at origin, peak at (1,0,0)
	void capU(int n);
	// cone==cap with closed bottom
	void coneU(int n);
	// radius 1, #latitudes(>1), #longitudes(>2)
	void sphere(int nlat, int nlong);
	// radius 1, #latitudes(>1), #longitudes(>2)
	void hemisphere(int nlat, int nlogn);
	// centered at centroid, edge=1 height=sqrt(2/3)
	void tetra();
	// bottom face centroid at origin, top at (0,0,sqrt(2/3))
	void tetraU();
  public:
	Mk3d& mk;
  private:
	int ismooth;
	Stack<int> Ssmooth;
	void gsphere(int nlat, int nlong, int hemi);
	DISABLECOPY(Mklib);
};


#endif
