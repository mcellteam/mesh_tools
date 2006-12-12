// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Mk3d_h
#define Mk3d_h

#include "A3dStream.h"
#include "Stack.h"

#define MKNEW(x) mk.push(); { x } mk.pop()

class Mk3d {
  public:
	Mk3d(WA3dStream& pos);
	~Mk3d();
	WA3dStream& oa3d();
	void push();
	void pop();
	void translate(double x, double y, double z);
	void translate(const Vector& v);
	void rotate(int axis, double angle);
	void scale(double x, double y, double z);
	void scale(double v);
	void apply(const Frame& t);
	Point transform(const Point& point);
	Vector transform(const Vector& normal);
	void pushcolor();
	void popcolor();
	void diffuse(double r, double g, double b);
	void diffuse(const A3dColor& col);
	void specular(double r, double g, double b);
	void specular(const A3dColor& col);
	void phong(double p);
	void color(const A3dVertexColor& color);
	void scalecolor(double sr, double sg, double sb);
	void point(double x, double y, double z);
	void point(const Point& p);
	void normal(double x, double y, double z);
	void normal(const Vector& normal);
	void beginForcePolyline(int force);
	void endForcePolyline();
	void beginForce2Sided(int force);
	void endForce2Sided();
	void beginForceFlip(int force);
	void endForceFlip();
	void endpolygon();
	void end2polygon();	// two sided polygon
	void endpolyline();
	void endpoint();
  private:
	WA3dStream& os;
	int ntrans;
	int maxntrans;
	int tottrans;
	Stack<Frame*> Sframe;
	Stack<Frame*> Sframei;
	Frame ctm;
	Frame ctmi;
	int ncolor;
	Stack<A3dVertexColor*> Scolor;
	A3dVertexColor cc;
	int maxvertices;
	int totvertices;
	A3dElem a3de;
	int totpolygons,totpolylines,totpoints;
	int forcepolyline,force2sided,forceflip;
	Stack<int> Sforcepolyline;
	Stack<int> Sforce2sided;
	Stack<int> Sforceflip;
	void outputpoly();
	void flippoly();
	DISABLECOPY(Mk3d);
};

#endif
