// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Mk3d.h"

Mk3d::Mk3d(WA3dStream& pos)
: os(pos),
  ntrans(0), maxntrans(0), tottrans(0),
  ncolor(0), cc(A3dColor(1,.16,0),A3dColor(1,.5,.3),A3dColor(3,0,0)),
  maxvertices(0), totvertices(0),
  totpolygons(0), totpolylines(0), totpoints(0),
  forcepolyline(0), force2sided(0), forceflip(0)
{
	os.writeComment(" Init Mk3d");
	ctm.ident();
	ctmi.ident();
}

Mk3d::~Mk3d()
{
	os.writeComment(hform(" End Mk3d: %dgons %dlines %dpoints\
, %dverts(max %d), %d transf(max %d)",
			      totpolygons,totpolylines,totpoints,
			      totvertices,maxvertices,tottrans,maxntrans));
	os.flush();
	assertw(!ntrans);
	assertw(Sframe.empty()); // redundant
	assertw(Sframei.empty()); // redundant
	assertw(!ncolor);
	assertw(Sforcepolyline.empty());
	assertw(Sforce2sided.empty());
	assertw(Sforceflip.empty());
	assertw(!a3de.num());
}

WA3dStream& Mk3d::oa3d()
{
	return os;
}

void Mk3d::push()
{
	Sframe.push(new Frame(ctm));
	Sframei.push(new Frame(ctmi));
	ntrans++;
	tottrans++;
	if (ntrans>maxntrans) maxntrans=ntrans;
}

void Mk3d::pop()
{
	assertx(ntrans--);
	Frame* f=Sframe.pop();
	ctm=*f;
	delete f;
	Frame* fi=Sframei.pop();
	ctmi=*fi;
	delete fi;
}

void Mk3d::translate(const Vector& v)
{
	apply(Frame::translation(v));
}

void Mk3d::translate(double x, double y, double z)
{
	translate(Vector(x,y,z));
}

void Mk3d::rotate(int axis, double angle)
{
	apply(Frame::rotation(axis,angle));
}

void Mk3d::scale(double x, double y, double z)
{
	if (x<=0 && x!=-1 || y<=0 && y!=-1 || z<=0 && z!=-1)
		if (Warning("mk3d: strange scale()"))
			SHOWN(x),SHOWN(y),SHOWN(z);
	apply(Frame::scaling(x,y,z));
}

void Mk3d::scale(double v)
{
	scale(v,v,v);
}

void Mk3d::apply(const Frame& t)
{
	ctm=t*ctm;
	assertx(invert(ctm,ctmi));
}

Point Mk3d::transform(const Point& point)
{
	return point*ctm;
}

Vector Mk3d::transform(const Vector& normal)
{
	return ctmi*normal;
}

void Mk3d::pushcolor()
{
	Scolor.push(new A3dVertexColor(cc));
	ncolor++;
}

void Mk3d::popcolor()
{
	assertx(ncolor--);
	A3dVertexColor* c=Scolor.pop();
	cc=*c;
	delete c;
}

void Mk3d::diffuse(double r, double g, double b)
{
	cc.d=A3dColor(r,g,b);
}

void Mk3d::diffuse(const A3dColor& col)
{
	cc.d=col;
}

void Mk3d::specular(double r, double g, double b)
{
	cc.s=A3dColor(r,g,b);
}

void Mk3d::specular(const A3dColor& col)
{
	cc.s=col;
}

void Mk3d::phong(double p)
{
	cc.g=A3dColor(p,0,0);
}

void Mk3d::color(const A3dVertexColor& color)
{
	cc=color;
}

void Mk3d::scalecolor(double sr, double sg, double sb)
{
	cc.d.c[0]*=sr; cc.d.c[1]*=sg; cc.d.c[2]*=sb;
	cc.s.c[0]*=sr; cc.s.c[1]*=sg; cc.s.c[2]*=sb;
}

void Mk3d::point(const Point& p)
{
	a3de+=A3dVertex(transform(p),Vector(0,0,0),cc);
}

void Mk3d::point(double x, double y, double z)
{
	point(Point(x,y,z));
}

void Mk3d::normal(const Vector& normal)
{
	int mod=0;
	assertx(a3de.num());
	Vector tn=transform(normal);
	assertw(tn.normalize());
	for (int i=0;i<3;i++) if (fabs(tn[i])<1e-5) tn[i]=0,mod=1;
	if (mod) assertw(tn.normalize());
	a3de[a3de.num()-1].n=tn;
}

void Mk3d::normal(double x, double y, double z)
{
	normal(Vector(x,y,z));
}

void Mk3d::beginForcePolyline(int force)
{
	Sforcepolyline.push(forcepolyline);
	forcepolyline=force;
}

void Mk3d::endForcePolyline()
{
	assertx(!Sforcepolyline.empty());
	forcepolyline=Sforcepolyline.pop();
}

void Mk3d::beginForce2Sided(int force)
{
	Sforce2sided.push(force2sided);
	force2sided=force;
}

void Mk3d::endForce2Sided()
{
	assertx(!Sforce2sided.empty());
	force2sided=Sforce2sided.pop();
}

void Mk3d::beginForceFlip(int force)
{
	Sforceflip.push(forceflip);
	forceflip=force;
}

void Mk3d::endForceFlip()
{
	assertx(!Sforceflip.empty());
	forceflip=Sforceflip.pop();
}

void Mk3d::endpolygon()
{
	assertx(a3de.num()>=3);
	if (forceflip && !force2sided) flippoly();
	if (forcepolyline) {
		a3de+=a3de[0];
		endpolyline();
	} else if (force2sided) {
		end2polygon();
	} else {
		outputpoly();
		a3de.init(A3dElem::TPolygon); // clear it
	}
}

void Mk3d::end2polygon()
{
	if (forcepolyline) {
		endpolygon();
	} else {
		outputpoly();
		flippoly();
		outputpoly();
	}
	a3de.init(A3dElem::TPolygon);
}

void Mk3d::endpolyline()
{
	assertx(a3de.num()>=2);
	a3de.update(A3dElem::TPolyline);
	os.write(a3de);
	totpolylines++;
	totvertices+=a3de.num();
	if (a3de.num()>maxvertices) maxvertices=a3de.num();
	a3de.init(A3dElem::TPolygon);
}

void Mk3d::endpoint()
{
	assertx(a3de.num()==1);
	a3de.update(A3dElem::TPoint);
	os.write(a3de);
	totpoints++;
	totvertices+=1;
	if (a3de.num()>maxvertices) maxvertices=a3de.num();
	a3de.init(A3dElem::TPolygon);
}

void Mk3d::outputpoly()
{
	assertx(a3de.num()>=3);
	assertx(a3de.type()==A3dElem::TPolygon); // optional
	os.write(a3de);
	totpolygons++;
	totvertices+=a3de.num();
	if (a3de.num()>maxvertices) maxvertices=a3de.num();
}

void Mk3d::flippoly()
{
        int i;
	A3dVertex t;
	int n=a3de.num();
//	for (int i=1;i<=(n-1)/2;i++)
	for (i=1;i<=(n-1)/2;i++)
		t=a3de[i], a3de[i]=a3de[n-i], a3de[n-i]=t;
	for (i=0;i<n;i++)
		if (!a3de[i].n.iszero()) a3de[i].n=-a3de[i].n;
}
