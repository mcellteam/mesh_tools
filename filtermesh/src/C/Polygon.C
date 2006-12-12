// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Polygon.h"
#include "A3dStream.h"
#include "Bbox.h"
#include "Array.h"

using std::ostream;

ALLOCATEPOOL(Polygon);

Polygon::Polygon(const A3dElem& el) : Array<Point>(0)
{
	assertx(el.type()==A3dElem::TPolygon);
	init(el.num());
	for (int i=0;i<num();i++) a[i]=el[i].p;
}

// used to be inline.  problem for GNUG 2.5.8; complained as for Frame.
Polygon::~Polygon() { }

void Polygon::copy(const Polygon& poly)
{
	init(poly.num());
	for (int i=0;i<num();i++) a[i]=poly.a[i];
}

void Polygon::getbbox(Bbox& bb) const
{
	assertx(num()>=1);
	bb.clear();
	for (int i=0;i<num();i++) bb.takeunion(a[i]);
}

Vector Polygon::getnormaldir() const
{
	if (num()==3) return cross(a[0],a[1],a[2]); // short-cut
	assertx(num()>=3);
	Vector nor(0,0,0);
	for (int i=1;i<num()-1;i++) nor+=cross(a[0],a[i],a[i+1]);
	return nor;
}

Vector Polygon::getnormal() const
{
	return oknormalize(getnormaldir());
}

double Polygon::getplanec(const Vector& pnor) const
{
	assertx(num()>=3);
	double sumd=0;
	for (int i=0;i<num();i++) sumd+=pvdot(a[i],pnor);
	return sumd/num();
}

double Polygon::gettolerance(const Vector& pnor, double d) const
{
	assertx(num()>=3);
	double tol=0;
	for (int i=0;i<num();i++) {
		double od=fabs(pvdot(a[i],pnor)-d);
		if (od>tol) tol=od;
	}
	return tol;
}

double Polygon::getarea() const
{
	assertx(num()>=3);
	double sum=0;
	for (int i=1;i<num()-1;i++) sum+=mysqrt(area2(a[0],a[i],a[i+1]));
	return sum;
}

int Polygon::intersectHyperplane(const Point& hp, const Vector& hn)
{
	assertx(num()>=3);
	static Array<double> sa; sa.init(num());
	int nin=0;
	for (int i=0;i<num();i++) {
		sa[i]=dot(a[i]-hp,hn);
		if (sa[i]>=0) nin++;
	}
	if (nin==num()) return 0;
	if (nin==0) { clear(); return 1; }
	static Polygon np; np.init(0);
	for (int vc=0;vc<num();vc++) {
		int vp=vc?vc-1:num()-1;
		int inc=sa[vc]>=0;
		int inp=sa[vp]>=0;
		if (inp+inc==1)
			np+=interp(a[vp],a[vc],sa[vc]/(sa[vc]-sa[vp]));
		if (inc)
			np+=a[vc];
	}
	this->copy(np);
	return 1;
}

int Polygon::intersectBbox(const Bbox& bb)
{
	assertx(num()>=3);
	int m=0;		// polygon_is_modified
	m|=intersectHyperplane(bb[0],Vector(+1,+0,+0)); if (!num()) return 1;
	m|=intersectHyperplane(bb[1],Vector(-1,+0,+0)); if (!num()) return 1;
	m|=intersectHyperplane(bb[0],Vector(+0,+1,+0)); if (!num()) return 1;
	m|=intersectHyperplane(bb[1],Vector(+0,-1,+0)); if (!num()) return 1;
	m|=intersectHyperplane(bb[0],Vector(+0,+0,+1)); if (!num()) return 1;
	m|=intersectHyperplane(bb[1],Vector(+0,+0,-1)); if (!num()) return 1;
	return m;
}

int Polygon::intersectSegment(const Point& p1, const Point& p2,
			      Point& pint) const
{
	assertx(num()>=3);
	Vector n=getnormal();
	assertx(!n.iszero());
	if (!IntersectPlaneSegment(n,getplanec(n),p1,p2,pint))
		return 0;
	return pointinside(n,pint);
}

static Vector gvint;

static int cmpinter(const Point* p1, const Point* p2)
{
	double a1=pvdot(*p1,gvint);
	double a2=pvdot(*p2,gvint);
	return a1<a2?-1:a1>a2?1:0;
}

void Polygon::intersectPlane(const Vector& polynor, const Vector& planenor,
			     double planed, double planetol,
			     Array<Point>& pa) const
{
        int i;
	assertx(num()>=3);
	static Array<double> sa; sa.init(num());
//	for (int i=0;i<num();i++) {
	for (i=0;i<num();i++) {
		double sc=pvdot(a[i],planenor)-planed;
		if (fabs(sc)<=planetol) sc=0;
		sa[i]=sc;
	}
	double sp=0;
	// make points lying in plane fall off to the side using propagation
	for (i=0;i<2*num();i++) {
		int i0=i%num();
		double sc=sa[i0];
		if (!sc && sp) sc=sa[i0]=1e-15*sign(sp);
		sp=sc;
	}
	pa.init(0);
	if (!sp) return;	// polygon lies in plane
	for (i=0;i<num();i++) {
		assertx(sa[i]);
		int i0=(i+0)%num();
		int i1=(i+1)%num();
		if (sa[i0]*sa[i1]>0) continue;
		pa+=interp(a[i0],a[i1],sa[i1]/(sa[i1]-sa[i0]));
	}
	assertx(pa.num()%2==0);
	if (!pa.num()) return;
	Vector vint=cross(polynor,planenor);
	assertx(vint.normalize()); // was 'if (!...) return'
	VectorStandardDirection(vint);
	gvint=vint;		// not reentrant
	// hqsort(pa,np,cmpinter);
	qsort(static_cast<Point*>(pa),pa.num(),sizeof(Point),
	      reinterpret_cast<int(*)(const void*,const void*)>(cmpinter));
}

inline double adjusttolerance(double tol)
{
	if (tol<1e-6) tol=1e-6;
	tol*=1.02;
	return tol;
}

void IntersectPolyPoly(const Polygon& p1, const Polygon& p2,
		       Array<Point>& pa)
{
	assertx(p1.num()>=3 && p2.num()>=3);
	Vector n1=p1.getnormal();
	Vector n2=p2.getnormal();
	double d1=p1.getplanec(n1);
	double d2=p2.getplanec(n2);
	double t1=adjusttolerance(p1.gettolerance(n1,d1));
	double t2=adjusttolerance(p2.gettolerance(n2,d2));
	Array<Point> pa1, pa2;
	p1.intersectPlane(n1, n2,d2,t2, pa1);
	p2.intersectPlane(n2, n1,d1,t1, pa2);
	int in1=0,in2=0,i1=0,i2=0,wasin=0;
	pa.init(0);
	Point* cp;
	for (;;) {
		if (i1==pa1.num() && i2==pa2.num()) break;
		// gvint is set in the last intersectPlane() call
		//  hence not reentrant
		if (i2==pa2.num() || (i1<pa1.num() &&
				      cmpinter(&pa1[i1],&pa2[i2])<=0)) {
			cp=&pa1[i1];
			in1=!in1;
			i1++;
		} else {
			cp=&pa2[i2];
			in2=!in2;
			i2++;
		}
		int in=in1&&in2;
		if (in!=wasin) {
			pa+=*cp;
			int pn=pa.num();
			if (!in && !compare(pa[pn-2],pa[pn-1],1e-6))
				pa.need(pn-2); // remove zero-length segment
		}
		wasin=in;
	}
	assertx(pa.num()%2==0 && !in1 && !in2);
}

int Polygon::pointinside(const Vector& pnor, const Point& point) const
{
	assertx(num()>=3);
	int axis=-1;
	double maxd=0;
	for (int c=0;c<3;c++) {
		double d=fabs(pnor[c]);
		if (d>maxd) maxd=d,axis=c;
	}
	assertx(maxd);
	int ax0=(axis+1)%3;
	int ax1=(axis+2)%3;
	double py=point[ax0];
	double pz=point[ax1];
	double y0=a[num()-1][ax0]-py,y1;
	double z0=a[num()-1][ax1]-pz,z1;
	int nint=0;
	for (int i=0;i<num();i++,y0=y1,z0=z1) {
		y1=a[i][ax0]-py;
		z1=a[i][ax1]-pz;
		if (z0>=0 && z1>=0) continue;
		if (z0<0 && z1<0) continue;
		if (y0<0 && y1<0) continue;
		if (y0>=0 && y1>=0) { nint++; continue; }
		if (y0-(y1-y0)/(z1-z0)*z0>=0) nint++;
	}
	return nint%2;
}

int Polygon::isconvex() const
{
	assertx(num()>=3);
	if (num()==3) return 1;
	Vector dir=getnormaldir();
	for (int i=0;i<num();i++) {
		Vector v=cross(a[i],a[(i+1)%num()],a[(i+2)%num()]);
		if (dot(v,dir)<0) return 0;
	}
	return 1;
}

ostream& operator<<(ostream& s, const Polygon& poly)
{
	s << "Polygon (" << poly.num() << " vertices) {\n";
	for (int i=0;i<poly.num();i++)
		s << "  " << poly.a[i];
	s << "}\n";
	return s;
}

int IntersectPlaneSegment(const Vector& normal, double d,
			  const Point& p1, const Point& p2,
			  Point& pint)
{
	double s1=pvdot(p1,normal)-d;
	double s2=pvdot(p2,normal)-d;
	if (s1<0 && s2<0 || s1>0 && s2>0) return 0;
	// what to do when segment lies in plane?  report nothing?
	if (!s1 && !s2) return 0;
	pint=interp(p1,p2,s2/(s2-s1));
	return 1;
}

Vector OrthogonalVector(const Vector& v)
{
	Vector vaxis(0,0,0);
	int minc=0;
	double mina=1e30;
	for (int c=0;c<3;c++) {
		double a=fabs(v[c]);
		if (a<mina) mina=a,minc=c;
	}
	vaxis[minc]=1;
	Vector vo=cross(v,vaxis);
	assertx(vo.normalize());
	return vo;
}

void VectorStandardDirection(Vector& v)
{
	int maxc=0;
	double maxa=fabs(v[0]);
	for (int c=1;c<3;c++) {
		double a=fabs(v[c]);
		if (a>maxa) maxa=a,maxc=c;
	}
	assertx(maxa);
	if (v[maxc]<0) v=-v;
}
