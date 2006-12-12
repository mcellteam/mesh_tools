// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef A3dStream_h
#define A3dStream_h

#include "Geometry.h"
#include "Array.h"
#include <iostream>

class A3dColor {
  public:
	A3dColor() { }
	A3dColor(double r, double g, double b) { c[0]=r; c[1]=g; c[2]=b; }
	double c[3];		// r,g,b
	friend int compare(const A3dColor& c1, const A3dColor& c2);
	friend std::ostream& operator<<(std::ostream& s, const A3dColor& c);
};

class A3dVertexColor {
  public:
	A3dVertexColor() : d(), s(), g() { }
	A3dVertexColor(const A3dColor& pd); // specular=white, phong=1
	A3dVertexColor(const A3dColor& pd, const A3dColor& ps,
		       const A3dColor& pg);
	A3dColor d;
	A3dColor s;
	A3dColor g;
	friend int compare(const A3dVertexColor& c1, const A3dVertexColor& c2);
	static const A3dVertexColor White,Black,Red,Green,Blue;
};

class A3dVertex {
  public:
	A3dVertex() : p(), n(), c() { }
	A3dVertex(const Point& pp, const Vector& pn, const A3dVertexColor& pc);
	Point p;
	Vector n;
	A3dVertexColor c;
};

// I thought about making this a base class for different types of elements,
// but thought it would be too inefficient to allocate/deallocate every time.
class A3dElem {
  public:
	enum Type { TPolygon='P', TPolyline='L', TPoint='p',
		    TComment='#',
		    TEndObject='o', TEndFrame='f', TEndFile='q',
		    TEditObject='O' };
	// also reserved: 'v','E','n','d','s','g'
	A3dElem();
	A3dElem(Type type, int binary=0, int nv=0); // allocates AND init()
	~A3dElem();
	// both A3dElem(..nv) and init() allocate and initialize for nv
	void init(Type type, int binary=0, int nv=0);
	void update(Type type, int binary=0); // polygon<>polyline<>point
	void copy(const A3dElem& e); // deep copy, shallow copy not allowed
	Type type() const;
	void setbinary(int b);
	int binary() const;
	static int statustype(Type type);
	static int commandtype(Type type);
// for TPolygon || TPolyline || TPoint:
	int num() const;
	void addvertices(int i);
	A3dElem& operator+=(const A3dVertex& vertex);
	A3dVertex& operator[](int i);
	const A3dVertex& operator[](int i) const;
	Vector pnormal() const;	// may be degenerate (zero)!
// for TComment:
	void setcomment(const char* str); // str may start with ' '
	const char* comment() const; // user's responsibility to dup string
// for commandtype():
	double& f(int i);
	double f(int i) const;
  private:
	Type itype;
	int ibinary;
	Array<A3dVertex> v;
	union {
		struct {	// TComment
			const char* s;
		} comment;
		struct {	// commandtype()
			double f[3];
		} other;
	} u;
	void realloc(int nv);
	DISABLECOPY(A3dElem);
};

class A3dStream {
  protected:
	A3dStream();
	virtual ~A3dStream();
	A3dVertexColor curcol;
	void setcurcolor(int type, const double f[3]);
	enum { LINELENGTH=1023 }; // max length of ascii line
	enum { BINARYCODE=3 };
  private:
	DISABLECOPY(A3dStream);
};


class RA3dStream : public A3dStream {
  public:
	RA3dStream();
	virtual ~RA3dStream();
	void read(A3dElem& e);
  protected:
	virtual int readline(int& binary, int& type, double f[3],
			     const char*& comment)=0; // ret success
};

class RSA3dStream : public RA3dStream {	// Read from stream
  public:
	RSA3dStream(std::istream& pis);
	~RSA3dStream();
	std::istream& is();
  protected:
	int readline(int& binary, int& type, double f[3], const char*& comment);
  private:
	std::istream& iis;
};


class WA3dStream : public A3dStream {
  public:
	WA3dStream();
	virtual ~WA3dStream();
	void write(const A3dElem& e);
	void writeComment(const char* string); // can contain newlines
	void writeEndObject(int binary=0, double f0=1, double f1=1);
	void writeClearObject(int binary=0, double f0=1, double f1=0);
	void writeEndFrame(int binary=0);
	virtual void flush()=0;
  protected:
	int oldformat;		// old a3d format
	virtual void output(int binary, int type,
			    double f1, double f2, double f3)=0;
	virtual void outputcomment(const char* s)=0;
	virtual void blankline()=0;
  private:
	int first;		// first write?
	int forcechoicebinary;	// force choice one way or the other
	int choicebinary;	// which way is forced
	int pblank;		// previous element left a blank line
	void writeoldformat(const A3dElem& e);
};

class WSA3dStream : public WA3dStream {	// Write to stream
  public:
	WSA3dStream(std::ostream& pos);
	~WSA3dStream();
	void flush();
	std::ostream& os();
  protected:
	void output(int binary, int type, double f1, double f2, double f3);
	void outputcomment(const char* s);
	void blankline();
  private:
	std::ostream& ios;
};

//----------------------------------------------------------------------------

//*** A3dVertexColor

inline A3dVertexColor::A3dVertexColor(const A3dColor& pd)
: d(pd), s(1,1,1), g(1,0,0) { }

inline A3dVertexColor::A3dVertexColor(const A3dColor& pd, const A3dColor& ps,
				      const A3dColor& pg)
: d(pd), s(ps), g(pg) { }

//*** A3dVertex

inline A3dVertex::A3dVertex(const Point& pp, const Vector& pn,
			    const A3dVertexColor& pc)
: p(pp), n(pn), c(pc) { }

//*** A3dElem

inline A3dElem::A3dElem() : itype(TPolygon), ibinary(0) { }

inline int A3dElem::statustype(Type type)
{
	return type=='d' || type=='s' || type=='g';
}

inline int A3dElem::commandtype(Type type)
{
	return type==A3dElem::TEndObject || type==A3dElem::TEndFrame ||
		type==A3dElem::TEndFile || type==A3dElem::TEditObject;
}

inline A3dElem::Type A3dElem::type() const { return itype; }

inline int A3dElem::num() const { return v.num(); }

inline A3dElem& A3dElem::operator+=(const A3dVertex& vertex)
{ v+=vertex; return *this; }

inline A3dVertex& A3dElem::operator[](int i) { return v[i]; }

inline const A3dVertex& A3dElem::operator[](int i) const { return v[i]; }

#endif
