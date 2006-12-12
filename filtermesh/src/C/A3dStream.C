// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "A3dStream.h"

using std::istream;
using std::ostream;
using std::ws;

const A3dColor ColorZERO(0,0,0);
const A3dColor ColorUNDEF(-1,0,0);


//*** A3dColor

int compare(const A3dColor& c1, const A3dColor& c2)
{
	return c1.c[0]!=c2.c[0] || c1.c[1]!=c2.c[1] || c1.c[2]!=c2.c[2];
}

ostream& operator<<(ostream& s, const A3dColor& c)
{
	return s << "Color(r=" << c.c[0] << ",g=" << c.c[1] <<
		",b=" << c.c[2] << ")\n";
}

//*** A3dVertexColor

// added =A3dVertexColor to work around bug in GNU, it still complains though
const A3dVertexColor A3dVertexColor::White=A3dVertexColor(A3dColor(1,1,1));
const A3dVertexColor A3dVertexColor::Black=A3dVertexColor(A3dColor(0,0,0));
const A3dVertexColor A3dVertexColor::Red=A3dVertexColor(A3dColor(1,0,0));
const A3dVertexColor A3dVertexColor::Green=A3dVertexColor(A3dColor(0,1,0));
const A3dVertexColor A3dVertexColor::Blue=A3dVertexColor(A3dColor(0,0,1));

int compare(const A3dVertexColor& c1, const A3dVertexColor& c2)
{
	return compare(c1.d,c2.d) || compare(c1.s,c2.s) || compare(c1.g,c2.g);
}

//*** A3dElem

A3dElem::~A3dElem() {}		// GNUG 2.5.8

void A3dElem::addvertices(int i)
{
	assertx(itype==TPolygon || itype==TPolyline || itype==TPoint);
	if (itype==TPoint) assertx(v.num()+i==1);
	v.need(v.num()+i);
}

void A3dElem::init(Type type, int binary, int nv)
{
	itype=type;
	ibinary=binary;
	switch (itype) {
	  case TPolygon:
	  ocase TPolyline:
	  ocase TPoint:
		v.init(nv);
	  bcase TComment:
		u.comment.s=0;
	  bdefault:
		u.other.f[0]=u.other.f[1]=u.other.f[2]=0;
	}
}

A3dElem::A3dElem(Type type, int binary, int nv)
{
	init(type,binary,nv);
}

void A3dElem::update(Type type, int binary)
{
	assertx(itype==TPolygon || itype==TPolyline || itype==TPoint);
	itype=type;
	assertx(itype==TPolygon || itype==TPolyline || itype==TPoint);
	ibinary=binary;
	if (itype==TPoint) assertx(v.num()<=1);
}

void A3dElem::setcomment(const char* str)
{
	assertx(itype==TComment);
	u.comment.s=str;
}

const char* A3dElem::comment() const
{
	assertx(itype==TComment);
	return u.comment.s;
}

void A3dElem::setbinary(int b)
{
	ibinary=b;
}

int A3dElem::binary() const
{
	return ibinary;
}

double& A3dElem::f(int i)
{
	assertx(commandtype(itype));
	assertx(i>=0 && i<3);
	return u.other.f[i];
}

double A3dElem::f(int i) const
{
	assertx(commandtype(itype));
	assertx(i>=0 && i<3);
	return u.other.f[i];
}

void A3dElem::copy(const A3dElem& e)
{
	assertx(this!=&e);
	init(e.itype,e.ibinary);
	if (itype==TPolygon || itype==TPolyline || itype==TPoint) {
		v.init(e.num());
		for (int i=0;i<e.num();i++) v[i]=e.v[i];
	} else if (itype==TComment) {
		u.comment.s=e.u.comment.s;
	} else {
		u.other=e.u.other;
	}
}

Vector A3dElem::pnormal() const
{
	Vector vt(0,0,0);
	assertx(itype==TPolygon && v.num()>=3);
	for (int i=1;i<v.num()-1;i++)
		vt+=cross(v[0].p,v[i].p,v[i+1].p);
	(void)vt.normalize();
	return vt;
}

//*** A3dStream

A3dStream::A3dStream() { }

A3dStream::~A3dStream() { }

void A3dStream::setcurcolor(int type, const double f[3])
{
	switch(type) {
	  case 'd':  curcol.d=A3dColor(f[0],f[1],f[2]);
	  bcase 's': curcol.s=A3dColor(f[0],f[1],f[2]);
	  bcase 'g': curcol.g=A3dColor(f[0],f[1],f[2]);
	  bdefault: assertnever("");
	    }
}

//*** RA3dStream

RA3dStream::RA3dStream()
{
	curcol.d=curcol.s=curcol.g=ColorZERO;
}

RA3dStream::~RA3dStream() { }

void RA3dStream::read(A3dElem& e)
{
	int i,binary,type;
	Vector normal(0,0,0);
	double f[3];
	const char* commentstr;
	for (;;) {
		i=readline(binary,type,f,commentstr);
		if (!i) binary=0,type=A3dElem::TEndFile,f[0]=f[1]=f[2]=0;
		if (type=='n') {
			assertx(normal.iszero());
			normal=Vector(f[0],f[1],f[2]);
			continue;
		} else if (A3dElem::statustype(A3dElem::Type(type))) {
			setcurcolor(type,f);
			continue;
		}
		e.init(A3dElem::Type(type),binary);
		if (type==A3dElem::TComment) {
			assertx(normal.iszero());
			e.setcomment(commentstr);
			return;
		} else if (A3dElem::commandtype(A3dElem::Type(type))) {
			assertx(normal.iszero());
			e.f(0)=f[0]; e.f(1)=f[1]; e.f(2)=f[2];
			return;
		} else if (type==A3dElem::TPoint) {
			e+=A3dVertex(Point(f[0],f[1],f[2]),normal,curcol);
			normal=Vector(0,0,0);
			return;
		} else if (type==A3dElem::TPolygon ||
			   type==A3dElem::TPolyline) {
			break;
		} else {
			SHOWF("RA3dStream: unknown type '%c'\n",char(type));
			assertnever("type unknown");
		}
	}
	assertx(normal.iszero());
	int etype,nv=0;
	for (;;) {
		i=readline(binary,etype,f,commentstr);
		if (!i) assertnever("RA3dStream: EOF within poly");
		if (etype=='E') {
			break;
		} else if (etype=='n') {
			assertx(normal.iszero());
			normal=Vector(f[0],f[1],f[2]);
		} else if (etype=='v') {
			e+=A3dVertex(Point(f[0],f[1],f[2]),normal,curcol);
			normal=Vector(0,0,0);
			nv++;
		} else if (A3dElem::statustype(A3dElem::Type(etype))) {
			setcurcolor(etype,f);
		} else {
			SHOWF("RA3dStream, within polygon, type '%c'\n",
			      char(etype));
			assertnever("that type does not belong there");
		}
	}
	assertx(normal.iszero());
	switch(type) {
	  case A3dElem::TPolygon:   assertx(nv>=3);
	  bcase A3dElem::TPolyline: assertx(nv>=2);
	  bdefault: assertnever("?internal");
	    }
}

//*** RSA3dStream

RSA3dStream::RSA3dStream(istream& pis) : iis(pis) { }

RSA3dStream::~RSA3dStream() { }

int RSA3dStream::readline(int& binary, int& type, double f[3],
			  const char*& comment)
{
	iis >> ws;
	char c=iis.peek();
	if ((binary=c==A3dStream::BINARYCODE) != 0) { // binary
		char str[16];
		unsigned short utype;
		if (!iis.read(str,16)) return 0;
		StdToShort(&str[2],&utype),type=utype;
		for (int i=0;i<3;i++) StdToFloat(&str[i*4+4],&f[i]);
		return 1;
	}
	iis.ignore();
	type=c;
	if (type==A3dElem::TComment) {
		static char str[A3dStream::LINELENGTH];	// not reentrant
		iis.get(str,sizeof(str));
		comment=str;
	} else {
		if (iis) assertw1(iis.peek()==' ');
		iis >> f[0] >> f[1] >> f[2];
	}
	int isgood=iis?1:0;
	iis.ignore(9999,'\n');
	return isgood;
}

istream& RSA3dStream::is()
{
	return iis;
}

//*** WA3dStream

WA3dStream::WA3dStream() : first(1), pblank(0) { }

WA3dStream::~WA3dStream() { }

void WA3dStream::write(const A3dElem& e)
{
	int binary=e.binary();
	A3dElem::Type type=e.type();
	if (first) {
		first=0;
		curcol.d=curcol.s=curcol.g=ColorUNDEF;
		if ((forcechoicebinary=getenv("A3DBINARY")?1:0) != 0)
			choicebinary=GetenvValue("A3DBINARY");
		oldformat=GetenvValue("A3DOLD");
		writeComment(hform(" Created by WA3dStream on %s",CTime()));
	}
	if (forcechoicebinary) binary=choicebinary;
	if (type==A3dElem::TPolygon || type==A3dElem::TPolyline ||
	    type==A3dElem::TPoint) {
		// common case, skip for now
	} else if (A3dElem::commandtype(type)) {
		if (!binary && !pblank) blankline();
		output(binary,type,e.f(0),e.f(1),e.f(2));
		if (!binary) blankline(),pblank=1;
		curcol.d=curcol.s=curcol.g=ColorUNDEF;
		flush();
		return;
	} else if (type==A3dElem::TComment) {
		outputcomment(e.comment()),pblank=0;
		return;
	} else {
		assertnever("Unrecognized type");
	}
	if (type==A3dElem::TPolygon) assertx(e.num()>=3);
	else if (type==A3dElem::TPolyline) assertx(e.num()>=2);
	else if (type==A3dElem::TPoint) assertx(e.num()==1);
	if (oldformat) { writeoldformat(e); return; }
	if (type==A3dElem::TPolygon || type==A3dElem::TPolyline) {
		if (!binary && !pblank) blankline();
		output(binary,type,0,0,0);
	}
	for (int i=0;i<e.num();i++) {
		if (compare(curcol.d,e[i].c.d)) {
			curcol.d=e[i].c.d;
			output(binary,'d',curcol.d.c[0],curcol.d.c[1],
			       curcol.d.c[2]);
		}
		if (compare(curcol.s,e[i].c.s)) {
			curcol.s=e[i].c.s;
			output(binary,'s',curcol.s.c[0],curcol.s.c[1],
			       curcol.s.c[2]);
		}
		if (compare(curcol.g,e[i].c.g)) {
			curcol.g=e[i].c.g;
			output(binary,'g',curcol.g.c[0],curcol.g.c[1],
			       curcol.g.c[2]);
		}
		if (!e[i].n.iszero()) {
			Vector nl=e[i].n;
			if (mag2(nl)>.99)
				for (int c=0;c<3;c++)
					if (nl[c] && fabs(nl[c])<1e-5) nl[c]=0;
			output(binary,'n',nl[0],nl[1],nl[2]);
		}
		output(binary,type==A3dElem::TPoint?'p':'v',
		       e[i].p[0],e[i].p[1],e[i].p[2]);
	}
	pblank=0;
	if (type==A3dElem::TPolygon || type==A3dElem::TPolyline) {
		output(binary,'E',0,0,0);
		if (!binary) blankline(),pblank=1;
	}
}

void WA3dStream::writeoldformat(const A3dElem& e)
{
	int binary=e.binary();
	A3dElem::Type type=e.type();
	for (int i=0;i<e.num();i++) {
		if (compare(curcol.d,e[i].c.d)) {
			curcol.d=e[i].c.d;
			output(binary,'d',curcol.d.c[0],curcol.d.c[1],
			       curcol.d.c[2]);
		}
		if (compare(curcol.s,e[i].c.s)) {
			curcol.s=e[i].c.s;
			output(binary,'s',curcol.s.c[0],curcol.s.c[1],
			       curcol.s.c[2]);
		}
		if (compare(curcol.g,e[i].c.g)) {
			curcol.g=e[i].c.g;
			output(binary,'p',curcol.g.c[0],curcol.g.c[1],
			       curcol.g.c[2]);
		}
		output(binary,type==A3dElem::TPoint?'P':
		       type==A3dElem::TPolygon?(i?'l':'m'):(i?'1':'0'),
		       e[i].p[0],e[i].p[1],e[i].p[2]);
		Vector nl=e[i].n;
		if (mag2(nl)>.99)
			for (int c=0;c<3;c++)
				if (nl[c] && fabs(nl[c])<1e-5) nl[c]=0;
		if (!nl.iszero())
			output(binary,'n',nl[0],nl[1],nl[2]);
	}
}

void WA3dStream::writeComment(const char* string)
{
	A3dElem e(A3dElem::TComment,0);
	for (;;) {
		const char* p=strchr(string,'\n');
		if (!p) break;
		char s[200];
		assertx(p-string<static_cast<int>(sizeof(s)));
		s[0]=0; strncat(s,string,p-string);
		e.setcomment(s);
		write(e);
		string=p+1;
	}
	e.setcomment(string);
	write(e);
}

void WA3dStream::writeEndObject(int binary, double f0, double f1)
{
	A3dElem e(A3dElem::TEndObject,binary);
	e.f(0)=f0; e.f(1)=f1; e.f(2)=0;
	write(e);
}

void WA3dStream::writeClearObject(int binary, double f0, double f1)
{
	A3dElem e(A3dElem::TEditObject,binary);
	e.f(0)=f0; e.f(1)=f1; e.f(2)=0;
	write(e);
}

void WA3dStream::writeEndFrame(int binary)
{
	A3dElem e(A3dElem::TEndFrame,binary);
	e.f(0)=0; e.f(1)=0; e.f(2)=0;
	write(e);
}

//*** WSA3dStream

WSA3dStream::WSA3dStream(ostream& pos) : ios(pos) { }

WSA3dStream::~WSA3dStream() { flush(); }

void WSA3dStream::output(int binary, int type, double f1, double f2, double f3)
{
	if (binary) {
		char s[16];
		double lf1=f1,lf2=f2,lf3=f3; // make sure not double
		unsigned short t=type;
		s[0]=A3dStream::BINARYCODE;
		s[1]=0;
		ShortToStd(&t,&s[2]);
		FloatToStd(&lf1,&s[4]);
		FloatToStd(&lf2,&s[8]);
		FloatToStd(&lf3,&s[12]);
		ios.write(s,16);
	} else {
		if (oldformat) {
			ios << f1 << " " << f2 << " " << f3 <<
				" " << char(type) << '\n';
		} else {
			ios << char(type) << " " <<
				f1 << " " << f2 << " " << f3 << '\n';
		}
	}
	assertx(ios);
}

void WSA3dStream::outputcomment(const char* s)
{
	ios << "#" << s << "\n";
}

void WSA3dStream::blankline()
{
	if (oldformat) return;
	ios << '\n';
}

void WSA3dStream::flush()
{
	ios.flush();
}

ostream& WSA3dStream::os()
{
	return ios;
}
