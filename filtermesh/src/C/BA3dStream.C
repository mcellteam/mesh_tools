// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include <string>
#include <sstream>

#include "Hh.h"
#include "BA3dStream.h"

using std::string;
using std::istringstream;
using std::ostringstream;

//*** RBA3dStream

RBA3dStream::RBA3dStream(RBuffer& b) : buf(b) { }

RBA3dStream::~RBA3dStream() { }

int RBA3dStream::readline(int& binary, int& type, double f[3],
			  const char*& comment)
{
        int i;
	// read first character.  Read past optional newlines
//	for (int i=0;;i++) {
	for (i=0;;i++) {
		assertx(i<buf.num());
		if (buf[i]!='\n') break;
	}
	if (i) buf.extract(i);
	if ((binary=buf[0]==A3dStream::BINARYCODE) != 0) { // binary
		assertx(buf.num()>=16);
		type=buf.bgetshort(2);
		for (int i=0;i<3;i++) f[i]=buf.bgetfloat(i*4+4);
		buf.extract(16);
		return 1;
	}
	type=buf[0];
	static char str[A3dStream::LINELENGTH];	// not reentrant
	assertx(buf.extractline(str, sizeof(str)));
	if (type==A3dElem::TComment) {
		str[strlen(str)-1]=0; // remove trailing '\n'
		comment=&str[1];
		return 1;
	}
	assertx(str[1]==' ');
        string the_line(&str[2]);
	istringstream istr(the_line); // ok since trailing '\n'
	istr >> f[0] >> f[1] >> f[2];
	assertx(istr);
	return 1;
}

int RBA3dStream::recognize() const // -1=parse_err, 0=no, 1=partial, 2=yes
{
        int i;
	// Skip leading newlines
//	for (int i=0;i<buf.num();i++)
	for (i=0;i<buf.num();i++)
		if (buf[i]!='\n') break;
	if (i==buf.num()) return 0;
	char c=buf[i];
	int havespace=i+1<buf.num() && buf[i+1]==' ';
	if (!(c==A3dElem::TComment ||
	      c==A3dStream::BINARYCODE ||
	      havespace &&
	      (A3dElem::commandtype(A3dElem::Type(c)) || c==A3dElem::TPoint ||
	       c==A3dElem::TPolygon || c==A3dElem::TPolyline ||
	       A3dElem::statustype(A3dElem::Type(c)) || c=='n'))) return 0;
	for (;;) {
		if (i>=buf.num()) return 1;
		c=buf[i];
		if (c=='\n') { i++; continue; }
		if (c==A3dStream::BINARYCODE) { // binary record
			if (i+16>buf.num()) return 1;
			c=buf.bgetshort(i+2);
			i+=16;
		} else {
			for (i++;;) {
				if (i>=buf.num()) return 1;
				if (buf[i++]=='\n') break;
			}
		}
		if (c==A3dElem::TComment ||
		    A3dElem::commandtype(A3dElem::Type(c)) ||
		    c==A3dElem::TPoint || c=='E') return 2;
		if (!(c==A3dElem::TPolygon || c==A3dElem::TPolyline ||
		      A3dElem::statustype(A3dElem::Type(c)) ||
		      c=='n' || c=='v'))
			return -1;
	}
}

//*** WBA3dStream

WBA3dStream::WBA3dStream(WBuffer& b) : buf(b) { }

WBA3dStream::~WBA3dStream() { flush(); }

void WBA3dStream::output(int binary, int type, double f1, double f2, double f3)
{
	if (binary) {
		buf.bput(char(A3dStream::BINARYCODE));
		buf.bput(char(0));
		buf.bput(short(type));
		buf.bput(f1); buf.bput(f2); buf.bput(f3);
	} else {
		ostringstream ostr;
		if (oldformat) {
			ostr << f1 << " " << f2 << " " << f3 <<
				" " << char(type) << '\n';
		} else {
			ostr << char(type) << " " <<
				f1 << " " << f2 << " " << f3 << '\n';
		}
		string s(ostr.str());
		buf.bput(s.c_str(), s.length());
	}
}

void WBA3dStream::outputcomment(const char* s)
{
	buf.bput("#",1);
	buf.bput(s,strlen(s));
	buf.bput("\n",1);
}

void WBA3dStream::blankline()
{
	if (oldformat) return;
	buf.bput('\n');
}

void WBA3dStream::flush()
{
	buf.flush();
}
