// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include <string>
#include "Hh.h"
#include "FrameIO.h"
#include "Buffer.h"

#include <iostream>
#include <sstream>

using std::ios;
using std::istream;
using std::ostream;
using std::istringstream;

const int BinaryCode=2;
const int ASCBUFLEN=400;
const int BINBUFLEN=200;

const double BIGF=1e30f;

int FrameIO::recognize(RBuffer& b) // -1=parse_err, 0=no, 1=partial, 2=yes
{
	if (b.num()<2) return 0;
	if (b[0]==BinaryCode) return b.num()>=14*4?2:1; // binary record
	if (!(b[0]==' ' || b[0]=='F' && b[1]==' ')) return 0;
	return b.hasline()?2:1;
}

int FrameIO::read(istream& is, Frame& f, int& obn, double& zoom, int& bin)
{
	int c=is.peek();
	if (c<0) return 0;
	if ((bin=c==BinaryCode) != 0) { // binary
		char s[BINBUFLEN];
		if (!is.read(s,1+1+2+4*3*sizeof(double)+sizeof(double)))
                    return 0;

		unsigned short shortobn;
		StdToShort(&s[2], &shortobn);
		obn = shortobn;
                char const *buffer = s+4;
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 3; ++j)
                        {
				StdToFloat(buffer,&f[i][j]);
                                buffer += sizeof(double);
                        }
		StdToFloat(buffer,&zoom);
		return 1;
	} else {
		decode(is,f,obn,zoom);
		int isgood=is?1:0;
		is.ignore(9999,'\n'); // read past newline
		return isgood;
	}
}

int FrameIO::read(RBuffer& b, Frame& f, int& obn, double& zoom, int& bin)
{
	if (!b.num()) return 0;
	if ((bin=b[0]==BinaryCode) != 0) {	// binary
		if (b.num()<14*4) return 0;
		obn=b.bgetshort(2);
		for (int i=0;i<4;i++)
			for (int j=0;j<3;j++)
				f[i][j]=b.bgetfloat((1+(i*3)+j)*4);
		zoom=b.bgetfloat(13*4);
		b.extract(14*4);
		return 1;
	} else {
		if (b[0]!='F' && b[0]!=' ') return 0;
		char s[ASCBUFLEN];
		if (!b.extractline(s,sizeof(s))) return 0;
		std::string str(s);
		istringstream istr(str); // ok since trailing '\n'
		decode(istr,f,obn,zoom);
		// return istr?1:0;
		int ret=istr?1:0;	// cxx compiler bug?? 
		return ret;
	}
}

int FrameIO::write(ostream& os, const Frame& f, int obn, double zoom, int bin)
{
	if (bin) {
		char s[BINBUFLEN];
		s[0]=BinaryCode;
		s[1]=0;
		assertw(obn>=-32768 && obn<32768);
		unsigned short shortobn=obn;
		ShortToStd(&shortobn,&s[2]);
                char *buffer = s+4;
		for (int i=0;i<4;i++)
			for (int j=0;j<3;j++) {
				double a=f[i][j];
				FloatToStd(&a, buffer);
                                buffer += sizeof(double);
			}
		double fzoom=zoom;
		FloatToStd(&fzoom, buffer);
		os.write(s, 1 + 1 + 2 + 3*4*sizeof(double) + sizeof(double));
	} else {
		os << string(f,obn,zoom);
	}
	return os?1:0;
}

int FrameIO::write(WBuffer& b, const Frame& f, int obn, double zoom, int bin)
{
	if (bin) {
		b.bput(char(BinaryCode));
		b.bput(char(0));
		assertw(obn>=-32768 && obn<32768);
		b.bput(short(obn));
		for (int i=0;i<4;i++)
			for (int j=0;j<3;j++)
				b.bput(f[i][j]);
		b.bput(zoom);
	} else {
		const char* s=string(f,obn,zoom);
		b.bput(s,strlen(s));
	}
	return 1;
}

void FrameIO::decode(istream& is, Frame& f, int& obn, double& zoom)
{
	char c=0;
	is.get(c);		// is>>c would eat up white space
	assertx(c!=BinaryCode);	// unsafe
	if (c!='F' && c!=' ') is.clear(ios::badbit);
	is >> obn;
	f.ident();
	for (int i=0;i<4;i++)
		for (int j=0;j<3;j++)
			is >> f[i][j];
	is >> zoom;
}

const char* FrameIO::string(const Frame& f, int obn, double zoom)
{
	return hform("F %d  %g %g %g  %g %g %g  %g %g %g  %g %g %g  %g\n",
		     obn,
		     f[0][0],f[0][1],f[0][2],
		     f[1][0],f[1][1],f[1][2],
		     f[2][0],f[2][1],f[2][2],
		     f[3][0],f[3][1],f[3][2],
		     zoom);
}

int FrameIO::isNaF(const Frame& f)
{
	return f[0][0]==BIGF;
}

void FrameIO::MakeNaF(Frame& f)
{
	f.ident();
	f[0][0]=BIGF;
}
