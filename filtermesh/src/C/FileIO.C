// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "FileIO.h"
#include "fdstream.hpp"
#include <sys/types.h>
#include <sys/stat.h>

using std::cin;
using std::cout;
using std::istream;
using std::ifstream;
using std::ostream;
using std::ofstream;

///*** WFile

WFile::WFile(const char* s) : ios(0), file(0)
{
	assertx(s && *s);
	if (s[0]=='|') {
		file=assertv(popen(&s[1],"w"));
		ios=new boost::fdostream(fileno(file));
	} else if (strlen(s)>2 && !strcmp(&s[strlen(s)-2],".Z")) {
		file=assertv(popen(hform("compress >%s",s),"w"));
		ios=new boost::fdostream(fileno(file));
	} else if (strlen(s)>3 && !strcmp(&s[strlen(s)-3],".gz")) {
		file=assertv(popen(hform("gzip >%s",s),"w"));
		ios=new boost::fdostream(fileno(file));
	} else if (!strcmp(s,"-")) {
		ios=&cout;
	} else {
		ios=new ofstream(s);
	}
	assertx(*ios);
}

WFile::~WFile()
{
	if (ios!=&cout) delete ios;
	if (file) assertw(!pclose(file));
}

ostream& WFile::operator()() const
{
	return *ios;
}

//*** RFile

RFile::RFile(const char* s) : iis(0), file(0)
{
	assertx(s && *s);
	if (strlen(s)>2 && s[strlen(s)-1]=='|') {
		char const* t = hform("%s",s);
		file=assertv(popen(t,"r"));
		iis=new boost::fdistream(fileno(file));
		return;
	} else if (strlen(s)>2 && (!strcmp(&s[strlen(s)-2],".Z") ||
			    !strcmp(&s[strlen(s)-2],".z") ||
			    !strcmp(&s[strlen(s)-3],".gz"))) {
		file=assertv(popen(hform("zcat %s",s),"r"));
		iis=new boost::fdistream(fileno(file));
		return;
	} else if (!strcmp(s,"-")) {
		iis=&cin;
		return;
	}
	iis=new ifstream(s);
	if (*iis) return;
	delete iis;
	struct stat dummystat;
	if (!stat(hform("%s.Z",s),&dummystat)) {
		file=assertv(popen(hform("zcat %s.Z",s),"r"));
		iis=new boost::fdistream(fileno(file));
		return;
	}
	if (!stat(hform("%s.gz",s),&dummystat)) {
		file=assertv(popen(hform("zcat %s.gz",s),"r"));
		iis=new boost::fdistream(fileno(file));
		return;
	}
	SHOWN(s);
	assertnever("Could not open file");
}

RFile::~RFile()
{
	if (iis!=&cin) delete iis;
	if (file) assertw(!pclose(file));
}

istream& RFile::operator()() const{
	return *iis;
}
