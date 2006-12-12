// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Options.h"
#include <iostream>

using std::cerr;
using std::ostream;

Options* Options::curopt=0;

Options::Options(const char* pname)
: name(pname), iallmustparse(0)
{
	assertx(name);
	ename[0]=0; strncat(ename,name,100);
	iadd("-h",parsequestion,0,0,": print available options");
}

Options::~Options()
{
	while (!qoptions.empty()) delete qoptions.dequeue();
}

void Options::iadd(const char* st, PARSEF f, void* argp, int narg,
		   const char* doc)
{
	option* o=new option;
	o->str=st; o->f=f; o->argp=argp; o->narg=narg; o->doc=doc;
	qoptions.enqueue(o);
	if (*doc && narg>=0) {	// enforce my "conventions".
		if (narg>0 && assertw(*doc!=':')) SHOWN(st);
		if (!narg  && assertw(*doc==':')) SHOWN(st);
	}
}

Options::option* Options::match(const char* s)
{
	option* omatch=0;
	int nmatches=0,minlfound=999;
	int ls=strlen(s);
	ForQueue(qoptions,option*,o) {
		if (!o->f) continue;
		int lo=strlen(o->str);
		const char* t=strchr(o->str,'[');
		int minfit=t?t-o->str:0;
		int nchar=min(max(ls,2),lo);
		if (minfit) nchar=minfit;
		if (!strncmp(s,o->str,nchar)) {
			nmatches++;
			if (nchar==ls && lo==ls) { // exact match
				// bad Options if 2 exact matches
				assertx(minlfound);
				minlfound=0,omatch=o;
			}
			if (lo<minlfound) minlfound=lo,omatch=o;
		}
	} EndFor;
	if (!minlfound) nmatches=1;
	if (nmatches>1)
		SHOWF("Argument '%s' ambiguous, '%s' assumed\n",s,omatch->str);
	return omatch;
}

void Options::allmustparse() { iallmustparse=1; }

int Options::parse(int& argc, char const** argv)
{
	int nargc=1,skip=0;
	int i=1;
	ename[0]=0; strncat(ename,argv[0],100);
	if (name[0]) strcat(ename," "),strncat(ename,name,100);
	for (;i<argc;) {
		if (!strcmp(argv[i],"--")) skip=1;
		option* o=match(argv[i]);
		if (!o || skip) { // not found
			if (iallmustparse) {
				SHOWF("%s: Can't parse '%s'\n",ename,argv[i]);
				return 0;
			}
			argv[nargc++]=argv[i++]; // may be aliased
			continue;
		}
		curopt=this;
		curarg=o->argp;
		curnarg=o->narg;
		if (curnarg>=0 && curnarg+1>argc-i) {
			cerr << ename << ": Not enough options after '" <<
				argv[i] << "' (it takes " <<
				curnarg << " args)\n";
			return 0;
		}
		int nread=1;
		if (curnarg==-2) {
			(*HANDLEF(o->f))();
		} else {
			nread=(*o->f)(argc-i,&argv[i]);
		}
		if (!nread && !curnarg) {
			// allow simple flag '-h' to not recognize itself
			argv[nargc++]=argv[i++];
			if (iallmustparse) return 0;
			continue;
		}
		if (nread<1) {
			SHOWF("%s: Error parsing '%s'\n",ename,argv[i]);
			return 0;
		}
		i+=nread;
		assertx(i<=argc);
	}
	assertx(nargc<=argc);
	argv[nargc]=0;
	argc=nargc;
	return 1;
}

void Options::print_help(ostream& s) const
{
	s << ename << " Options:\n";
	ForQueue(qoptions,option*,o) {
		option& op=*o;
		const char* p=strchr(op.doc,':');
		int pref=p?p-op.doc:0;
		const char* s1; const char* s2;
		s1=hform("%s%s",op.str,
			 strchr(op.str,'[') && !strchr(op.str,']')?"]":"");
		s2=hform("%s %.*s",s1,pref,op.doc);
		s << hform("  %-20s %s\n",s2,op.doc+pref);
	} EndFor;
}

static int check_integer(const char* s)
{
	for (;*s;s++)
		if (!(*s>='0' && *s<='9') && *s!='-' && *s!='+')
		return 0;
	return 1;
}

static int check_float(const char* s)
{
	for (;*s;s++)
		if (!(*s>='0' && *s<='9') &&
		    *s!='-' && *s!='+' && *s!='.' && *s!='e')
			return 0;
	return 1;
}

int fint(int argc, char const** argv)
{
	int* argp = reinterpret_cast<int *>(Options::curopt->curarg);
	int n=Options::curopt->curnarg;
	if (!n) { (*argp)++; return 1; } // set a flag variable
	argv++; argc--;
	for (int i=0;i<n;i++) {
		if (!check_integer(argv[0])) {
			SHOWF("Argument not integer: '%s'\n",argv[0]);
			return -1;
		}
		argp[i]=atoi(argv[0]);
		argv++; argc--;
	}
	return n+1;
}

int ffloat(int argc, char const** argv)
{
	double* argp = reinterpret_cast<double *>(Options::curopt->curarg);
	int n=Options::curopt->curnarg;
	argv++; argc--;
	for (int i=0;i<n;i++) {
		if (!check_float(argv[0])) {
			SHOWF("Argument not double: '%s'\n",argv[0]);
			return -1;
		}
		argp[i]=atof(argv[0]);
		argv++; argc--;
	}
	return n+1;
}

int fstring(int argc, char const** argv)
{
	char const** argp = reinterpret_cast<char const **>(Options::curopt->curarg);
	int n=Options::curopt->curnarg;
	if (!n) { *argp=argv[0]; return 1; } // set a flag string
	argv++; argc--;
	for (int i=0;i<n;i++) {
		argp[i]=argv[0];
		argv++; argc--;
	}
	return n+1;
}

int parsequestion(int, char const**)	// argc, argv
{
	Options::curopt->print_help(cerr);
	return 0;
}

int Options::noflags(int& argc, char const** argv)
{
	int nargc=1,ok=0;
	for (int i=1;i<argc;i++) {
		if (!strcmp(argv[i],"--")) {
			ok=1;	// and do not copy into final argv
		} else if (!ok && argv[i][0]=='-') {
			if (strcmp(argv[i],"-h")) {
				cerr << ename << ": Flag '" << argv[i] <<
					"' not recognized\n";
			}
			return 0;
		} else {
			argv[nargc++]=argv[i]; // may be aliased
		}
	}
	argv[nargc]=0;
	argc=nargc;
	return 1;
}

int Options::noargs(int argc, char const** argv)
{
	if (argc==1) return 1;
	if (strcmp(argv[1],"-h")) {
		cerr << ename << ": Argument '" << argv[1] <<
			"' not recognized\n";
	}
	return 0;
}

void Options::problem(int argc, char const** argv)
{
	int alreadyprinted=0;
	for (int i=1;i<argc;i++)
		if (!strcmp(argv[i],"-h")) alreadyprinted=1;
	if (!alreadyprinted) print_help(cerr);
}

void Options::f(const char* st, int& arg, const char* doc)
{
	iadd(st, fint, &arg, 0, doc);
}

void Options::f(const char* st, const char*& arg, const char* doc)
{
	iadd(st, fstring, &arg, 0, doc);
}

void Options::p(const char* st, int& arg, const char* doc)
{
	iadd(st, fint, &arg, 1, doc);
}

void Options::p(const char* st, int* argp, int narg, const char* doc)
{
	assertx(narg>0);
	iadd(st, fint, argp, narg, doc);
}

void Options::p(const char* st, double& arg, const char* doc)
{
	iadd(st, ffloat, &arg, 1, doc);
}

void Options::p(const char* st, double* argp, int narg, const char* doc)
{
	assertx(narg>0);
	iadd(st, ffloat, argp, narg, doc);
}

void Options::p(const char* st, char const*& arg, const char* doc)
{
	iadd(st, fstring, &arg, 1, doc);
}

void Options::p(const char* st, char const** argp, int narg,
		       const char* doc)
{
	assertx(narg>0);
	iadd(st, fstring, argp, narg, doc);
}

void Options::p(const char* st, PARSEF f, const char* doc)
{
	iadd(st, f, 0, -1, doc);
}

void Options::p(const char* st, HANDLEF f, const char* doc)
{
	iadd(st, PARSEF(f), 0, -2, doc);
}

void Options::c(const char* st, const char* doc)
{
	iadd(st,0,0,0,doc);
}
