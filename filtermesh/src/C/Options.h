// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Options_h
#define Options_h

#include "Queue.h"
#include <iostream>

class Options {
  public:
	// sample names (when multiple Options):  "Main", "Hx", "Starbase"
	Options(const char* pname="");
	~Options();
	
	// Flag:
	void f(const char* st, int& arg,		const char* doc="");
	void f(const char* st, const char*& arg,	const char* doc="");
	
	// Parameter:
	void p(const char* st, int& arg,		const char* doc="");
	void p(const char* st, int* argp, int narg,	const char* doc="");
	void p(const char* st, double& arg,		const char* doc="");
	void p(const char* st, double* argp, int narg,	const char* doc="");
	void p(const char* st, char const*& arg,	const char* doc="");
	void p(const char* st, char const** argp, int narg,
	       const char* doc="");
	// Comment:
	void c(const char* st="", const char* doc="");
	
	// Do not skip over unrecognized args
	void allmustparse();
	
	// Function f is given (argc,argv) starting with the matched option.
	// It must return the number of arguments removed (incl. option).
	// A return value less than 0 indicates a parsing error.
	typedef int (*PARSEF)(int argc, char const** argv);
	void p(const char* st, PARSEF f,		const char* doc="");

	// Function for handling a flag
	typedef void (*HANDLEF)();
	void p(const char* st, HANDLEF f,		const char* doc="");
	
	// Only one parse() may be active at any one time.
	// Returns 0 if either not enough arguments after last flag.
	//        or if user-defined flag handler returns value<0 .
	// Returns 1 upon success.
	// Note: argv[0] is unchanged
	// Note: changes argv[] in place!
	int parse(int& argc, char const** argv);

	int noflags(int& argc, char const** argv); // ret 0 if any flags remain (--)
	int noargs(int argc, char const** argv);   // ret 0 if any args remain
	void problem(int argc, char const** argv); // call if cannot parse

	void print_help(std::ostream& s) const;

  public:			// for stupid cxx
	struct option {
		const char* str;
		PARSEF f;
		void* argp;
		int narg;	// number of arguments if f is built-in
		const char* doc;
	};
  private:
	static Options* curopt;
	const char* name;	// name of options (hx,starbase)
	char ename[210];	// execution name: argv[0]+name
	Queue<option*> qoptions;
	void* curarg;
	int curnarg;
	int iallmustparse;
	option* match(const char* s);
	friend int fint(int argc, char const** argv);
	friend int ffloat(int argc, char const** argv);
	friend int fstring(int argc, char const** argv);
	friend int parsequestion(int argc, char const** argv);
	void iadd(const char* st, PARSEF f, void* argp, int narg,
		  const char* doc);
	DISABLECOPY(Options);
};

#define OPTSF(opts,var) opts.f("-" #var, var, ": set flag '" #var "'")

#define OPTSFC(opts,var,comment) opts.f("-" #var, var, comment)

#define OPTSP(opts,var) opts.p("-" #var, var, "v : set variable '" #var "'")

#define OPTSPC(opts,var,comment) opts.p("-" #var, var, comment)

#define OPTSA(opts,var,n) opts.p("-" #var, var, n, \
#n "_args : set variables for '" #var "'")

#define OPTSAC(opts,var,n,comment) opts.p("-" #var, var, n, comment)

#define OPTSD(opts,var) opts.p("-" #var, do##var, ": execute do" #var)

#define OPTSDC(opts,var,comment) opts.p("-" #var, do##var, comment)

#endif
