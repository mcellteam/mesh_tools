// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Timer_h
#define Timer_h

#include <iostream>

#if 0
{
	NEST { TIMER(__atimer2); stat; }
	NEST { TIMER(atimer2); stats; ETIMER(atimer2); ...; }
}
#endif

// GetenvValue("SHOWTIMES")==-1 -> all->Noprint
// GetenvValue("SHOWTIMES")==1 -> Diagnostic->Default, Abbrev->Default

class Timer {
  public:
	enum mode { Default, Diagnostic, Abbrev, Noprint };
	// Default:	SHOWDF every time
	// Diagnostic:	SHOWF every time
	// Abbrev:	SHOWF first time
	// Noprint:	never print, do not keep stats
	Timer(const char* pname=0, mode=Noprint); // zero but not started
	~Timer();
	void setname(const char* pname);
	void setmode(mode pmode);
	void zero();
	void start();
	void stop();
	void terminate();
	const char* name() const;
	double real() const;
	double user() const;
	friend std::ostream& operator<<(std::ostream& s, const Timer& t);
	static int showtimes(int val);	// ret pval
  private:
	const char* iname;
	mode imode;
	int started;
	int tmspu, tmsps;	// process{user,system}
	double ireal;
	static int ishow;
	static int firstprint;
	DISABLECOPY(Timer);
};

#define TIMERN(id) Timer_##id
#define TIMER(id) Timer TIMERN(id)(#id,Timer::Default); TIMERN(id).start()
#define DTIMER(id) TIMER(id); TIMERN(id).setmode(Timer::Diagnostic);
#define ATIMER(id) TIMER(id); TIMERN(id).setmode(Timer::Abbrev);
#define ETIMER(id) TIMERN(id).terminate()

#endif
