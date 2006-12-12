// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Timer.h"
#include "Map.h"
#include "Stat.h"
#include "Queue.h"

#include <sys/time.h>		// for gettimeofday()
#include <ctime>		// for gettimeofday()
#include <sys/times.h>		// for times()
#include <iostream>

using std::ostream;

extern "C" {
#ifdef __DECCXX
	int gettimeofday(struct timeval*, struct timezone*);
	clock_t times(struct tms* buffer);
#elif __GNUG__==2 && __GNUC_MINOR__==5
	int gettimeofday(struct timeval*, struct timezone*);
#endif
}

// See time.h, sysconf(1)
#if defined(CLK_TCK)
const int NTICKS=CLK_TCK;
#elif defined(__alpha) && !defined(__osf__)
const int NTICKS=1000000;
#elif defined(sgi)
const int NTICKS=100;
#else
const int NTICKS=60;
#endif

// Overhead for use of Timer:
//  on sgi:  timer off by +39e-6 sec, overhead is  81e-6 sec
//  on mips: timer off by +82e-6 sec, overhead is 168e-6 sec

static struct Timerclass {
	~Timerclass() { flush(); }
	void flush();
	Map<const char*,Stat*> m;
	Queue<Stat*> queue;	// in order encountered during execution
} Timers;

void Timerclass::flush()
{
	int havesome=0;
	ForQueue(queue,Stat*,stat) {
		if (stat->num()>1) havesome=1;
	} EndFor;
	if (havesome)
		SHOWDF("Summary of timers (cpu=%s host=%s):\n",
		       getenv("CPU")?getenv("CPU"):"?",
		       getenv("LHOST")?getenv("LHOST"):
		       getenv("HOST")?getenv("HOST"):"?");
	ForQueue(queue,Stat*,stat) {
		if (!havesome) continue;
		int n=stat->num();
		SHOWDF(" %-20.20s(%-6i)%s:%s av=%9.2lf   sum=%9.2lf\n",
		       hform("%.19s:",stat->name()), n,
		       n>1?hform("%8.2lf",stat->min()):"        ",
		       n>1?hform("%-8.2lf",stat->max()):"        ",
		       stat->avg(), stat->sum());
	} EndFor;
	ForQueue(queue,Stat*,stat) {
		delete stat;
	} EndFor;
	m.clear();
	queue.clear();
}

void FlushTimers() { Timers.flush(); }

int Timer::ishow=getenv("SHOWTIMES")?GetenvValue("SHOWTIMES"):0;
int Timer::firstprint=1;

Timer::Timer(const char* pname, mode pmode)
: iname(pname), imode(pmode)
{
	if (iname && ishow>=2)
		SHOWF(" (%-20.20s started)\n",hform("%.19s:",iname));
	zero();
}

Timer::~Timer() { terminate(); }

void Timer::terminate()
{
	if (ishow>0 && imode!=Noprint) imode=Default;
	if (ishow<0) imode=Noprint;
	mode cmode=imode;
	imode=Noprint;
	if (cmode==Noprint || assertw1(iname)) return;
	if (!assertw1(started)) stop();
	double u=user();
	Stat* stat=Timers.m.retrieve(iname);
	if (!stat) {
		stat=new Stat(iname,0);
		Timers.m.enter(iname,stat);
		Timers.queue.enqueue(stat);
	}
	stat->enter(u);
	if (cmode==Abbrev && stat->num()>1) return;
	if (firstprint) {
		firstprint=0;
		SHOWFF("(Timing on cpu=%s host=%s)\n",
		       getenv("CPU")?getenv("CPU"):"?",
		       getenv("LHOST")?getenv("LHOST"):
		       getenv("HOST")?getenv("HOST"):"?");
	}
	const char* s=hform(" (%-20.20s %8.2lf)\n", hform("%.19s:",iname),u);
	if (cmode==Default) {
		SHOWDF("%s",s);
	} else {
		SHOWF("# %s",s);
	}
}

int Timer::showtimes(int val)
{
	int pval=ishow; ishow=val; return pval;
}

void Timer::setname(const char* pname)  { iname=pname; }

void Timer::setmode(mode pmode) { imode=pmode; }

void Timer::zero()
{
	started=0;
	tmspu=tmsps=0;
	ireal=0;
}

void Timer::start()
{
	struct tms buf2; struct timeval ti; struct timezone tz;
	assertx(!started++);
	(void)times(&buf2);
	tmspu -= static_cast<int>(buf2.tms_utime);
        tmsps -= static_cast<int>(buf2.tms_stime);
	assertw(!gettimeofday(&ti,&tz));
	ireal-=double(ti.tv_sec)+ti.tv_usec/1e6;
}

void Timer::stop()
{
	struct tms buf2; struct timeval ti; struct timezone tz;
	assertx(!--started);
	(void)times(&buf2);
	tmspu += static_cast<int>(buf2.tms_utime);
        tmsps += static_cast<int>(buf2.tms_stime);
	assertw(!gettimeofday(&ti,&tz));
	ireal+=double(ti.tv_sec)+ti.tv_usec/1e6;
}

const char* Timer::name() const { return iname; }

double Timer::real() const { assertx(!started); return ireal; }

double Timer::user() const
{
	assertx(!started);
	return double(tmspu+tmsps)/NTICKS;
}

ostream& operator<<(ostream& s, const Timer& t)
{
	s << "Timer " << (t.iname?t.iname:"?") << ": ";
	if (t.started) return s << "running\n";
	return s << hform("real=%-.6lf  u=%-5.2lf  s=%-5.2lf\n",
			  t.ireal,double(t.tmspu)/NTICKS,double(t.tmsps)/NTICKS);
}
