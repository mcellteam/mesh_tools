// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Stat.h"
#include "Queue.h"

using std::ostream;
using std::ofstream;

static struct Statclass {
	// destructor here may not assume Stats in queue still exist!
	void flush();
	Queue<Stat*> queue;
} Stats;

void Statclass::flush()
{
	if (queue.empty()) return;
	SHOWDF("Summary of statistics:\n");
	ForQueue(queue,Stat*,stat) {
		stat->terminate();
	} EndFor;
	queue.clear();
}

void FlushStats() { Stats.flush(); }

Stat::Stat(const char* pname, int pprint, int isstatic)
: iname(pname), iprint(pprint), fos(0)
{
	zero();
	if (iname && GetenvValue("STATFILES"))
		fos=new ofstream(hform("Statfile.%.90s",iname?iname:""));
	if (isstatic) Stats.queue.enqueue(this);
}

Stat::~Stat() { terminate(); }

void Stat::terminate()
{
	if (iprint && num()) SHOWDF("%s",namestring());
	iprint=0;
	delete fos,fos=0;	// may be 0, fos->close() called by destructor
}

void Stat::zero()
{
	assertw(!fos);		// just a warning
	n=0;
	isum=isum2=0;
	imin=1e30;
	imax=-1e30;
}

void Stat::setname(const char* pname)
{
	iname=pname;
}

Stat& Stat::operator+=(const Stat& st)
{
	assertw(!fos);		// just a warning
	n+=st.n;
	if (st.imin<imin) imin=st.imin;
	if (st.imax>imax) imax=st.imax;
	isum+=st.isum;
	isum2+=st.isum2;
	return *this;
}

double Stat::var() const
{
	if (assertw1(n>1)) return 0;
	return ::max((isum2-isum*isum/n)/(double(n-1)),0.);
}

const char* Stat::string() const
{
	double tavg=n>0?avg():0, tsdv=n>1?sdv():0;
	return hform("(%-7i)%11g:%-11g av=%-11g sd=%g",
		     n,imin,imax,tavg,tsdv);
}

const char* Stat::namestring() const
{
	return hform("%-12.20s%s\n",hform("%.19s:",iname?iname:"?"),string());
}

void Stat::output(double f) const
{
	*fos << f << '\n';
}

ostream& operator<<(ostream& s, const Stat& st)
{
	return s << st.namestring();
}
