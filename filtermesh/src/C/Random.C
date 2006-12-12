// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Hh.h"
#include "Random.h"

#if defined(__DECCXX) || defined(__alpha) || __GNUG__==2 && __GNUC_MINOR__==5
extern "C" {
	char* setstate(char*);
	char* initstate(unsigned, char*, int);
	long random(void);
	void srandom(int);
}
#endif

Random Random::G;		// default seed

Random::Random(int seed)
{
	setseed(seed);
	uniffactor=1/pow(2,31);
	// bug in GNUG 2.5.8: int/double/double(-NGAUSS) is huge 4.29e9
	gaussoffset=-int(NGAUSS)*.5*pow(2,31);
	gaussfactor=mysqrt(12)/pow(2,31)/mysqrt(double(NGAUSS));
}

Random::~Random() { }

int Random::getint()
{
	char* oldstate=setstate(reinterpret_cast<char *>(state));
	int i=int(random());
	setstate(oldstate);
	return i;
}

double Random::unif()
{
	char* oldstate=setstate(reinterpret_cast<char *>(state));
	double f=random()*uniffactor;
	setstate(oldstate);
	return f;
}

double Random::dunif()
{
	char* oldstate=setstate(reinterpret_cast<char *>(state));
	double f=random()*uniffactor;
	setstate(oldstate);
	return f;
}

double Random::gauss()
{
	char* oldstate=setstate(reinterpret_cast<char *>(state));
	double acc=0;
	for (int i=0;i<NGAUSS;i++) acc+=random();
	double f=(acc+gaussoffset)*gaussfactor;
	setstate(oldstate);
	return f;
}

double Random::dgauss()
{
	char* oldstate=setstate(reinterpret_cast<char *>(state));
	double acc=0;
	for (int i=0;i<NGAUSS;i++) acc+=random();
	double f=(acc+gaussoffset)*gaussfactor;
	setstate(oldstate);
	return f;
}

void Random::showstate() const
{
	SHOWS("Random state is {");
	for (int i=0;i<SIZE;i++) SHOWF("%3d %d\n",i,state[i]);
	SHOWS("} End Random state");
}

void Random::setseed(int seed)
{
	assertx(sizeof(int)==4);
	setstate(initstate(seed,reinterpret_cast<char *>(state),SIZE*sizeof(int)));
}
