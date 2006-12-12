// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Random_h
#define Random_h

class Random {
  public:
	static Random G;
	Random(int seed=1);
	~Random();
	int getint();
	double unif();
	double dunif();
	double gauss();
	double dgauss();
	void setseed(int seed);
  private:
	enum { SIZE=32, NGAUSS=10 };
	int state[SIZE];	// int should be 32 bit
	double uniffactor;
	double gaussoffset,gaussfactor;
	void showstate() const;
	DISABLECOPY(Random);
};

#endif
