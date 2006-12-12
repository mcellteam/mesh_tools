// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef FileIO_h
#define FileIO_h
#include <iostream>

class WFile {
  public:
	WFile(const char *filename); // supports "-", ".Z", "|..."
	~WFile();
	std::ostream& operator()() const;
  private:
	std::ostream* ios;
	FILE* file;
	DISABLECOPY(WFile);
};

class RFile {
  public:
	RFile(const char *filename); // supports "-", ".Z", "...|"
	~RFile();
	std::istream& operator()() const;
  private:
	std::istream* iis;
	FILE* file;
	DISABLECOPY(RFile);
};

#endif
