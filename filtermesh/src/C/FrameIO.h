// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef FrameIO_h
#define FrameIO_h

#include "Geometry.h"
#include <iostream>

class RBuffer;
class WBuffer;

class FrameIO {
  public:
// input
	static int recognize(RBuffer& b); // ret -1=err, 0=no, 1=partial, 2=yes
	static int read(std::istream& is, Frame& f, int& obn, double& zoom,
			int& bin); // ret is_success
	static int read(RBuffer& b, Frame& f, int& obn, double& zoom,
			int& bin); // ret is_success
// output
	static const char* string(const Frame& f, int obn, double zoom);
	static int write(std::ostream& os, const Frame& f, int obn, double zoom,
			 int bin); // ret is_success
	static int write(WBuffer& b, const Frame& f, int obn, double zoom,
			 int bin); // ret is_success
// special frames
	static int isNaF(const Frame& f);
	static void MakeNaF(Frame& f);
  private:
	static void decode(std::istream& is, Frame& f, int& obn, double& zoom);
};

#endif
