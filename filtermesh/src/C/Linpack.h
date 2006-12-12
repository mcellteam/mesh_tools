// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#ifndef ANSI
#define ANSI
#define ANSIWASNOTDEFINED
#endif

extern "C" {
#include "linpack.h"
}

#ifdef ANSIWASNOTDEFINED
#undef ANSIWASNOTDEFINED
#undef ANSI
#endif
