// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Principal_h
#define Principal_h

#include "Geometry.h"

// Given n points pa[], compute the principal component frame f,
// and the (redundant) eigenvalues eimag[3].
// The frame f is guaranteed to be invertible and orthogonal, as the axis
// will always have non-zero (albeit very small) lengths.
// The values eimag[] will be zero for the axes that should be zero.
// The frame f is also guaranteed to be right-handed.
extern void PrincipalComponents(const Point pa[], int n, Frame& f,
				double eimag[3]);

#endif
