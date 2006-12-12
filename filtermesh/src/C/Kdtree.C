// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#include "Kdtree.h"

int KDTREE_STATS=GetenvValue("KDTREESTATS");
Stat SKDsearchnel("SKDsearchnel",KDTREE_STATS);	/* UNSAFE! */
ALLOCATEPOOL(Kdentry);
