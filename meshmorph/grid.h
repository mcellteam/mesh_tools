// Author: Justin Kinney
// Date: Sep 2008

#ifndef GRID_H
#define GRID_H 1

#include "meshmorph.h"

struct Grid
{
  vector3 unit_u;
  vector3 unit_v;
  double uv_vert1_u;
  double uv_vert2[2];
  Grid(Face*);
  void computeBarycenter(vector3 &,int,int,Face*);
};

#endif
