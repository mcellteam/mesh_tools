// Author: Justin Kinney
// Date: Sep 2008

#ifndef DISTANCE_H
#define DISTANCE_H 1

#include "meshheal.h"

#include "vertex.h"

class Distance
{
public:
  double d; // distance between vertices with indices vA and vB
  Vertex *vA;
  Vertex *vB;
  bool deleteme;
  Distance(double,void_list*,void_list*);
};

#endif    
