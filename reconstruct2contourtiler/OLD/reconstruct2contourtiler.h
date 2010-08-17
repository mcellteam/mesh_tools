// Author: Justin Kinney
// Date: Feb 2009

#ifndef RECONSTRUCT2CONTOURTILER_H
#define RECONSTRUCT2CONTOURTILER_H 1

class Contour;
class Histogram;
class Object;
class Point;
class SplinePoint;

#include "container.h"
#include "controls.h"

class SplinePoint
{
public:
  double t;
  double x;
  double y;
  double r;
  double intfac;
  SplinePoint (void)
  :t(0),x(0),y(0),r(0),intfac(0)
  {
  }
};

#endif
