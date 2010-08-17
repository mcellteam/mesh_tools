// Author: Justin Kinney
// Date: Feb 2009

#ifndef RECONSTRUCT2CONTOURTILER_H
#define RECONSTRUCT2CONTOURTILER_H 1

#include <vector>

class Contour;
class Histogram;
class Object;
class Point;
class SplinePoint;

#include "container.h"
#include "controls.h"

bool distinguishable (double a,double b,double epsilon);
bool distinguishable (double a,double b);

#endif
