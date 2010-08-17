// Author: Justin Kinney
// Date: Feb 2009

#ifndef POINT_H
#define POINT_H 1

class Point
{
private:
  double x,y,z;
public:
  Point (char *str,int section,double thickness,double *transform);
  Point (double xval,double yval,double zval);
  double getX (void) const { return x; }
  double getY (void) const { return y; }
  double getZ (void) const { return z; }
};

#endif
