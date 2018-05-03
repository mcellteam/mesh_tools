// Author: Justin Kinney
// Date: Feb 2009

#ifndef POINT_H
#define POINT_H 1

#include <math.h>
#include <stdio.h>

class Point
{
private:
  double x,y;
public:
  Point (void);
  Point (char const * str, int const dim, double * const transform);
  Point (double xval, double yval);
  double getX (void) const { return x; }
  double getY (void) const { return y; }
public:
  void setX (double i)
  {
    x=i;
  }
  void setY (double i)
  {
    y=i;
  }
  void print (char * const line) const
  {
    sprintf(line,"%g %g",x,y);
  }
  Point operator- ( const Point & p) const
  {
    return Point( x - p.x,
                  y - p.y );
  }
  double dot ( const Point & p) const
  {
    return x*p.x+y*p.y;
  }
  double distance_squared ( const Point & p) const
  {
    double a = x - p.x;
    double b = y - p.y;
    return a*a+b*b;
  }
  double length (void) const
  {
    return sqrt(x*x+y*y);
  }
  Point& operator*= ( double a )
  {
    x *= a;
    y *= a;
    return *this;
  }
};

#endif
