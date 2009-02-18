// Author: Justin Kinney
// Date: Feb 2009

#ifndef POINT_H
#define POINT_H 1

class Point
{
public:
  double x,y,z; // current vertex
  double a,b,c,L; // closest point and squared distance
  void add (double j,double k,double l)
  {
    double t=(x-j)*(x-j)+(y-k)*(y-k)+(z-l)*(z-l);
    if (t<L){a=j;b=k;c=l;L=t;}
  }
  Point(double j,double k,double l)
    :x(j),y(k),z(l),a(0.0),b(0.0),c(0.0),L(1e300)
  {
  }
};

#endif
