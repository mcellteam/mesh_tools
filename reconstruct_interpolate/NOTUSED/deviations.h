// Author: Justin Kinney
// Date: July 2009

#ifndef DEVIATIONS_H
#define DEVIATIONS_H 1

#include <vector>

#include "point.h"

typedef std::vector<Point>::const_iterator  c_p_iterator;

class Deviations
{
private:
  std::vector<Point>     center;
  std::vector<double>    radius;
  std::vector<double> distances; // deviation distances
public:
  Deviations (void);
  void initialize          (const std::vector<Point> & raw_points);
  void calculateDeviations (const std::vector<Point> & spline_samples);
  double computeDeviation  (const double & x,const double & y);
};

#endif
