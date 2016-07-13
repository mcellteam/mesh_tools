#include "deviations.h"

#include <math.h>

Deviations::Deviations (void)
  :center(),radius(),distances()
{
}

void Deviations::initialize (const std::vector<Point> & raw_points)
{
  center.reserve(raw_points.size());
  radius.reserve(raw_points.size());
  c_p_iterator i = raw_points.begin();
  c_p_iterator j = i;
  j++;
  while (j!=raw_points.end())
  {
    // compute center
    double x = 0.5*((*i).getX()+(*j).getX());
    double y = 0.5*((*i).getY()+(*j).getY());
    center.push_back(Point(x,y));
    // compute radius
    double delx = (*i).getX()-(*j).getX();
    double dely = (*i).getY()-(*j).getY();
    double r = sqrt(delx*delx + dely*dely);
    radius.push_back(r);
    // advance iterators
    i++;
    j++;
  }
  i = raw_points.end();
  i--;
  j = raw_points.begin();
  // compute center
  double x = 0.5*((*i).getX()+(*j).getX());
  double y = 0.5*((*i).getY()+(*j).getY());
  center.push_back(Point(x,y));
  // compute radius
  double delx = (*i).getX()-(*j).getX();
  double dely = (*i).getY()-(*j).getY();
  double r = sqrt(delx*delx + dely*dely);
  radius.push_back(r);
}

double Deviations::computeDeviation (const double & x,const double & y)
{
  std::vector<double> d;
  d.reserve(center.size());
  double min_deviation = 1E30;
  int min_index = -1;
  // for each line segment
  for (int i=0;i<static_cast<int>(center.size());i++)
  {
    double dx = center[i].getX()-x;
    double dy = center[i].getY()-y;
    d.push_back(sqrt(dx*dx+dy*dy));
    // if sample point is inside circumscribing line segment circle
    if (d[i]<=radius[i])
    {
      // compute min distance from sample point to line segment
    }
  }
  // for each line segment
  for (int i=0;i<static_cast<int>(center.size());i++)
  {
    // if sample point is outside circumscribing line segment circle
    if (d[i]>radius[i])
    {
      double min = d[i]-radius[i];
      if (min<=min_deviation)
      {
        // compute min distance from sample point to line segment
      }
    }
  }

  //  > identify subset of line segments to check
  //  > for each line segment
  //    + compute perpendicular distance between sample point and line segment
  //    + keep minimum distance
  //  > if minimum distance is larger than threshold
  //    + duplicate control point (which one?)
  //
  // NOTE THIS RETURN VALUE IS A PLACE HOLDER
  return min_index*x;
}

void Deviations::calculateDeviations (const std::vector<Point> & spline_samples)
{
  distances.reserve(spline_samples.size());
  // for each sample point
  for (c_p_iterator i = spline_samples.begin();i!=spline_samples.end();i++)
  {
    distances.push_back(computeDeviation((*i).getX(),(*i).getY()));
  }
}

