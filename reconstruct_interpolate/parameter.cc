#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "parameter.h"

#include "controls.h"

/** Write path parameter values, spline indices,
 *  and path parameter distances to left-hand side sample point to file.
 * \param[in] filename Output file name.
 */

void Parameter::print (char const * const filename) const
{
  assert(s.size()==spline_index.size());
  assert(s.size()==lhs_s.size());
  FILE *F = fopen(filename,"a");
  if (!F) { printf("Could not open output file %s\n",filename); return; }
  char line[2048];
  sprintf(line,"\n# spline samples = %i\n",static_cast<int>(s.size()));
  fputs(line,F);
  for (int i=0;i<static_cast<int>(s.size());i++)
  {
    sprintf(line,"%.5g %i %.5g\n",s[i],spline_index[i],lhs_s[i]);
    fputs(line,F);
  }
  fclose(F);
}

/** Calculate perpendicular distance from sample point
 *  to nearest raw contour point line segment.
 * \param[in] spline_sample Spline sample point of interest.
 * \param[in] contour_point_1 One end of line segment.
 * \param[in] contour_point_2 Other end of line segment.
 * \return Perpendicular distance from sample_point to line segment.
 */

double Parameter::calcMinDistPtLineSegment (const Point & spline_sample,
                                            const Point & contour_point_1,
                                            const Point & contour_point_2) const
{
  // Adapted from
  // Distance Between Point and Line, Ray, or Line Segment
  // David Eberly
  // Geometric Tools, LLC
  // http://www.geometrictools.com/
  // Copyright 
  // c 1998-2008. All Rights Reserved.
  // Created: March 2, 1999
  // Last Modified: March 1, 2008
  // see doc/DistancePointLine.pdf
  Point  M(contour_point_2-contour_point_1);
  Point PB(spline_sample-contour_point_1);
  double a = M.dot(PB);
  if (a<=0.0)
  {
    return PB.length();
  }
  else
  {
    double b = M.dot(M);
    if (a>=b)
    {
      Point C(PB-M);
      return C.length();
    }
    else
    {
      M *= a/b;
      Point C(PB-M);
      return C.length();
    }
  }
}

/** Determine which control points of splines to duplicate
 *  to reduce deviation between spline samples
 *  and linearly-interpolated raw contour points.
 * \param[in] raw_points All original points in contour.
 * \param[out] raw_points_to_duplicate Raw points that
 *  should be duplicated as control points.
 */

void Parameter::getRawPointsToDuplicate (const vec_p & raw_points,
                                         vec_i & raw_points_to_duplicate) const
{
  assert(deviations.size()==closest_raw.size());
  raw_points_to_duplicate.clear();
  double deviation_threshold = Controls::instance().getDeviationThreshold();
  int num_raw_points = static_cast<int>(raw_points.size());
  // for each deviation
  for (int i=0;i<static_cast<int>(deviations.size());i++)
  {
    if (deviations[i]>deviation_threshold)
    {
      int closest = closest_raw[i];
      std::vector<int>::iterator j;
      j = find(raw_points_to_duplicate.begin(),raw_points_to_duplicate.end(),closest);
      // if closest raw point not already in vector
      if (j==raw_points_to_duplicate.end())
      {
        assert(closest<num_raw_points);
        raw_points_to_duplicate.push_back(closest);
      }
      // if left-hand side closest raw point not already in vector
      int lhs = closest-1;
      if (lhs<0) lhs=num_raw_points-1;
      j = find(raw_points_to_duplicate.begin(),raw_points_to_duplicate.end(),lhs);
      if (j==raw_points_to_duplicate.end())
      {
        assert(lhs<num_raw_points);
        raw_points_to_duplicate.push_back(lhs);
      }
      // if right-hand side closest raw point not already in vector
      int rhs = closest+1;
      if (rhs==num_raw_points) rhs=0;
      j = find(raw_points_to_duplicate.begin(),raw_points_to_duplicate.end(),rhs);
      if (j==raw_points_to_duplicate.end())
      {
        assert(rhs<num_raw_points);
        raw_points_to_duplicate.push_back(rhs);
      }
    }
  }
}

/** Calculate perpendicular distance between each spline samples
 *  and nearest linearly-interpolated raw contour points.
 * \param[in] raw_points All original points in contour.
 */

void Parameter::computeDeviations (const std::vector<Point> & raw_points)
{
  // for each spline sample
  int sample_index = 0;
  for (c_p_iterator k = spline_samples.begin();k!=spline_samples.end();k++)
  {
    closest_raw.push_back(-1);
    double min_d = 1E30;
    // for each line segment
    c_p_iterator i = raw_points.begin();
    c_p_iterator j = i;
    j++;
    int raw_point_index = 0;
    while (j!=raw_points.end())
    {
      double a = calcMinDistPtLineSegment(*k,*i,*j); 
      if (a < min_d)
      {
        min_d = a;
        // identify closest raw point
        double b = (*k).distance_squared(*i);
        double c = (*k).distance_squared(*j);
        if (b<c) closest_raw[sample_index]=raw_point_index;
        else     closest_raw[sample_index]=raw_point_index+1;
      }
      i++;
      j++;
      raw_point_index++;
    }
    i = raw_points.end();
    i--;
    j = raw_points.begin();
    double a = calcMinDistPtLineSegment(*k,*i,*j); 
    if (a < min_d)
    {
      min_d = a;
      // identify closest raw point
      double b = (*k).distance_squared(*i);
      double c = (*k).distance_squared(*j);
      if (b<c) closest_raw[sample_index]=raw_point_index;
      else     closest_raw[sample_index]=0;
    }
    deviations.push_back(min_d);
    sample_index++;
  }
}
