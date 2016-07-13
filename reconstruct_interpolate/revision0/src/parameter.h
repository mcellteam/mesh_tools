// Author: Justin Kinney
// Date: July 2009

#ifndef PARAMETER_H
#define PARAMETER_H 1

#include <cassert>
#include <vector>

#include "controls.h"
#include "point.h"

typedef std::vector<Point>                         vec_p;
typedef std::vector<Point>::iterator          p_iterator;
typedef std::vector<Point>::const_iterator  c_p_iterator;
typedef std::vector<double>::const_iterator c_d_iterator;
typedef std::vector<int>                           vec_i;

class Parameter
{
private:
  std::vector<double>                s; // path parameter value, 0<=s<1.0
  std::vector<int>        spline_index; // cubic curve segment index, 0<=index
  std::vector<double>            lhs_s; // path parameter distance to lhs sample
  std::vector<Point>    spline_samples; // spline evaluated at s
  std::vector<double>       deviations; // distance between spline_sample
                                        // and linearly-interpolated control point segment
  std::vector<int>         closest_raw; // index of closest raw point to spline sample
                                        // that forms line segment for deviation measurment
public:
  void print                      (char const * const filename) const;
  void computeDeviations          (const std::vector<Point> & raw_points);
  void getRawPointsToDuplicate    (const std::vector<Point> & raw_points,
                                   std::vector<int> & raw_points_to_duplicate) const;
  double calcMinDistPtLineSegment (const Point & spline_sample,
                                   const Point & contour_point_1,
                                   const Point & contour_point_3) const;
public:
  void addDeviation (double d)
  {
    deviations.push_back(d);
  }
  Parameter (void)
    :s(),spline_index(),lhs_s(),spline_samples(),deviations(),closest_raw()
  {
  }

  /** Determine if any deviation of sample points exceeds threshold.
   */

  bool deviationsLarge (void)
  {
    double deviation_threshold = Controls::instance().getDeviationThreshold();
    // for each deviation
    for (c_d_iterator i = deviations.begin();i!=deviations.end();i++)
    {
      if (*i>deviation_threshold) return true;
    }
    return false;
  }

  /** Clear all class data and allocate memory.
   * \param[in] target_num_samples Amount of memory to reserve for future data.
   */

  void reserve (const int & target_num_samples)
  {
    s.clear();
    s.reserve(target_num_samples);
    spline_index.clear();
    spline_index.reserve(target_num_samples);
    lhs_s.clear();
    lhs_s.reserve(target_num_samples);
    spline_samples.clear();
    spline_samples.reserve(target_num_samples);
    deviations.clear();
    deviations.reserve(target_num_samples);
    closest_raw.clear();
    closest_raw.reserve(target_num_samples);
  }

  /** Get sample point.
   * \param[in] sample_index Index of sample point of interest.
   * \return Sample point.
   */

  Point const * getPoint (const int & sample_index) const
  {
    int a = static_cast<int>(spline_samples.size())-1;
    if (sample_index<0)
    {
      return &spline_samples[0];
    }
    else if (sample_index>a)
    {
      return &spline_samples[a];
    }
    else
    {
      return &spline_samples[sample_index];
    }
  }

  /** Clear samples and reallocate memory.
   */

  void clearAndReserveSamples (void)
  {
    spline_samples.clear();
    spline_samples.reserve(s.size());
  }

  /** Add sample point.
   * \param[in] path_parameter_value Path parameter value of sample point for spline.
   * \param[in] index_of_spline Index of spline on which sample point lies.
   * \param[in] delta_s_to_lhs Path parameter distance to nearest sample point to the left.
   */

  void add (const double & path_parameter_value,
            const int & index_of_spline,
            const double & delta_s_to_lhs)
  {
    s.push_back(path_parameter_value);
    spline_index.push_back(index_of_spline);
    lhs_s.push_back(delta_s_to_lhs);
  }

  /** Move spline sample point.
   * \param[in] new_s New path parameter value.
   * \param[in] new_spline_index Index of new spline on which point sits.
   * \param[in] sample_index Index of sample point being moved.
   * \param[in] s_disp Change in path parameter value from old value. Signed value.
   * \param[in] p New location of sample point.
   */

  void set (double new_s,int new_spline_index,int sample_index,double s_disp,Point p)
  {
    // now abs refers to lhs s parameter separation 
    assert(sample_index>=0);
    assert(sample_index<static_cast<int>(s.size()));
    s[sample_index] = new_s;
    spline_index[sample_index] = new_spline_index;
    int i0 = sample_index+1;
    if (i0==getNumSamples()) i0=0;
    assert(i0>=0);
    assert(i0<getNumSamples());
    double a = lhs_s[i0]-s_disp;
    assert(a>0.0);
    lhs_s[i0]=a;
    a = lhs_s[sample_index]+s_disp;
    assert(a>0.0);
    lhs_s[sample_index]=a;
    spline_samples[sample_index]=p;
  }
  int getNumSamplePoints (void) const
  {
    return spline_samples.size();
  }
  int getNumSamples (void) const
  {
    assert(s.size()==spline_index.size());
    return s.size();
  }
  double getParameterValue (int sample_index) const
  {
    assert(sample_index>=0);
    assert(sample_index<static_cast<int>(s.size()));
    return s[sample_index];
  }
  int getSplineIndex (int sample_index) const
  {
    assert(sample_index>=0);
    assert(sample_index<static_cast<int>(spline_index.size()));
    return spline_index[sample_index];
  }
  double getLhsS (int sample_index) const
  {
    assert(sample_index>=0);
    assert(sample_index<static_cast<int>(lhs_s.size()));
    assert(lhs_s[sample_index]>0.0);
    return lhs_s[sample_index];
  }
  void addSplineSample (Point p)
  {
    spline_samples.push_back(p);
  }
  c_d_iterator getFirstDeviationConst (void) const
  {
    return deviations.begin();
  }
  c_d_iterator getOnePastLastDeviationConst (void) const
  {
    return deviations.end();
  }
  c_p_iterator getFirstSplineSampleConst (void) const
  {
    return spline_samples.begin();
  }
  c_p_iterator getOnePastLastSplineSampleConst (void) const
  {
    return spline_samples.end();
  }
};

#endif
