// Author: Justin Kinney
// Date: Feb 2009

#ifndef CONTOUR_H
#define CONTOUR_H 1

#include <string>
#include <vector>
#include <list>

#include "histogram.h"
#include "point.h"
#include "reconstruct2contourtiler.h"

typedef std::vector<double>::iterator d_iterator;
typedef std::vector<double>::const_iterator c_d_iterator;
typedef std::vector<Point>::iterator p_iterator;
typedef std::vector<Point>::const_iterator c_p_iterator;
typedef std::list<Point>::iterator p_l_iterator;
typedef std::list<Point>::const_iterator c_p_l_iterator;

class Contour
{
private:
  int  section;
  std::string         name;
  std::list<Point>    raw_points;
  std::vector<Point>  spline_samples;
  std::vector<double> deviations;
public:
  Contour (char* str, int sec);
  Contour &  operator =  (Contour const &);
  Contour                (Contour const &);
  void printDiagnostics (std::vector<SplinePoint> & sp,int num,char const * const outdir,int tag);
  void computeHistogram (Histogram & h);
  void computeDeviation (std::vector<SplinePoint> & sp,const int count,const int num);
  void sampleSplines (std::vector<SplinePoint> & sp,double dt,int limit,int tag);
  void fitSplines (Histogram & h);
  void removeDuplicates (void);
  void linearlyInterp (void);
  void addPreviousRaw (void);
  void addPreviousInterp (void);
  void interpPoints (p_l_iterator & lhs,p_l_iterator & rhs);
  void allocateDeviations (int count)
  {
    deviations.reserve(count);
    for (int j=0;j<count;j++)
    {
      deviations[j]=0;
    }
  }
  void addRawPoint (Point p)
  {
    raw_points.push_back(p);
  }
  int getNumRawPoints (void) const
  {
    return raw_points.size();
  }
  bool pointsAreIdentical (Point const & lhs,Point const & rhs)
  {
    return lhs.getX()==rhs.getX() && lhs.getY()==rhs.getY();
  }
  void clearSpline (void)
  {
    spline_samples.clear();
  }
  char const * getName (void) const
  {
    return name.c_str();
  }
  c_p_iterator getFirstSplineSample (void)
  {
    return spline_samples.begin();
  }
  c_p_iterator getOnePastLastSplineSample (void)
  {
    return spline_samples.end();
  }
  c_p_l_iterator getFirstRawPoint (void)
  {
    return raw_points.begin();
  }
  c_p_l_iterator getOnePastLastRawPoint (void)
  {
    return raw_points.end();
  }
  int getSection (void) const
  {
    return section;
  }
  int getNumSplineSamples (void) const
  {
    return spline_samples.size();
  }
  double getDeviation (int i) const
  {
    //assert(i<static_cast<int>(deviations.size()));
    return deviations[i];
  }
  void setDeviation (int i,double d)
  {
    deviations[i] = d;
  }
};

#endif
