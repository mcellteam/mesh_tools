// Author: Justin Kinney
// Date: Feb 2009

#ifndef CONTROL_POINTS_H
#define CONTROL_POINTS_H 1

#include <cassert>
#include <list>

#include "controls.h"

class Contour;

class Control_Points
{
private:
  std::vector<int> index;
public:
  Control_Points      (void);
private:
  void loadMatrix     (std::list<int> & closed_path);
  bool checkDeviation (Contour const *,int num_raw_points);
public:
  /** Calculate number of splines defined by control points.
   * \return Number of splines.
   */

  int getNumSplines (void) const
  {
    int i = static_cast<int>(index.size());
    assert(!(i%4));
    return i/4;
  }

  /** Create matrix of control points describing all splines in contour.
   * \param[in] control_points Sequence of raw points to use as control_points.
   */

  void initialize (std::list<int> & control_points)
  {
    int num_points = control_points.size();
    // at least enough points for a single spline
    assert(num_points>3);
    index.clear();
    index.reserve(4*num_points);
    loadMatrix(control_points);
  }

  /** Get control point.
   * \param[in] i Index of control point of interest.
   * \return Index of raw contour point stored as control point.
   */

  int getControlPointIndex (int i) const
  {
#ifndef NOASSERT
    assert(i>=0); 
    assert(i<static_cast<int>(index.size())); 
#endif
    return index[i];
  }

  /** Write control points (e.g. raw contour point indices) to file.
   * \param[in] filename Output file name.
   */

  void print (char * filename) const
  {
    // assert numbert indeces is multiple of 4
    int num_indeces = static_cast<int>(index.size());
    assert(!(num_indeces%4));
    FILE *F = fopen(filename,"a");
    if (!F) { printf("Could not open output file %s\n",filename); return; }
    char line[2048];
    sprintf(line,"\n# num splines = %i\n",num_indeces/4);
    fputs(line,F);
    for (int i=0;i<num_indeces;i+=4)
    {
      assert((i+3)<static_cast<int>(index.size()));
      sprintf(line,"%i %i %i %i\n",index[i+0],index[i+1],index[i+2],index[i+3]);
      fputs(line,F);
    }
    fclose(F);
  }
};

#endif
