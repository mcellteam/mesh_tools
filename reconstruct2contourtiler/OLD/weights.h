// Author: Justin Kinney
// Date: Feb 2009

#ifndef WEIGHTS_H
#define WEIGHTS_H 1

#include "reconstruct2contourtiler.h"

typedef std::vector<SplinePoint> v_sp;

class Weights
{
private:
  std::vector<double> bx;
  std::vector<double> by;
public:
  Weights (void);
  void loadWeightArrays (p_l_iterator,p_l_iterator,int const & limit);
  void computeSplines (v_sp & sp,double dt,int limit);
  bool checkDeviation (Contour const *,int limit);
  void freeX (void)
  {
    bx.clear();
  }
  void freeY (void)
  {
    by.clear();
  }
  void allocateX (int i)
  {
    bx.reserve(i);
  }
  void allocateY (int i)
  {
    by.reserve(i);
  }
};

#endif
