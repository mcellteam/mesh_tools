// Author: Justin Kinney
// Date: Feb 2009

#ifndef STATS_H
#define STATS_H 1

#include "controls.h"
#include "meshalyzer.h"

class Stats
{
private:
  Stats &  operator =            (Stats const &);
  Stats                          (Stats const &);
public:
  //	int n;
  double total;
  double min;
  double max;
  double sum;
  double sum2;
  int * count;
  double * bins;
  vec_d x;
  void clear(void);
  Stats(void);
  ~Stats(void);
  double mean(void);
  double median(void);
  double variance(void);
  void createHistogram(void);
  void createAdjacentFaceHistogram(void);
  void createAspectRatioHistogram(void);
  void printHistogram(void);
  void printAspectRatioHistogram(void);
  void printAdjacentFaceHistogram(void);
  void printStats(void);
};

#endif
