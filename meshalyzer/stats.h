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
  //	int n;
  double total;
  double min;
  double max;
  double sum;
  double sum2;
  int * count;
  double * bins;
  vec_d x;
public:
  Stats(void);
  ~Stats(void);
  double getMin (void)
  {
    return min;
  }
  double getMax (void)
  {
    return max;
  }
  void setMin (double i)
  {
    min = i;
  }
  void setMax (double i)
  {
    max = i;
  }
  void add2total (double i)
  {
    total += i;
  }
  void add2sum2 (double i)
  {
    sum2 += i;
  }
  void add2sum (double i)
  {
    sum += i;
  }
  void insertElementRange (c_d_iterator i,c_d_iterator j)
  {
    x.insert(x.end(),i,j);
  }
  void insertElement (double i)
  {
    x.push_back(i);
  }
  int getSize (void)
  {
    return x.size();
  }
  c_d_iterator first_element (void)
  {
    return x.begin();
  }
  c_d_iterator one_past_last_element (void)
  {
    return x.end();
  }
  double getSum (void)
  {
    return sum;
  }
  double mean(void);
  double median(void);
  double variance(void);
  void clear(void);
  void printStats(void);
  void printHistogram(void);
  void createHistogram(void);
  void printAspectRatioHistogram(void);
  void createAspectRatioHistogram(void);
  void printAdjacentFaceHistogram(void);
  void createAdjacentFaceHistogram(void);
};

#endif
