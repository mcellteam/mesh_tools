// Author: Justin Kinney
// Date: Feb 2009

#ifndef HISTOGRAM_H
#define HISTOGRAM_H 1

#include "reconstruct2contourtiler.h"

class Histogram
{
private:
  double min;
  double max;
  double mean;
  double stddev;
  double sum;
  int count[16];
  int num;
public:
  Histogram(void);
  void update (double const &);
  void printStatistics (void);
  void incrementCount (int i)
  {
    count[i]++;
  }
  void setMin (double d)
  {
    min = d;
  }
  void setMax (double d)
  {
    max = d;
  }
  void setMean (double d)
  {
    mean = d;
  }
  void setStddev (double d)
  {
    stddev = d;
  }
  void setSum (double d)
  {
    sum = d;
  }
  void setNum (int i)
  {
    num = i;
  }
  double getMin (void)
  {
    return min;
  }
  double getMax (void)
  {
    return max;
  }
  double getMean (void)
  {
    return mean;
  }
  double getStddev (void)
  {
    return stddev;
  }
  double getSum (void)
  {
    return sum;
  }
  int getNum (void)
  {
    return num;
  }
};

#endif
