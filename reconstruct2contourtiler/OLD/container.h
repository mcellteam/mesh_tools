// Author: Justin Kinney
// Date: Feb 2009

#ifndef CONTAINER_H
#define CONTAINER_H 1

#include <vector>

#include "object.h"

typedef std::vector<Object>::iterator o_iterator;
typedef std::vector<Object>::const_iterator c_o_iterator;

class Container
{
private:
  std::vector<Object> o;
  std::vector<int> points_per_contour;
  int num_files_read;
public:
  Container (void);
  int getMatchingObject (char const *);
  int newObject (char * myname,int i);
  int getNumRawPoints (void);
  void clearPtsFiles (void);
  void computeHistogram (Histogram & h);
  void fitSplines (Histogram & h);
  void interpolateRawPoints (void);
  void purgeBadContours (void);
  void removeDuplicates (void);
  void addPreviousRaw (void);
  void getContours (void);
  void initTransform (double *t);
  void setTransform (double *t,char *p);
  void addContour2Object (char * name,int section);
  void countPtsPerContour (void);
  c_l_iterator getLastContourFromObject (char * name);
  o_iterator firstObjectNonConst (void)
  {
    return o.begin();
  }
  o_iterator onePastLastObjectNonConst (void)
  {
    return o.end();
  }
  c_o_iterator firstObject (void)
  {
    return o.begin();
  }
  c_o_iterator onePastLastObject (void)
  {
    return o.end();
  }
};

#endif
