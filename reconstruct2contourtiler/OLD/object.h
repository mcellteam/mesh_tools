// Author: Justin Kinney
// Date: Feb 2009

#ifndef OBJECT_H
#define OBJECT_H 1

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "contour.h"
#include "controls.h"
#include "histogram.h"
#include "reconstruct2contourtiler.h"

using std::cout;
using std::endl;

typedef std::list<Contour>::iterator c_l_iterator;
typedef std::list<Contour>::const_iterator c_c_l_iterator;

class Object
{
private:
  std::string name;
  int min_section;
  int max_section;
  std::list<Contour> contours;
public:
  Object (char *str,int sec);
  Object &  operator =  (Object const &);
  Object      (Object const &);
  int getNumRawPoints (void);
  void printConfigFile (void);
  void appendScriptFile (void);
  void printPtsFiles (void);
  void addPreviousRaw (void);
  void interpolateRawPoints (void);
  void purgeBadContours (void);
  void removeDuplicates (void);
  bool degenerateObject (void) const;
  void addContour (Contour mycontour);
  void printCaps (void) const;
  void printRawPtsFiles (void);
  void computeHistogram (Histogram & h);
  void fitSplines (Histogram & h);
  void countPtsPerContour (std::vector<int> & points_per_contour);
  c_l_iterator getLastContour (void)
  {
    return --contours.end();
  }
  c_l_iterator firstContourNonConst (void)
  {
    return contours.begin();
  }
  c_l_iterator onePastLastContourNonConst (void)
  {
    return contours.end();
  }
  c_c_l_iterator firstContour (void) const
  {
    return contours.begin();
  }
  c_c_l_iterator onePastLastContour (void) const
  {
    return contours.end();
  }
  void addContourByRef (Contour & c)
  {
    contours.push_back(c);
  }
  void setMinSection (int i)
  {
    min_section = i;
  }
  void setMaxSection (int i)
  {
    max_section = i;
  }
  int getMaxSection (void) const
  {
    return max_section;
  }
  int getMinSection (void) const
  {
    return min_section;
  }
  char const * getName (void) const
  {
    return name.c_str();
  }
};

#endif
