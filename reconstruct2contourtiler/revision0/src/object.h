// Author: Justin Kinney
// Date: Feb 2009

#ifndef OBJECT_H
#define OBJECT_H 1

#include <cassert>
#include <iostream>
#include <list>
#include <stdlib.h>
#include <string>
#include <vector>

#include "contour.h"
#include "controls.h"
#include "histogram.h"

using std::cout;
using std::endl;

typedef std::list<Contour>                 c_list;
typedef std::list<Contour>::iterator       c_l_iterator;
typedef std::list<Contour>::const_iterator c_c_l_iterator;

class Object
{
private:
  std::string name;
  int  min_section;
  int  max_section;
  c_list  contours;
public:
  Object (char const * const str,const int & sec);
  Object & operator = (Object const &);
  Object (Object const &);
  int  getNumRawPoints     (void) const;
  int  getRanges           (std::vector<int> & ranges) const;
  void printRawPoints      (void) const;
  void purgeBadContours    (void);
  void removeDuplicates    (void);
  void processContour      (Histogram & h,Histogram & si,Histogram & si_before);
  void addContour          (Contour mycontour);
  void computeHistogram    (Histogram & h,Histogram & si);
  void writeOutputContours (const int & num_parts,const std::vector<int> & ranges);
private:
  bool isDegenerate        (void) const;
  void printConfigFile     (char  const * myname,
                            const int & min_sec,
                            const int & max_sec) const;
  void appendScriptFile    (char const * myname) const;
  void printPtsFiles       (char const * myname,
                            const int & min_sec,
                            const int & max_sec) const;
  void printCaps           (char const * myname,
                            const int & min_sec,
                            const int & max_sec) const;
  void printRawPtsFiles    (char const * myname,
                            const int & min_sec,
                            const int & max_sec) const;
public:
  /** Initialize all output pts files to empty.
   */
  
  void clearPtsFiles (const int & num_parts,
                      const std::vector<int> & ranges) const
  {
    Controls & cs(Controls::instance());
    char filename[256];
    if (num_parts==1)
    {
      for (int j = min_section;j<=max_section;j++)
      {
        sprintf(filename,"%s%s%d.pts",cs.getOutputDir(),name.c_str(),j);
        FILE * F = fopen(filename,"w");
        if (!F) { printf("Couldn't open output file %s\n",filename); exit(0); }
        assert(F);
        fclose(F);
      }
    }
    else
    {
      char temp[256];
      for (int i=0;i<num_parts;i++)
      {
        assert((2*i+1)<static_cast<int>(ranges.size()));
        for (int j=ranges[2*i];j<=ranges[2*i+1];j++)
        {
          sprintf(temp,"%s%s%s%d.pts",cs.getOutputDir(),name.c_str(),
                                          cs.getMultiPartSuffix(),j);
          sprintf(filename,temp,i);
          FILE * F = fopen(filename,"w");
          if (!F) { printf("Couldn't open output file %s\n",filename); exit(0); }
          assert(F);
          fclose(F);
        }
      }
    }
  }

  int getNumContours (void) const
  {
    return contours.size();
  }
  c_c_l_iterator firstContour (void) const
  {
    return contours.begin();
  }
  c_c_l_iterator onePastLastContour (void) const
  {
    return contours.end();
  }
  c_l_iterator getLastContour (void)
  {
    return --contours.end();
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
