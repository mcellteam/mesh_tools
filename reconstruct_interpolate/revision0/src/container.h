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
  std::vector<Object> o; // separate objects in input files
  int    num_files_read; // number of input files
private:
  int  getMatchingObject   (char const * const) const;
  bool parseTransform      (FILE * stream,double * const transform);
  bool parseContour        (char const * line,
                            const int & section,char * const name);
  void setTransform        (char const * str,double * const transform);
  void addContour2Object   (char * const name,const int & section);
  void createCallingScript (char const * const outdir,char const * script);
public:
  void processContour (Histogram & h,Histogram & si,Histogram & si_before);
  void clearOutputScripts  (void);
  void getContours         (void);
  void writeOutputContours (void);
private:
  /** Create a new object.
   * \param[in] object_name Name of new object.
   * \param[in] section Section number of new object.
   * \return Index of new object.
   */

  int newObject (char const * object_name,const int & section)
  {
    assert(o.size()<o.capacity());
    o.push_back(Object(object_name,section));
    assert(o.size()>0);
    return o.size()-1;
  }

  /** Initialize matrix for transforming contour data.
   * \param[out] t Matrix of transform coefficients.
   */

  void initTransform (double * const t)
  {
    t[0]=0.0; t[1]=1.0; t[2]=0.0; t[3]=0.0; t[4]=0.0; t[5]=0.0;
    t[6]=0.0; t[7]=0.0; t[8]=1.0; t[9]=0.0; t[10]=0.0; t[11]=0.0;
  }

public:
  /** Create new instance of Container class. 
   */

  Container (void)
  :o(),num_files_read(0)
  {
    o.reserve(Controls::instance().getNumReserve());
  }

  /** Get number of input files read.
   * \return Number of input files read.
   */

  int getNumFiles (void) const
  {
    return num_files_read;
  }

  /** Get number of different input objects found.
   * \return Number of objects.
   */

  int getNumObjects (void) const
  {
    return o.size();
  }

  /** Get cumulative number of contours of all objects. 
   */

  int getNumContours (void) const
  {
    int sum = 0;
    for (c_o_iterator i = o.begin();i!=o.end();i++)
    {
      sum += i->getNumContours();
    }
    return sum;
  }

  /** Get cumulative number of points in all contours of all objects. 
   */

  int getNumRawPoints (void) const
  {
    int sum = 0;
    for (c_o_iterator i = o.begin();i!=o.end();i++)
    {
      sum += i->getNumRawPoints();
    }
    return sum;
  }

  /** Get iterator to last contour of object with specified name. 
   * \param[in] object_name Object of interest.
   * \return Iterator to last contour in object.
   */

  c_l_iterator getLastContourFromObject (char const * const object_name)
  {
    int index = getMatchingObject(object_name);
    assert(index>=0);
    assert(index<static_cast<int>(o.size()));
    return o[index].getLastContour();
  }

  /** Get iterator to first object. 
   */

  o_iterator firstObjectNonConst (void)
  {
    return o.begin();
  }

  /** Get iterator to one past last object. 
   */

  o_iterator onePastLastObjectNonConst (void)
  {
    return o.end();
  }

  /** Get iterator to first object. 
   */

  c_o_iterator firstObject (void) const
  {
    return o.begin();
  }

  /** Get iterator to one past last object. 
   */

  c_o_iterator onePastLastObject (void) const
  {
    return o.end();
  }

  /** Remove contours with too few points in each object in container.
   */

  void purgeBadContours (void)
  {
    // for each object
    for (o_iterator i=o.begin();i!=o.end();i++)
    {
      i->purgeBadContours();
    }
  }

  /** Remove duplicate points in each object in container.
   */

  void removeDuplicates (void)
  {
    // for each object
    for (o_iterator i=o.begin();i!=o.end();i++)
    {
      i->removeDuplicates();
    }
  }

  /** For each object in container compute histogram of deviaton distance
   *  between spline samples and linearly interpolated input contour points.
   * \param[in] h Instance of Histogram class for sample point deviations.
   * \param[in] si Instance of Histogram class for sample intervals
   *  after simulated annealing.
   */

  void computeHistogram (Histogram & h,Histogram & si)
  {
    for (o_iterator i = o.begin();i!=o.end();i++)
    {
      i->computeHistogram(h,si);
    }
  }

  /** For each object in container write input contour points to file.
   */

  void printRawPoints (void) const
  {
    for (c_o_iterator i = o.begin();i!=o.end();i++)
    {
      i->printRawPoints();
    }
  }
};

#endif
