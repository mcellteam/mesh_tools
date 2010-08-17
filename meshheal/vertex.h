// Author: Justin Kinney
// Date: June 2010

#ifndef VERTEX_H
#define VERTEX_H 1

#include <cassert>
#include <iostream>

#include "controls.h"

class Vertex
{
private:
  Vertex &  operator =            (Vertex const &);
  Vertex                          (Vertex const &);
  int index;
  vector3 p;
  double projection; // projection of vertex onto random vector
public:
  Vertex(char* triplet);
  double get3Ddistance (Vertex const * v) const
  {
    const vector3 diff = p - v->getpN_ptr();
    return diff.length();
  }
  void project (void)
  {
    projection = p.dot(Controls::instance().get_random_vector());
  }
  double getProjection (void) const
  {
    return projection;
  }
  double getpN (int i) const
  {
    assert(i<3);
    return p.p[i];
  }
  vector3 getpN_ptr () const
  {
    return p;
  }
  int getIndex (void) const
  {
    return index;
  }
  void print (std::ostream & target) const
  {
    target.precision(12);
    target << "Vertex "
          << " <ind> " << index << " "
          << " ["
          << p.p[0] << " "
          << p.p[1] << " "
          << p.p[2] << "]";
  }
  void print4mesh (std::ostream & target) const
  {
    target.precision(12);
    target << "Vertex "
          << index << " "
          << p.p[0] << " "
          << p.p[1] << " "
          << p.p[2] << std::endl;
  }
};

struct my_ltv
{
  bool operator () (const Vertex * s1, const Vertex * s2) const
  {
    return s1->getIndex()<s2->getIndex();
  }
};

#endif
