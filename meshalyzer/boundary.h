// Author: Justin Kinney
// Date: Feb 2009

#ifndef BOUNDARY_H
#define BOUNDARY_H 1

#include "meshalyzer.h"

class Boundary
{
private:
  Boundary & operator = (Boundary const &);
  Boundary              (Boundary const &);
  vec_e e;
  Vertex *end,*begin;
  bool open;
public:
  Boundary     (void);
  void init    (Edge*);
  bool add     (Edge*);
  void print   (void);
  bool closed  (void);
  void setOpen (bool b)
  {
    open=b;
  }
  bool isOpen  (void)
  {
    return open;
  }
  bool edgeExtendsBoundary (Edge*);
};

#endif
