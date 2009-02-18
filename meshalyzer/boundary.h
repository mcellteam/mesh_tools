// Author: Justin Kinney
// Date: Feb 2009

#ifndef BOUNDARY_H
#define BOUNDARY_H 1

#include "meshalyzer.h"

class Boundary
{
private:
  Boundary &  operator =            (Boundary const &);
  Boundary                          (Boundary const &);
public:
  vec_e e;
  Vertex *end,*begin;
  bool open;
  Boundary (void);
  void init (Edge*);
  bool edgeExtendsBoundary (Edge*);
  bool add (Edge*);
  void print (void);
  bool closed (void);
};

#endif
