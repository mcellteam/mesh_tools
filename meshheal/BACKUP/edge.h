// Author: Justin Kinney
// Date: Sep 2008

#ifndef EDGE_H
#define EDGE_H 1

#include "face.h"
#include "meshheal.h"
#include "vertex.h"

class Edge
{
public:
  int hashval;
  int  v1; // vertex index
  int  v2; // vertex index
  int c12; // count of times edge was traversed from va to vb
  int c21; // count of times edge was traversed from vb to va
  void_list *f;
  Edge(void);
  void addFace(Face*);
  ~Edge(void);
public:
  Edge              (Edge const &);
  Edge & operator = (Edge const &);
  bool isValid (void);
  double getEdgeLengthSq (Vertex ** vert_array, const double & epsilon) const;
};

#endif    
