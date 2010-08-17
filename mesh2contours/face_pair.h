// Author: Justin Kinney
// Date: Feb 2009

#ifndef FACE_PAIR_H
#define FACE_PAIR_H 1

#include "vertex.h"

class Face_Pair
{
private:
  Vertex *a,*b,*c,*d;  // vertex*s
  Edge     * e;        // shared edge
  Face *f1,*f2;        // bad faces
  Object   * o;        // parent object of faces
  int   next_i;        // next available face index
public:
  Face_Pair (Object*);
  bool aspectRatiosImprove(void);
  bool existingEdge(void);
  void clear(void);
  void processBadFace(Face*);
  void print(void);
  void analyzeF1(void);
  void findF2(void);
};

#endif
