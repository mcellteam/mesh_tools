#ifndef FACE_H
#define FACE_H 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vertex.h"

class Face
{
public:
  int index;	// Face index
  int v1,v2,v3;	// vertex indices
  Face (char const * triplet);
  int vertexInFace (Vertex const * const V)
  {
    return ((v1==V->index)||(v2==V->index)||(v3==V->index));
  }

};

#endif
