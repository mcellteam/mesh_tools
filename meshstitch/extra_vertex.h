#ifndef EXTRA_VERTEX_H
#define EXTRA_VERTEX_H 1

#include "vertex.h"

class ExtraVertex
{
public:
  int g;	// group number
  Vertex *v;
  Vertex *la,*lb;
  ExtraVertex(void);
};

#endif
