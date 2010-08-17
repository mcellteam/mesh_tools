#ifndef VERTEX_H
#define VERTEX_H 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

class Vertex
{
public:
  double x,y,z;
  int index;
  Vertex(char *triplet);
};

#endif
