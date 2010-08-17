// Author: Justin Kinney
// Date: Apr 2010

#ifndef CONTOUR_H
#define CONTOUR_H 1

#include <vector>

#include "controls.h"

#include "vertex.h"

class Contour
{
public:
  //vec_v p;			// Vertex* points along contour
  std::vector<Vertex*> p;			// Vertex* points along contour
  //void print (Controls &cs,std::string);
  void print (std::string);
  Contour (void);
};


#endif


