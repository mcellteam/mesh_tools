// Author: Justin Kinney
// Date: Sep 2008

#include "edge.h"
#include "void_list.h"

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>

Edge::Edge(void)
:hashval(-1),v1(0),v2(0),c12(0),c21(0),f(NULL)
{
}


Edge::Edge (const Edge& rhs)           
:hashval(rhs.hashval),v1(rhs.v1),v2(rhs.v2),c12(rhs.c12),c21(rhs.c21),f(rhs.f)
{                                                        
}                                             

Edge& Edge::operator = (const Edge& rhs)
{                                                                              
  std::cout << "Assignment operator prohibited on instances of Edge class.\n";    
  std::cout << "Edge " << rhs.v1 << std::endl;                                      
  exit(1);                                               
}                                                        

Edge::~Edge (void)
{
  void_list *q,*p;
  q=f;
  while (q!=NULL)
  {
    p=q->next;
    delete q;
    q=p;
  }
}

void Edge::addFace (Face *fc)
{
  //add face to edge
  void_list *qq;
  qq = new void_list();
  qq->next = f;
  qq->data = (void*)fc;
  f = qq;
}

bool Edge::isValid (void)
{
  return v1!=0 && v2!=0;
}

double Edge::getEdgeLengthSq (Vertex ** vert_array,
                              const double & epsilon) const
{
  const double scale = 1.0/epsilon;
  Vertex * va = vert_array[v1];
  Vertex * vb = vert_array[v2];
  const double ax = va->x*scale;
  const double bx = vb->x*scale;
  const double ay = va->y*scale;
  const double by = vb->y*scale;
  const double az = va->z*scale;
  const double bz = vb->z*scale;
  return (ax-bx)*(ax-bx)+(ay-by)*(ay-by)+(az-bz)*(az-bz);
}


