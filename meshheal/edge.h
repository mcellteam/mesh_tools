// Author: Justin Kinney
// Date: June 2010

#ifndef EDGE_H
#define EDGE_H 1

#include "vertex.h"

class Edge
{
private:
  Edge & operator = (Edge const &);
  Edge              (Edge const &);
  Vertex *vv1,*vv2; // pointers to vertices on edge
  Face *f1,*f2;	    // pointers to adjacent faces (i.e. faces that contain edge)
public:
  Edge (Face*,Vertex*,Vertex*);
public:
  void getVertices (Vertex*&,Vertex*&,Vertex*&,Vertex*&) const;
  void addFace     (Face*);
  void print       (std::ostream&) const;
public:
  bool isBorder (void) const
  {
    assert(f1!=NULL);
    return f2==NULL;
  }
  Face * ptr_f1 (void) const
  {
    return f1;
  }
  Face * ptr_f2 (void) const
  {
    return f2;
  }
  Vertex * ptr_vv1 (void) const
  {
    return vv1;
  }
  Vertex * ptr_vv2 (void) const
  {
    return vv2;
  }
private:
  bool isMatch (Vertex const * va,
                Vertex const * vb) const
  {
    if (va==vv1 && vb==vv2) return true;
    if (va==vv2 && vb==vv1) return true;
    return false;
  }
};

#endif
