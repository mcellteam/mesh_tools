// Author: Justin Kinney
// Date: Feb 2009

#ifndef EDGE_H
#define EDGE_H 1

#include "vertex.h"

class Edge
{
private:
  Edge &  operator =            (Edge const &);
  Edge                          (Edge const &);
  Vertex *vv1,*vv2; // pointers to vertices on edge
  Face *f1,*f2;	    // pointers to adjacent faces (i.e. faces that contain edge)
  vec_f fvec;	    // pointers to additional adjacent faces
  double l;	    // original edge length
public:
  Edge (Face*,Vertex*,Vertex*);
  bool isConsistent(void);
  bool isManifold(void);
  bool getStartingFace(Face*&);
  bool valid(void);
  void getVertices(Vertex*&,Vertex*&,Vertex*&,Vertex*&);
  void update(Face*);
  void print   (std::ostream&);
  void printCP (std::ostream&);
  Face* getNewFace(Face*);
  double getSqLength(void);
  double getAngle(void);
  double getOrigLength (void)
  {
    return l;
  }
  c_f_iterator first_extra_face (void)
  {
    return fvec.begin();
  }
  c_f_iterator one_past_last_extra_face (void)
  {
    return fvec.end();
  }
  bool noExtraFaces (void)
  {
    return fvec.empty();
  }
  Face * ptr_f1 (void)
  {
    return f1;
  }
  Face * ptr_f2 (void)
  {
    return f2;
  }
  Vertex * ptr_vv1 (void)
  {
    return vv1;
  }
  Vertex * ptr_vv2 (void)
  {
    return vv2;
  }
};

#endif
