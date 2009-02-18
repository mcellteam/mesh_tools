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
public:
  Vertex *vv1,*vv2; // pointers to vertices on edge
  Face *f1,*f2;	    // pointers to adjacent faces (i.e. faces that contain edge)
  vec_f fvec;	    // pointers to additional adjacent faces
  double l;	    // original edge length
  Edge(Face*,Vertex*,Vertex*);
  void update(Face*);
  double getSqLength(void);
  double getAngle(void);
  void printEdge(std::string);
  void printEdgeCP();
  bool isConsistent(void);
  bool isManifold(void);
  bool getStartingFace(Face*&);
  Face* getNewFace(Face*);
  void getVertices(Vertex*&,Vertex*&,Vertex*&,Vertex*&);
  bool valid(void);
};

#endif
