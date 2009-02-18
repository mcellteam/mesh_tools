// Author: Justin Kinney
// Date: Feb 2009

#ifndef FACE_H
#define FACE_H 1

#include <vector>

#include "meshalyzer.h"

class Face
{
private:
  Face &  operator =            (Face const &);
  Face                          (Face const &);
public:
  int index;	// face index
  int vi[3];	// vertex indices
  Vertex *v[3];	// pointers to vertices
  Edge   *e[3];	// pointers to edges
  vec_b b;	// pointers to boxes
  Face(char*,Object*); 
  Face(Vertex*,Vertex*,Vertex*); 
  void addEdge(Edge*);
  void recordBoxes(std::vector<Box*>&); 
  void getNormal(double[3]);
  void getVertexCoordinates(double *[3]);
  double getAngle(Vertex *v);
  void printFace(Object*);
  void printFaceCP(void);
  void addFaceToTable_intf(void);
  bool faceInTable_intf(void);
  ff_iterator findFaceInTable_intf(void);
  void addFaceToVector(Face*);
  void clearFaceFromTable_intf(void);
  Edge* getNewEdge(Edge*,Vertex*);
  bool match(int);
  double getAspectRatio(void);
};

#endif
