// Author: Justin Kinney
// Date: Feb 2009

#ifndef FACE_H
#define FACE_H 1

#include <vector>

#include "meshalyzer.h"

class Face
{
private:
  Face &  operator =  (Face const &);
  Face                (Face const &);
  int     index;       // face index
  int     vi[3];       // vertex indices
  Vertex * v[3];       // pointers to vertices
  Edge   * e[3];       // pointers to edges
  vec_b       b;       // pointers to boxes
public:
  Face (char*,Object*); 
  Face (Vertex*,Vertex*,Vertex*); 
  bool faceInTable_intf(void);
  bool match(int);
  void addEdge(Edge*);
  void recordBoxes(std::vector<Box*>&); 
  void getNormal(double[3]);
  void getVertexCoordinates(double *[3]);
  void print(std::ostream &);
  void printCP(std::ostream &);
  void addFaceToTable_intf(void);
  void addFaceToVector(Face*);
  void clearFaceFromTable_intf(void);
  Edge* getNewEdge(Edge*,Vertex*);
  double getAspectRatio(void);
  double getAngle(Vertex *v);
  ff_iterator findFaceInTable_intf(void);
  Edge * ptr_edge (int i)
  {
    return e[i];
  }
  Vertex * ptr_vertex (int i)
  {
    return v[i];
  }
  int getIndex (void)
  {
    return index;
  }
};

#endif
