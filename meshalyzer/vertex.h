// Author: Justin Kinney
// Date: Feb 2009

#ifndef VERTEX_H
#define VERTEX_H 1

#include <vector>

#include "face.h"

class Vertex
{
private:
  Vertex &  operator =            (Vertex const &);
  Vertex                          (Vertex const &);
public:
  int index;
  double pN[3];		  // current position coordinates (x,y,z)
  double pC[3];		  // closest mesh position coordinates (x,y,z)
  Face *cl;		  // pointer to face on which closest mesh position lies
  Object *o;		  // pointer to parent object
  std::vector<Edge*> e;	  // pointers to adjacent edges
  std::vector<Face*> f;	  // pointers to adjacent faces
  std::vector<Face*> nf;  // pointers to neighborhood faces
  Vertex(char* triplet,Object*);
  void getAdjacentVertices(std::vector<Vertex*>&);
  void getAdjacentFaces(hset_f&);
  void getNormal(double*);
  void printVertex(std::string);
  void printVertexCP(void);
  bool vertexIsNice(void);
  int getVertexNiceness(void);
  void setVertexNiceness(int);
  bool isManifold(bool);
  bool scanAdjFaces(Edge*,Face*,bool&);
};

#endif
