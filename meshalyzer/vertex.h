// Author: Justin Kinney
// Date: Feb 2009

#ifndef VERTEX_H
#define VERTEX_H 1

#include <cassert>
#include <iostream>
#include <vector>

#include "face.h"
#include "object.h"

class Vertex
{
private:
  Vertex &  operator =            (Vertex const &);
  Vertex                          (Vertex const &);
  int index;
  double pN[3];		  // current position coordinates (x,y,z)
  double pC[3];		  // closest mesh position coordinates (x,y,z)
  Face *cl;		  // pointer to face on which closest mesh position lies
  Object *o;		  // pointer to parent object
public:
  Vertex(char* triplet,Object*);
  Object * getObject (void)
  {
    return o;
  }
  void setpC (int i,double j)
  {
    pC[i] = j;
  }
  double getpC (int i)
  {
    assert(i<3);
    return pC[i];
  }
  double getpN (int i)
  {
    assert(i<3);
    return pN[i];
  }
  double * getpN_ptr ()
  {
    return pN;
  }
  void setClosestFacePtr (Face * ptr)
  {
    cl = ptr;
  }
  Face * ptr_closest_face (void)
  {
    return cl;
  }
  int getIndex (void)
  {
    return index;
  }
  void print (std::ostream & target)
  {
    target.precision(12);
    target << "Vertex <obj> " << o->getName()
          << " <ind> " << index << " "
          << " ["
          << pN[0] << " "
          << pN[1] << " "
          << pN[2] << "]";
  }

  void printCP (std::ostream & target)
  {
    target.precision(12);
    target << pN[0] << " "
          << pN[1] << " "
          << pN[2] << " 1 0 0 1\n";
  }

  std::vector<Edge*> e;	  // pointers to adjacent edges
  std::vector<Face*> f;	  // pointers to adjacent faces
  std::vector<Face*> nf;  // pointers to neighborhood faces
  int  getVertexNiceness(void);
  void getAdjacentVertices(std::vector<Vertex*>&);
  void getAdjacentFaces(hset_f&);
  void getNormal(double*);
//  void printVertex(std::string);
  void printVertexCP(void);
  bool vertexIsNice(void);
  void setVertexNiceness(int);
  bool isManifold(bool);
  bool scanAdjFaces(Edge*,Face*,bool&);
};

#endif
