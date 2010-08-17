// Author: Justin Kinney
// Date: June 2010

#ifndef FACE_H
#define FACE_H 1

#include "meshheal.h"

class Face
{
private:
  Face &  operator =  (Face const &);
  Face                (Face const &);
  int     index; // face index
  Vertex * v[3]; // pointers to vertices
  Edge   * e[3]; // pointers to edges
public:
  Face (char*,mmap_iv & vp); 
  Face (int i, Vertex * v1, Vertex *v2, Vertex *v3);
public:
  void addEdge          (Edge*);
  void print            (std::ostream &) const;
  void print4mesh       (std::ostream &) const;
public:
  Edge * ptr_edge (int i) const
  {
    return e[i];
  }
  Vertex * ptr_vertex (int i) const
  {
    return v[i];
  }
  int getIndex (void) const
  {
    return index;
  }
  void replaceVertices (Vertex * bad_vertex,
                        Vertex * good_vertex)
  {
    for (int i=0;i<3;++i)
    {
      if (v[i]==bad_vertex) v[i]=good_vertex;
    }
  }
  void clearEdges (void)
  {
    for (int i=0;i<3;++i) e[i]=NULL;
  }
  bool vertSeqFound (Vertex *a,Vertex *b)
  {
    if (a==v[0] && b==v[1]) return true;
    if (a==v[1] && b==v[2]) return true;
    if (a==v[2] && b==v[0]) return true;
    return false;
  }
};

#endif
