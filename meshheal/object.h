// Author: Justin Kinney
// Date: Feb 2009

#ifndef OBJECT_H
#define OBJECT_H 1

#include "edge.h"
#include "face.h"
#include "vertex.h"
#include "meshheal.h"

class Object
{
private:
  list_v v;       // container of pointers to all vertices in object
  vec_f  f;       // container of pointers to all faces in object
  vec_e  e;       // container of pointers to all edges in object
  int num_digits; // number of digits in largest vertex index
  map_se hm;      // create map for finding edges
public:
  Object  (std::string);
  ~Object (void);
public:
  void projectVerts       (const vec_v & free_vertices);
  void gatherFreeVertices (vec_v & free_vertices) const;
  void rebuildEdges       (void);
  void removeVertex       (Vertex const * bad_vertex);
  void createEdges        (void);
  void fixFaces           (Vertex * bad_vertex,Vertex * good_vertex);
  void computeDistances   (Closest & distances,
                           const vec_v & free_vertices);
  void addFinalFace       (const vec_v & free_vertices);
private:
  void scanFile (const std::string & filename);
  void buildEdge (Face*,Vertex*,Vertex*);
  Edge* findEdge (Vertex*,Vertex*);
  void createEdge (Face*,Vertex*,Vertex*);
  void setNumDigits (void);
  std::string keyPair (int a,int b) const;
  void getClosestFreeVerts (Vertex * current_vertex,
                            const map_dv & projections,
                            Closest & distances);
  c_dv_iterator getIterator (const map_dv & projections,
                             Vertex const * current_vertex,
                             const double & current_projection) const;
public:
  void print (std::ostream & target)
  {
    for (vl_iterator i=v.begin();i!=v.end();i++) 
    {
      (*i)->print(target);
      target << "\n";
    }
    for (f_iterator i=f.begin();i!=f.end();i++) 
    {
      (*i)->print(target);
    }
    for (e_iterator i=e.begin();i!=e.end();i++) 
    {
      (*i)->print(target);
    }
  }
  void printMesh (std::ostream & target)
  {
    for (vl_iterator i=v.begin();i!=v.end();i++) 
    {
      (*i)->print4mesh(target);
    }
    for (f_iterator i=f.begin();i!=f.end();i++) 
    {
      (*i)->print4mesh(target);
    }
  }
private:
  void addFace (Face * ff)
  {
    f.push_back(ff);
  }
  void addVertex (Vertex * vv)
  {
    v.push_back(vv);
  }
};

#endif
