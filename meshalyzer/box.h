// Author: Justin Kinney
// Date: Feb 2009

#ifndef BOX_H
#define BOX_H 1

#include <vector>

#include "face.h"

//typedef std::vector<Face*>::iterator        vec_f_it;
//typedef std::vector<Face*>::const_iterator  vec_f_cit;

class Box
{
private:
  std::vector<Face*> f;
  int x,y,z;  // indices of box in space
public:
  void getFaceIntersection(Container*);
  bool faceIntersectionAlreadyKnown(Face*,Face*);
  void printBox(Space*);
  double xmin (double i,double sl) { return x*sl+i;     }
  double xmax (double i,double sl) { return (x+1)*sl+i; }
  double ymin (double i,double sl) { return y*sl+i;     }
  double ymax (double i,double sl) { return (y+1)*sl+i; }
  double zmin (double i,double sl) { return z*sl+i;     }
  double zmax (double i,double sl) { return (z+1)*sl+i; }

  Box (int a,int b,int c)
    :f(),x(a),y(b),z(c)
  {
  }
  c_f_iterator first_face (void)
  {
    return f.begin();
  }
  c_f_iterator one_past_last_face (void)
  {
    return f.end();
  }
  bool noFaces (void)
  {
    return f.empty();
  }
  void clearFaces (void)
  {
    f.clear();
  }
  void addFace (Face * ff)
  {
    f.push_back(ff);
  }
  int numFaces (void)
  {
    return f.size();
  }
};

#endif
