// Author: Justin Kinney
// Date: Sep 2008

#ifndef EDGE_H
#define EDGE_H 1

#include "meshmorph.h"

#include "vertex.h"

class Edge
{
private:
  Face   *f1,*f2; // pointers to adjacent faces
  float        l; // original edge length
  Vertex *v1,*v2; // Vertices at each end of edge
  Vertex *o1,*o2; // vertices of adjacent faces not part of this edge
public:
  Edge                          (Face * const ,
                                 Vertex * const ,
                                 Vertex * const ,
                                 Vertex * const);
  void    print                 (std::ostream &) const;
  void    update                (Face * const,Vertex * const);
  double  getAngle              (void) const;
  double  getCurvatureLength    (Vertex const * const o,
                                 Vertex const * const v1,
                                 Vertex const * const v2,
                                 double const & edge_length) const;
  double  getAngleForceEnergy   (int const &,vector3 &,bool) const;  
  void    getAngleReForceEnergy (vector3 &,bool) const;  
  double  getStretchForceEnergy (Vertex const * const,
                                 vector3 &,
                                 double const &,
                                 bool) const;
  Face *  getFace1              (void) const;
  Face *  getFace2              (void) const;

  /** Get the first vertex of this edge.
   * \return First vertex of this edge.
   */

  Vertex * getV1 (void)
  {
    return v1;
  }

  /** Get the second vertex of this edge.
   * \return Second vertex of this edge.
   */

  Vertex * getV2 (void)
  {
    return v2;
  }

  /** Get the non-edge vertex of first adjacent face.
   * \return Non-edge vertex of first adjacent face.
   */

  Vertex * getO1 (void)
  {
    return o1;
  }

  /** Get the non-edge vertex of second adjacent face.
   * \return Non-edge vertex of second adjacent face.
   */

  Vertex * getO2 (void)
  {
    return o2;
  }

  /** Get edge original length.
   * \return Edge original length.
   */

  double getOriginalLength (void) const
  {
    return l;
  }

  /** Get current squared length of this edge.
   * \return Current squared length of this edge.
   */

  double getSqLength (void) const
  {
    vector3 a(*v1->getPos()-*v2->getPos());
    return a.dot(a); 
  }

};

#endif
