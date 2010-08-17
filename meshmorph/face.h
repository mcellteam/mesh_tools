// Author: Justin Kinney
// Date: Sep 2008

#ifndef FACE_H
#define FACE_H 1

#include <cmath>

#include "meshmorph.h"

#include "Vector3r.h"
#include "vertex.h"

class Face
{
private:
  int      index; // from file, unaltered
  bool      flag; // true if face was returned from octree; false otherwise
  Edge *    e[3]; // pointers to edges
  Vertex *  v[3]; // pointers to vertices
  vector3      n; // normal vector
public:
  Face &  operator =            (Face const &);
//  Face                          (char *,vec_vp &); 
  Face                          (int const &,int const &,int const &,int const &,vec_vp &); 
  Face                          (Face const &);
  bool    isMatch               (int i, std::string const & name) const;
  void    print                 (std::ostream &) const;
  void    printCP               (std::ostream &) const;
  void    addEdge               (Edge *);
  void    getBoundingBox        (hxa7241_graphics::Vector3r & lower,
                                 hxa7241_graphics::Vector3r & upper) const;
  void    getVertexCoord        (vector3 const * &,
                                 vector3 const * &,
                                 vector3 const * &) const;
  double  getAngle              (Vertex const * const v) const;
  double  getAngleProxy         (Vertex const * const) const;
  double  getAspectRatio        (vector3 &,Vertex const *&) const;
  double  getAspectRatioForceEnergy (Vertex const * const v,
                                    vector3 & force,
                                    bool compute_force) const;
  void    updateNormal          (void);
  Object const * getObject (void) const;


  /** Get flag of this face.
   * \return Flag value.
   */

  bool getFlag (void) const
  {
    return flag;
  }

  /** Set flag of this face to true.
   */

  void setFlag (void)
  {
    flag = true;
  }

  /** Reset flag of this face to false.
   */

  void clearFlag (void)
  {
    flag = false;
  }

  /** Get face index.
   * \return The face index.
   */

  int getIndex (void) const
  {
    return index;
  }

  /** Get one of edges of this face.
   * \param[in] i Edge of interest.
   * \return The ith edge.
   */

  Edge * getEdge (int const & i)
  {
    return e[i];
  }

  /** Get pointer to face vertex.
   * \param[in] i Index (0,1,2) of face vertex of interest.
   * \return Pointer to vertex of this face.
   */

  Vertex * getVertex (int const & i) const
  {
    return v[i];
  }

  /** Return face normal (not guaranteed to be unit vector).
   * \return Face normal.
   */

  vector3 const *getNormal (void) const
  {
    return &n;
  }

  /** Calculate face area.
   * \return Face area.
   */

  double computeArea (void)
  {
    // compute face area = half normal vector length
    return sqrt(n.dot(n))/2.0;
  }
};

#endif
