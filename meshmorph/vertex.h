// Author: Justin Kinney
// Date: Sep 2008

#ifndef VERTEX_H
#define VERTEX_H 1

#include <algorithm>
#include <cassert>

#include "meshmorph.h"

#include "face.h"
#include "object.h"

class Vertex
{
private:
  int     index; // from file, unaltered
  Face *     cl; // pointer to face on which closest mesh position lies
  Object *    o; // pointer to parent object
  vec_fp      f; // pointers to adjacent faces
  vector3     n; // vertex normal
  vector3     p; // current position coordinates (x,y,z)
  int         r; // last refractory iteration
  int         m; // num different objects encountered in closest point search
public:
  Vertex                                (int const &,double const & x,double const & y,double const & z,Object * const);
//  Vertex                                (int,double x,double y,double z,Object * const);
  //Vertex                                (char const * triplet,Object * const);
  Vertex                                (Vertex const &);
  Vertex & operator =                   (Vertex const &);
  void    print                         (std::ostream &) const;
  void    printCP                       (std::ostream &) const;
  void    setNewPos                     (vector3 const * const);
  void    setNormal                     (void);
  void    getBoundingBox                (vector3 &,vector3 &) const;
  void    updateAdjFaceBoundingBoxes    (void);
  void    defineLocalRegion             (vector3 &,vector3 &);
  bool    surfaceIsConcave              (const vector3 & radial_vector) const;
  void    getAdjVertices                (vec_vp&) const;
  void    getAdjVerticesMulti           (vec_vp & expanded_verts,
                                         const int num_expansions) const;
  void    getAdjacentEdges              (vec_ep&) const;
  double  getArea                       (void) const;
  double  getSqSepDist                  (void) const; 
  double  getSqVirtualDisp              (double);
  double  getRadiusOfCurvature          (void) const;
  std::string  getRadiusOfCurvatureDEBUG     (void) const;
  vector3 getNewPos                     (double);
  void    recordAdjFaceBoundingBoxes    (hxa7241_graphics::Vector3r * const adjacent_face_lower,
                                         hxa7241_graphics::Vector3r * const adjacent_face_upper);
  // force
  void    getTotalForceEnergy           (vector3 &) const;
  void    getAdjFaceIntForce            (vector3 &) const;
  double  getEcwForceEnergy             (vector3 &,bool) const;
  double  getEdgeStretchForceEnergy     (vector3 &,bool) const;
  double  getFaceAspectRatioForceEnergy (vector3 &,bool) const;
  double  getEdgeAngleForceEnergy       (vector3 &,bool) const;
  Face *  getFaceNotAdjToVertex         (Vertex const * const) const;
  void    updateAdjacentFaceNormals     (void);

  /** Determine if the extracellular space in front of this vertex
   * is part of a sheet or tunnel based on the number of neighboring
   * objects.
   *
   *  \return True if vertex is part of a sheet; otherwise, false
   *    means vertex is part of a tunnel.
   */

  bool isSheet (void) const
  {
    assert (m>-1);
    return m<2;
  }

  int getLastRefractoryIter (void)
  {
    return r;
  }

  void setLastRefractoryIter (int i)
  {
    r = i;
  }

  /** Sort stored adjacent faces to this vertex by face pointer.
  */

  void sortAdjacentFaces (void)
  {
    sort(f.begin(),f.end());
  }

  /** Set number of different objects encountered during search for
   *  closest point to this vertex.
   *
   *  \param[in] i Number of different objects.
   */

  void setNeighborCount (int i)
  {
    m = i;
  }

  /** Get stored number of different objects encountered during search for
   *  closest point to this vertex.
   *
   *  \return Number of different objects.
   */

  int getNeighborCount (void)
  {
    return m;
  }

  /** Get stored vertex normal vector.
   *
   *  \return Vertex normal vector. Note vector length
   *  is not necessarily of unit length.
   */

  vector3 getNormal (void) const
  {
    return n;
  }

  /** Get object name.
   *
   *  \return Name of parent object of this vertex. 
   */

  std::string const & getObjectName (void) const
  {
    return o->getName();
  }

  /** Return the index of vertex as recorded from input file. 
   *
   * \return Vertex index.
   */

  int getIndex (void) const
  {
    return index;
  }

  /** Calculate and return the number of stored adjacent faces to this vertex.
   * \return Number of stored adjacent faces to this vertex.
   */

  int getNumAdjFaces (void) const
  {
    return f.size();
  }

  /** Determine if input face is adjacent to this vertex.
   * \param[in] face Face of interest.
   * \return True if input face is found among stored adjcent faces
   * of this vertes; false otherwise.
   */

  bool faceIsAdjacent (Face * face) const
  {
    return binary_search(f.begin(),f.end(),face);
  }

  /** Compare input vertex index and object name
   *  to this vertex's index and object name.
   * \param[in] i Input vertex index.
   * \param[in] name Input vertex parent object name.
   * \return True if indices and names match; false otherwise.
   */

  bool isMatch (int const & i, std::string const & name) const
  {
    return i==index && name==o->getName();
  }

  /** Record input face as adjacent to this vertex.
   * \param[in] face New adjacent face to this vertex.
   */

  void addAdjacentFace (Face * face)
  {
    f.push_back(face);
  }

  /** Retrieve one coordinate of position of this vertex.
   * \param[in] axis Axis of interest (0,1,2==x,y,z axis respectively).
   * \return The position coordinate of this vertex with index equal to axis.
   */

  double const * getCoord (int const & axis) const
  {
    return &p.p[axis];
  }

  /** Retrieve the stored face on which
   * the closest point to this vertex should lie.
   * \return Stored closest face pointer.
   */

  Face * getClosestFace (void)
  {
    return cl;
  }

  /** Record the stored face on which
   * the closest point to this vertex should lie.
   * \param[in] face Pointer to closet face to this vertex..
   */

  void setFace (Face * const face)
  {
    cl=face;
  }

  /** Retrieve the current position of this vertex.
   * \return Pointer to this vertex.
   */

  vector3 const * getPos (void) const
  {
    return &p;
  }

  /** Set the position of this vertex.
   * \param[in] x X coordinate.
   * \param[in] y Y coordinate.
   * \param[in] z Z coordinate.
   */

  void setPos (double const & x, double const & y, double const & z)
  {
    p.p[0] = x;
    p.p[1] = y;
    p.p[2] = z;
  }

  /** Retrieve a pointer to the parent Object of this vertex.
   * \return Pointer to parent object of this vertex.
   */

  Object const * getObject (void) const
  {
    return o;
  }

  /** Get an iterator pointing to the first in the collection
   * of adjacent faces to this vertex.
   * \return Iterator pointing to the first adjacent face.
   */

  fp_cit begin (void) const
  {
    return f.begin();
  }

  /** Get an iterator pointing to one past the last in the collection
   * of adjacent faces to this vertex.
   * \return Iterator pointing to one past the last adjacent face.
   */

  fp_cit end (void) const
  {
    return f.end();
  }
};

struct my_ltv
{
  bool operator () (const Vertex * s1, const Vertex * s2) const
  {
    return s1->getIndex()<s2->getIndex();
  }
};

#endif
