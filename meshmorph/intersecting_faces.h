// Author: Justin Kinney
// Date: Sep 2008

#ifndef INTERSECTING_FACES_H
#define INTERSECTING_FACES_H 1

#include "meshmorph.h"

//typedef __gnu_cxx::hash_map<Face*,vec_fp,f_hash,eqf>                 hmap_f_f;
//typedef __gnu_cxx::hash_map<Face*,vec_fp,f_hash,eqf>::iterator       htff_it;
//typedef __gnu_cxx::hash_map<Face*,vec_fp,f_hash,eqf>::const_iterator htff_cit;
//typedef __gnu_cxx::hash_set<Vertex*,v_hash,eqv>::iterator            hv_it;
typedef std::unordered_map<Face*,vec_fp,f_hash,eqf>                 hmap_f_f;
typedef std::unordered_map<Face*,vec_fp,f_hash,eqf>::iterator       htff_it;
typedef std::unordered_map<Face*,vec_fp,f_hash,eqf>::const_iterator htff_cit;
typedef std::unordered_set<Vertex*,v_hash,eqv>::iterator            hv_it;

class Intersecting_Faces
{
private:
  hmap_f_f intf; // store intersecting faces
                 // Note sequence of elements in hash_map
                 // is not necessarily repeatable.
  static Intersecting_Faces * only_one;
  Intersecting_Faces                     (void);
  Intersecting_Faces                     (Intersecting_Faces const &);
  Intersecting_Faces & operator =        (Intersecting_Faces const &);
  ~Intersecting_Faces                    (void);
public:
  static Intersecting_Faces & instance   (void);
  vec_fp *  getIntersectingFacesRHS      (Face const * const);
  int       getNumUniqueVerts            (Face const * const,
                                          Face const * const,
                                          int * const) const;
  int       getCountOfIntFaces           (bool);
  bool      checkFaceEdgeInt             (Face const * const,
                                          Face const * const) const;
  bool      checkEdgeEdgeIntersection    (Face const * const,
                                          Face const * const) const;
  bool      facesParallel                (Face const * const,
                                          Face const * const) const;
  bool      facesCoplanar                (Face const * const,
                                          Face const * const) const;
  bool      faceIntersectsFace           (Face const * const haystack,
                                          Face const * const needle);
  bool      faceIsIntersectedRHS         (Face const * const) const;
  bool      findAndRecordNewFaceInt      (Face * const);
  bool      vertAdjFacesHaveNewInt       (Vertex const * const);
  bool      checkFaceFaceInts            (Face const * const,
                                          Face const * const) const;
  bool      detectNewFaceInts            (Face * const,vec_fp &);
  void      getNiceSet                   (v_set &,hashset_v &);
  hashset_v getNiceCheckSet              (Vertex const * const);
  void      addFaceToFace                (Face * const,Face const * const);
  void      removeOldIntersections       (Vertex const * const,hashset_v &);
  void      updateNewIntersections       (Vertex const * const,hashset_v &);
  void      findAllFaceIntersections     (void);
  void      setFaceNotIntersectedLHS     (Face const * const);
  void      removeFaceFromFaceInt        (Face * const haystack,
                                          Face const * const needle);
  bool      intFacesAreSymmetric         (void);
  // intersection force
  void     getFaceIntersectionForce      (Face * const f,
                                          vector3 & total_force);
  void     getFaceFaceIntForce           (Face * const lhs,
                                          Face * const rhs,
                                          vector3 &);

  /** Weakly determine if input face is recorded in class as being intersected.
   * \param[in] face Face of interest.
   * \return True if face is used as key in intersecting faces container
   * (with no further checks); otherwise false.
   */

  bool faceIsIntersectedLHS (Face const * const face) const
  {
    return intf.find(const_cast<Face*>(face))!=intf.end();
  }

  /** Get an iterator pointing to first intersecting face.
   * \return Iterator pointing to beginning of intersecting faces container.
   */

  htff_cit begin (void)
  {
    return intf.begin();
  }

  /** Get an iterator pointing to one past the last intersecting face.
   * \return Iterator pointing to one past the last intersecting face.
   */

  htff_cit end (void)
  {
    return intf.end();
  }

};

#endif
