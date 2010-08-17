// Author: Justin Kinney
// Date: Sep 2008

#ifndef NICE_H
#define NICE_H 1

#include "meshmorph.h"

typedef std::map<double,int,ltd>            map_di;
typedef std::map<double,int,ltd>::iterator   di_it;
//typedef __gnu_cxx::hash_map<Vertex*,int,v_hash,eqv>::iterator        vhm_it;
typedef std::unordered_map<Vertex*,int,v_hash,eqv>::iterator        vhm_it;

struct face_grp
{
  vec_fp crossed_faces;
  vec_fp edge_faces;
  face_grp (void) :crossed_faces(),edge_faces() { }
};

class Nice
{
private:
  hmap_v nonnice; // vertex* is key to int nice value
                  // Note sequence of elements in hash_map
                  // is not necessarily repeatable.
  static Nice * only_one;
  Nice                                    (void);
  Nice                                    (Nice const &);
  Nice & operator =                       (Nice const &);
  ~Nice                                   (void);
public:
  static Nice & instance                  (void);
  int       getNonniceCount               (bool);
  int       getVertexNiceVal              (Vertex const * const) const;
  bool      findExtraPoint                (Vertex const * const,
                                           vector3 &,
                                           fp_cit const &) const;
  bool      updateVertexNiceness          (Vertex * const);
  bool      setVertexNiceness             (Vertex * const,vec_op &);	
  bool      getCrossedObjFromVertToExtra  (Vertex const * const,
                                           vector3 const &,
                                           vec_op &) const;
  bool      getCrossedObjFromExtraToLimit (vector3 const &,
                                           vector3 const &,
                                           vec_op &) const;
  bool      vertexIsNice                  (Vertex const * const) const;
  bool      getPointOutsideObject         (Vertex const * const,
                                           vector3 &,
                                           vec_op &) const;
  void      getCrossedObjects             (Vertex const * const,
                                           vec_op &) const;
  void      getVertAdjFaceRay             (vector3 &,
                                           vector3 &,
                                           fp_cit const &) const;
  void      getPenetratedObjs             (vec_fp &,
                                           vec_fp &,
                                           vec_op &) const;
  void      getObjectsFromEdgeHits        (vec_op & edge_hits,
                                           vec_op & objs) const;
  void      findOddObjects                (vec_op &,vec_op &) const;
  void      getRaysToWorldLimit           (vector3 const &,double [6][3]) const;
  void      setVertexNiceVal              (int const &,Vertex * const);
  void      findNonniceVertices           (void);
  face_grp  findIntFacesAlongRay          (Vertex const * const,
                                           vector3 const &,
                                           vector3 const &,
                                           bool) const;
  /** Get an iterator pointing to first nonnice vertex.
   * \return Iterator pointing to beginning of nonnice vertex container.
   */

  vhm_cit beginNice (void) const
  {
    return nonnice.begin();
  }

  /** Get an iterator pointing to one past the last nonnice vertex.
   * \return Iterator pointing to one past the last nonnice vertex.
   */

  vhm_cit endNice (void) const
  {
    return nonnice.end();
  }
};

#endif
