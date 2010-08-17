// Author: Justin Kinney
// Date: Sep 2008

#ifndef STATE_H
#define STATE_H 1

#include "meshmorph.h"

#include "Vector3r.h"

typedef std::unordered_map<Edge*,int,e_hash,eqe>                    hmap_e;
typedef std::unordered_map<Edge*,float,e_hash,eqe>                 hmap_e_d;
typedef std::unordered_map<Edge*,float,e_hash,eqe>::iterator       edhm_it;
typedef std::unordered_map<Edge*,float,e_hash,eqe>::const_iterator edhm_cit;

struct Search_Stats
{
  int fs_changed; // Number of vertices that required an update by full search.
  int ps_changed; // Number of vertices that required an update by partial search.
  int fs_changed_adj; // Number of vertices that required an update by full search
  // and new closest face is adjacent to vertex being moved.
  int ps_changed_adj; // Number of vertices that required an update by partial search
  // and new closest face is adjacent to vertex being moved.
};

class State
{
private:
  static State * only_one;
  vec_ep               ae; // edge*s from adjacent faces to current vertex
  hmap_e_d             ea; // edge* -> float edge angle in radians
                           // Note sequence of elements in hash_map
                           // is not necessarily repeatable.
  float               vd2; // square of actual virtual displacement of current vertex
  State                                           (void);
  State                                           (State const &);
  State & operator =                              (State const &);
  ~State                                          (void);
public:
  static        State & instance                  (void);
  v_set         getVertsForFullClosestPtSearch    (Vertex const * const,
                                                   vector3 const &,
                                                   vector3 const &);
  v_set         getVertsForPartialClosestPtSearch (vector3 const &,
                                                   vector3 const &);
  void          recordVertAdjFaceEdgeAngles       (Vertex const * const);
  Search_Stats  updateClosestFaceToVertices       (Vertex * const,
                                                   v_set &,
                                                   v_set &);
  void          updateAdjacentFacesInTree         (Vertex const * const);
  void          updateVertexVD                    (Vertex * const);
  bool          extremeAngleFound                 (void) const;
  bool          assignNewVertexCoords             (Vertex * const,
                                                   vector3 const * const,
                                                   vector3 const &,
                                                   bool &,
                                                   bool &,
                                                   bool &); 
  bool          angleChangeIsWrong                (float const &,float const &) const;
  bool          vertexOutsideOctreeBounds         (vector3 const * const);
  void          updateVertexNormals               (v_set & fs);
  void          updateAdjacentFacesInOctree       (Vertex * const v,
                                                   hxa7241_graphics::Vector3r * const adjacent_face_lower,
                                                   hxa7241_graphics::Vector3r * const adjacent_face_upper);
  void          updateOctree                      (Vertex * const v,
                                                    hxa7241_graphics::Vector3r * const old_adjacent_face_lower,
                                                    hxa7241_graphics::Vector3r * const old_adjacent_face_upper,
                                                    hxa7241_graphics::Vector3r * const new_adjacent_face_lower,
                                                    hxa7241_graphics::Vector3r * const new_adjacent_face_upper);
  
  /** Return stored virtual displacement squared.
   * \return Squared virtual displacement.
   */

  float getVD2 (void)
  {
    return vd2;
  }

  /** Set stored virtual displacement squared.
   * \param[in] sqd_disp Squared virtual displacement.
   */

  void setVD2 (float const & sqd_disp)
  {
    vd2 = sqd_disp;
  }
};

#endif
