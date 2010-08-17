// Author: Justin Kinney
// Date: Sep 2008

#ifndef REFRACTORY_H
#define REFRACTORY_H 1

#include <list>

#include "meshmorph.h"

typedef std::list<Vertex*>                   v_list;
typedef std::list<Vertex*>::iterator         vl_it;
typedef std::map<Vertex*,int,ltv>            map_i;
typedef std::map<Vertex*,int,ltv>::iterator  ti_it;

class Refractory
{
private:
  // Four containers (small_disp,n,face_int,ang) are used where a single
  // container would accomplish same function, but this way provides
  // more detailed diagnostics infomation if needed.
//  hmap_v  small_disp;  // vertex*s refracted because small move attempted
//  hmap_v  n;           // vertex*s refracted because moved max # times
//  hmap_v  face_int;    // vertex*s refracted because move results in face intersections
//  hmap_v  ang;         // vertex*s refracted because move results in really small or large angle
//  hmap_v  oct;         // vertex*s refracted because move results in breach of octree boundary
  vec_vp  small_disp;  // vertex*s refracted because small move attempted
  vec_vp  n;           // vertex*s refracted because moved max # times
  vec_vp  face_int;    // vertex*s refracted because move results in face intersections
  vec_vp  ang;         // vertex*s refracted because move results in really small or large angle
  vec_vp  oct;         // vertex*s refracted because move results in breach of octree boundary
//  hmap_v vert_move_distr; // imap of vertex * to # times moved
  map_i vert_move_distr; // imap of vertex * to # times moved
//  v_set   set_last_N_moved_verts;  // set  of last N vertices moved
//  v_list  list_last_N_moved_verts; // list of last N vertices moved
  static Refractory * only_one;
  Refractory                            (void);
  Refractory                            (Refractory const &);
  Refractory & operator =               (Refractory const &);
  ~Refractory                           (void);
public:
  static Refractory & instance          (void);
  int     numRecordedVertexMoves        (Vertex *);
  int     getNumVertRefractedN          (void) const; 
  int     getNumVertRefractedInt        (void) const; 
  int     getNumVertRefractedAng        (void) const; 
  int     getNumVertRefractedSmallDisp  (void) const; 
  int     getNumVertRefractedOctreeVio  (void) const; 
  bool    isRefracted                   (Vertex * const);
  bool    vertexIsMoveCandidate         (Vertex * const v);
  bool    vertexMovesAreRecorded        (Vertex *);
  void    refractVertexForSmallDispVio  (Vertex * const);
  void    refractVertexForIntVio        (Vertex * const);
  void    refractVertforAngleVio        (Vertex * const);
  void    refractVertforOctreeVio       (Vertex *v);
  void    refractVertForNumVio          (Vertex * const);
  void    updateVertMoveDistr           (Vertex * const,ti_it &);
//  void    updateLastNMovedVerts         (Vertex * const);
  void    enforceMaxdisplacement        (Vertex * const,vector3 &);
//  void    updateSuccessfulMove          (const int &,Vertex * const);
  void    updateSuccessfulMove          (Vertex * const);
  void    resetForNewGroup              (void);
//  vp_cit  getNextVertex                 (const int & group,vp_cit,bool const,bool const,bool const);
  vp_cit  getNextVertex                 (vp_cit,bool const,bool const,bool const);
//  vhm_cit beginN                        (void) const;
//  vhm_cit endN                          (void) const;
//  vhm_cit beginInt                      (void) const;
//  vhm_cit endInt                        (void) const;
//  vhm_cit beginAng                      (void) const;
//  vhm_cit endAng                        (void) const;
//  vhm_cit beginSmallDisp                (void) const;
//  vhm_cit endSmallDisp                  (void) const;
//  vhm_cit beginOctreeVio                (void) const;
//  vhm_cit endOctreeVio                  (void) const;
  vp_cit beginN                        (void) const;
  vp_cit endN                          (void) const;
  vp_cit beginInt                      (void) const;
  vp_cit endInt                        (void) const;
  vp_cit beginAng                      (void) const;
  vp_cit endAng                        (void) const;
  vp_cit beginSmallDisp                (void) const;
  vp_cit endSmallDisp                  (void) const;
  vp_cit beginOctreeVio                (void) const;
  vp_cit endOctreeVio                  (void) const;
};

#endif
