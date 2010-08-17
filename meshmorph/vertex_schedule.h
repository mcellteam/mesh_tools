// Author: Justin Kinney
// Date: Sep 2008

#ifndef VERTEX_SCHEDULE_H
#define VERTEX_SCHEDULE_H 1

#include "meshmorph.h"

class Vertex_Schedule 
{
private:
  static Vertex_Schedule * only_one;
  Vertex *    cv; // Vertex* moving now
  Vertex *  seed; // candidate vertex with largest vd
  vec_vp    vset; // set of vertex*s to move next
  vector3     pH; // destination of vertex move
  int      count; // number of moved vertices in group
  int  count_ref; // detect if any vertices from set were moved
  Vertex_Schedule                            (void);
  Vertex_Schedule                            (Vertex_Schedule const &);
  Vertex_Schedule& operator=                 (Vertex_Schedule const &);
  ~Vertex_Schedule                           (void);
public:
  static Vertex_Schedule & instance          (void);
  void          calculateMoveLocation        (Vertex * const);
  void          collectVerticesToMoveNext    (const int &);
  void          readVertexSequence           (const char *filename);

  /** Get an iterator pointing to the first in the collection
   * of vertices to move next.
   * \return Iterator pointing to the first vertex to move next.
   */

  vp_cit beginVset (void) const
  {
    return vset.begin();
  }

  /** Get an iterator pointing to one past the last in the collection
   * of vertices to move next.
   * \return Iterator pointing to one past the last vertex to move next.
   */

  vp_cit endVset (void) const
  {
    return vset.end();
  }

  /** Get number of vertices in collection to move next.
   * \return Number of vertices in collection.
   */

  int getVsetSize (void) const
  {
    return vset.size();
  }

  /** Increment the number of moved vertices in group.
  */

  void incrementNumMovedVertsGroup (void)
  {
    count++;
  }

  /** Retrieve the location to which current vertex is attempting to move.
   * \return Recorded move location of current vertex.
   */

  vector3 const * getVertexDestination (void) const
  {
    return &pH;
  }

  /** Retrieve vertex attempting to move now.
   * \return Pointer to vertex attempting to move now.
   */

  Vertex * getCurrentVertex (void) const
  {
    return cv;
  }

  /** Retrieve pointer to seed vertex of collection of vertices to move.
   * \return Pointer to seed vertex.
   */

  Vertex * getSeedVertex (void) const
  {
    return seed;
  }

  /** Retrieve the number of vertices moved so far in this group.
   * \return The number of vertices moved so far in this group.
   */

  int getNumMovedVertsGroup (void) const
  {
    return count;
  }

  /** Initialize number of moved vertices in group.
   * \param[in] init_value Initial value of number of moved vertices.
   */

  void setNumMovedVertsGroup (int const & init_value)
  {
    count = init_value;
  }

  /** Check if more vertices in collection to move.
   * \return True if no more vertices to move; false otherwise.
   */

  bool noMoreVerticesToMove (void) const
  {
    return vset.empty();
  }

  /** Retrieve the number of moved vertices in group
   * at start of last set of moved vertices.
   * \return Reference count of moved vertices in group.
   */

  int getReferenceCountMovedVerts (void) const
  {
    return count_ref;
  }

  /** Determine if any vertices were moved from last set of move candidates.
   * \return True if no vertices were moved; false otherwise.
   */

  bool noSetVerticesMoved (void) const
  {
    return count_ref==count;
  }

};

#endif
