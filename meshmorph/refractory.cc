// Author: Justin Kinney
// Date: Sep 2008

#include "refractory.h"

#include <cmath>
#include <iostream>

#include "log.h"
#include "edge.h"
#include "state.h"
#include "controls.h"
#include "container.h"
#include "gain_schedule.h"

using std::cout;
using std::endl;

Refractory * Refractory::only_one = NULL;

Refractory & Refractory::instance()
{
  // Not thread-safe.
  // -- lock mutex
  if (only_one == NULL)
    only_one = new Refractory();
  // -- unlock mutex
  return *only_one;
}

Refractory::Refractory (void)
:small_disp(),n(),face_int(),
  ang(),oct(),vert_move_distr()
//  set_last_N_moved_verts(),list_last_N_moved_verts()
{
}

/** Determine if input vertex is currently in refractory period.
 * \param[in] v Vertex of interest.
 * \return True if refracted; false otherwise.
 */

bool Refractory::isRefracted (Vertex * const v)
{
//  return (       n.find(v)!=n.end() ||
//                 small_disp.find(v)!=small_disp.end() ||
//                 face_int.find(v)!=face_int.end() ||
//                 oct.find(v)!=oct.end() ||
//                 //com.find(v)!=com.end() || 
//                 ang.find(v)!=ang.end() ||
//                 set_last_N_moved_verts.find(v)!=set_last_N_moved_verts.end());
  return v->getLastRefractoryIter() && 
         (v->getLastRefractoryIter() > Log::instance().getN());
           //set_last_N_moved_verts.find(v)!=set_last_N_moved_verts.end());
}

/** Put input vertex into refractory state (period when vertex cannot move)
 * for attempting such a small move.
 * \param[in] v Vertex to be refracted.
 */

void Refractory::refractVertexForSmallDispVio (Vertex * const v)
{
  Controls & cs(Controls::instance());
  //small_disp[v]=cs.get_group_size();
  small_disp.push_back(v);
  v->setLastRefractoryIter(Log::instance().getN()+cs.get_group_size());
  if (cs.get_disable_messages()==false)
  {
    cout << "\nVertex refracted because smallest virtual displacement ("
          << sqrt(State::instance().getVD2())
          << ") is unacceptably small.\n";
    v->print(cout);
    cout << endl;
  }
}

/** Put input vertex into refractory state (period when vertex cannot move)
 * for resulting in face intersections with every attempted move.
 * \param[in] v Vertex to be refracted.
 */

void Refractory::refractVertexForIntVio (Vertex * const v)
{
  Controls & cs(Controls::instance());
  // get_group_size() is duration of refractory period.
  //face_int[v]=cs.get_group_size();
  face_int.push_back(v);
  v->setLastRefractoryIter(Log::instance().getN()+cs.get_group_size());
  if (cs.get_disable_messages()==false)
  {
    cout << "\nVertex refracted because even smallest virtual displacement ("
          << sqrt(State::instance().getVD2())
          << ") results in face intersections.\n";
    v->print(cout);
    cout << endl;
  }
}

/** Put input vertex into refractory state (period when vertex cannot move)
 * for passing the octree boundary limits
 * with every attempted move.
 * \param[in] v Vertex to be refracted.
 */

void Refractory::refractVertforOctreeVio(Vertex *v)
{
  Controls & cs(Controls::instance());
  // get_group_size() is duration of refractory period.
  //oct[v]=cs.get_group_size();
  oct.push_back(v);
  v->setLastRefractoryIter(Log::instance().getN()+cs.get_group_size());
  if (cs.get_disable_messages()==false)
  {
    cout << "\nVertex refracted because octree boundary was breached.\n";
    v->print(cout);
    cout << endl;
  }
}

/** Put input vertex into refractory state (period when vertex cannot move)
 * for ulting in very small (near 0) or very large (near 2*pi) angles
 * with every attempted move.
 * \param[in] v Vertex to be refracted.
 */

void Refractory::refractVertforAngleVio (Vertex *v)
{
  Controls & cs(Controls::instance());
  // get_group_size() is duration of refractory period.
  //ang[v]=cs.get_group_size();
  ang.push_back(v);
  v->setLastRefractoryIter(Log::instance().getN()+cs.get_group_size());
  if (cs.get_disable_messages()==false)
  {
    cout << "\nVertex refracted because smallest virtual displacement ("
          << sqrt(State::instance().getVD2())
          << ") results in extreme angles.\n";
    v->print(cout);
    cout << endl;
  }
}

///** Update distribution of number of times each vertex has been moved.
// * \param[in] v Last moved vertex.
// */
//
//void Refractory::updateVertMoveDistr (Vertex * const v,ti_it & target)
//{
//  // if vertex in vert_move_distr, then increment
//  if (target!=vert_move_distr.end())
//  //if (vert_move_distr.find(v) != vert_move_distr.end())
//  {
//    //vert_move_distr[v]++;
//    (*target).second++;
//    //assert((*target).second < Controls::instance().get_max_touches())
//  }
//  // else init to one touch
//  else
//  {
//    vert_move_distr[v] = 1;
//  }
//}

/** Put input vertex into refractory state (period when vertex cannot move)
 * for being moved too many times.
 * \param[in] v Vertex to be refracted.
 */

void Refractory::refractVertForNumVio (Vertex * const v)
{
  Controls & cs(Controls::instance());
  // GROUP_SIZE is duration of refractory period.
  //n[v] = cs.get_group_size();
  n.push_back(v);
  v->setLastRefractoryIter(Log::instance().getN()+cs.get_group_size());
  if (cs.get_disable_messages()==false)
  {
    cout << "\nVertex refracted because it moved its quota: "
          << Controls::instance().get_max_touches() << ".\n";
    v->print(cout);
    cout << endl;
  }
}

///** Update collections of last N moved vertices.
// * \param[in] v Last moved vertex.
// */
//
//void Refractory::updateLastNMovedVerts (Vertex * const v)
//{
//  // erase oldest vertex
//  uint rp = static_cast<uint>(Controls::instance().get_refractory_period());
//  if (list_last_N_moved_verts.size()==rp)
//  {
//    set_last_N_moved_verts.erase(list_last_N_moved_verts.front());
//    list_last_N_moved_verts.pop_front();
//  }
//  // store new vertex
//  list_last_N_moved_verts.push_back(v);
//  set_last_N_moved_verts.insert(v);
//}

//bool isBetweenBoxes( Vertex * const v)
//{
//  if (
//      (
//       (*v->getCoord(0) < 60 && *v->getCoord(0) > 2700 ) &&
//       (*v->getCoord(1) < 80 && *v->getCoord(1) > 2860 ) &&
//       (*v->getCoord(2) < 0 && *v->getCoord(2) > 2900 )
//      )  
//      &&
//      (
//       (*v->getCoord(0) > -800 && *v->getCoord(0) < 3500 ) &&
//       (*v->getCoord(1) > -500 && *v->getCoord(1) < 3500 ) &&
//       (*v->getCoord(2) > -100 && *v->getCoord(2) < 3250 )
//      )  
//     )
// {
//    return true;
//  }
//  else
//  {
//    return false;
//  }
//}
//
///** Determine if vertex is allowed to move.
// * \param[in] v Vertex of interest.
// * \return True if allowed to move; false otherwise.
// */
//
//bool Refractory::vertexIsMoveCandidate (Vertex * const v)
//{
//  // if vertex is candidate, i.e. closest point was found for vertex
//  // and vertex not found in refractory
//  // and vertex is not frozen
//  return (v->getClosestFace()!=NULL) && (isRefracted(v)==false)
//        && (Container::instance().vertexIsFrozen(v)==false)
//        && (isBetweenBoxes(v)==true)
//        ;
//  //return (isRefracted(v)==false)
//  //      && (Container::instance().vertexIsFrozen(v)==false);
//}

/** Enforce maximum vertex displacment policy.
 * \param[in] v Vertex being moved.
 * \param[out] new_pos New Position of input vertex.
 */

bool Refractory::vertexIsMoveCandidate (Vertex * const v)
{
  // if vertex is candidate, i.e. closest point was found for vertex
  // and vertex not found in refractory
  // and vertex is not frozen
  return (v->getClosestFace()!=NULL) && (isRefracted(v)==false)
        && (Container::instance().vertexIsFrozen(v)==false)
        && (!Controls::instance().get_freeze_sheets()  || v->isSheet()==false)
        && (!Controls::instance().get_freeze_tunnels() || v->isSheet()==true)
        //&& (isBetweenBoxes(v)==true)
        ;
  //return (isRefracted(v)==false)
  //      && (Container::instance().vertexIsFrozen(v)==false);
}

/** Enforce maximum vertex displacment policy.
 * \param[in] v Vertex being moved.
 * \param[out] new_pos New Position of input vertex.
 */

void Refractory::enforceMaxdisplacement (Vertex * const v,
                                         vector3 & new_pos)
{
  // displacement vector
  vector3 d(new_pos-*v->getPos());
  // square of displacement
  double disp = sqrt(d.dot(d));
  // find shortest edge length from adjacent faces
  double min_edge_length = 1e30;
  // for each adjacent face of vertex
  for (fp_cit i=v->begin();i!=v->end();++i)
  {
    // find shortest edge
    for (int j=0;j<3;++j)
    {
      if ((*i)->getEdge(j)->getOriginalLength() < min_edge_length)
      {
        min_edge_length = (*i)->getEdge(j)->getOriginalLength();
      }
    }
  }
  // compute max allowed displacement
  double max_allowed_disp = min_edge_length
        *Controls::instance().get_max_actual_displ_fraction();
  // if too big
  if (disp>max_allowed_disp)
  {
    // scale down displacement vector
    d *= max_allowed_disp/disp;
    // compute new position of vertex after move
    new_pos = *v->getPos()+d;
  }
}

/** Decide whether to try moving the same vertex again
 * or give up and try moving another vertex.
 * \param[in] v Vertex of interest.
 * \param[in] int_flag True if vertex move resulted in new face intersections;
 * false otherwise.
 * \param[in] ang_flag True if vertex move resulted in extreme angles;
 * false otherwise.
 * \param[in] outside_octree True if vertex move breached octree boundary;
 * false otherwise.
 * \return Iterator pointing to next vertex move candidate.
 */

//vp_cit Refractory::getNextVertex (const int & group,
vp_cit Refractory::getNextVertex (vp_cit v,
                                  bool const int_flag,
                                  bool const ang_flag,
                                  bool const outside_octree)
{
  Controls & cs(Controls::instance());
  // NOTE THIS STRATEGY FOR HANDLING FAILED MOVES
  // IS QUITE SIMPLISTIC
  // first try shortening move
  if (State::instance().getVD2()>cs.get_min_displacement_sq())
  {
    // try halving gain to halve attempted virtual displacement and try again
    Gain_Schedule::instance().halveGain();
    return v;
  }
  // attempted move is already short
  bool detect = false;
  // refract if faces were intersected
  if (int_flag==true)
  {
    refractVertexForIntVio(*v);
    //Log::instance().writeRefractedNow(group,4);
    detect = true;
  }
  // refract if octree boundary was breached
  if (outside_octree==true)
  {
    refractVertforOctreeVio(*v);
    //Log::instance().writeRefractedNow(group,5);
    detect = true;
  }
  // refract if extreme angles were generated
  if (ang_flag==true)
  {
    refractVertforAngleVio(*v);
    //Log::instance().writeRefractedNow(group,2);
    detect = true;
  }
  // refract for small virtual displacement
  if (detect==false)
  {
    refractVertexForSmallDispVio(*v);
    //Log::instance().writeRefractedNow(group,3);
  }
  // reset gain
  Gain_Schedule::instance().resetGain();
  // move on to next vertex in set
  return ++v;
}

/** Update class after successful vertex move.
 * \param[in] v Last moved vertex.
 */

//void Refractory::updateSuccessfulMove (const int & group,Vertex * const v)
void Refractory::updateSuccessfulMove (Vertex * const v)
{
  Controls & cs(Controls::instance());
  ti_it target = vert_move_distr.find(v);
  // if vertex already in vert_move_distr
  if (target!=vert_move_distr.end())
  {
    // record vertex as having moved
    (*target).second++;
    // vertex has moved its allowed number of times
    if ((*target).second > cs.get_max_touches())
    {
      if (cs.get_disable_messages()==false)
      {
        cout << "vert_move_distr[v] = "
              << (*target).second
              << ", Controls::instance().get_max_touches() = "
              << Controls::instance().get_max_touches() << endl;
      }
      // place vertex into long refraction
      refractVertForNumVio(v);
      // recrd refraction
      //Log::instance().writeRefractedNow(group,1);
    }
    else
    {
      // refract vertex for having moved
      v->setLastRefractoryIter(Log::instance().getN()+cs.get_refractory_period());
    }
  }
  else
  {
    // else init to one touch
    vert_move_distr[v] = 1;
    // refract vertex for having moved
    v->setLastRefractoryIter(Log::instance().getN()+cs.get_refractory_period());
  }
}

/** Clear class members for new group of moved vertices.
*/

void Refractory::resetForNewGroup (void)
{
  small_disp.clear();
  n.clear();
  face_int.clear();
  ang.clear();
  vert_move_distr.clear();
  // for each object in container
  Container & c(Container::instance());
  for (o_it i=c.o.begin();i!=c.o.end();++i)
  {
    // for each vertex in object
    for (v_it j=i->v.begin();j!=i->v.end();++j)
    {
      j->setLastRefractoryIter(0);
    }
  }
}

/** Check if any moves by input vertex are recorded in this class.
 * \param[in] v Vertex of interest.
 * \return True if recorded moves found for vertex.
 */

bool Refractory::vertexMovesAreRecorded (Vertex * v)
{
  return vert_move_distr.find(v)!=vert_move_distr.end();
}

/** Get number of recorded moves for input vertex.
 * \param[in] v Vertex of interest.
 * \return Number of recorded moves found for vertex.
 */

int Refractory::numRecordedVertexMoves (Vertex * v)
{
  ti_it target = vert_move_distr.find(v);
  //return vert_move_distr[v];
  return (*target).second;
}

/** Get number of vertices currently in refraction
 * for filling their quota of moves per group.
 * \return Number of vertices refracted for moving MAX_TOUCHES times. 
 */

int Refractory::getNumVertRefractedN (void) const
{
  return n.size();
}

/** Get number of vertices currently in refraction
 * for breaching boundarsy of octree.
 * \return Number of vertices refracted for breaching octree boundary. 
 */

int Refractory::getNumVertRefractedOctreeVio (void) const
{
  return oct.size();
}

/** Get an iterator pointing to the first in the collection
 * of vertices refracted for breaching boundary of octree.
 * \return Iterator pointing to the first refracted vertex.
 */

vp_cit Refractory::beginOctreeVio (void) const
{
  return oct.begin();
}

/** Get an iterator pointing to one past the last in the collection
 * of vertices refracted for breaching boundary of octree.
 * \return Iterator pointing to one past the last refracted vertex.
 */

vp_cit Refractory::endOctreeVio (void) const
{
  return oct.end();
}

/** Get an iterator pointing to the first in the collection
 * of vertices refracted for moving MAX_TOUCHES times.
 * \return Iterator pointing to the first refracted vertex.
 */

vp_cit Refractory::beginN (void) const
{
  return n.begin();
}

/** Get an iterator pointing to one past the last in the collection
 * of vertices refracted for moving MAX_TOUCHES times.
 * \return Iterator pointing to one past the last refracted vertex.
 */

vp_cit Refractory::endN (void) const
{
  return n.end();
}

/** Get number of vertices currently in refraction
 * for resulting in intersected faces.
 * \return Number of vertices refracted for resulting in intersected faces. 
 */

int Refractory::getNumVertRefractedInt (void) const
{
  return face_int.size();
}

/** Get an iterator pointing to the first in the collection
 * of vertices refracted for resulting in intersected faces.
 * \return Iterator pointing to the first refracted vertex.
 */

vp_cit Refractory::beginInt (void) const
{
  return face_int.begin();
}

/** Get an iterator pointing to one past the last in the collection
 * of vertices refracted for resulting in intersected faces.
 * \return Iterator pointing to one past the last refracted vertex.
 */

vp_cit Refractory::endInt (void) const
{
  return face_int.end();
}

/** Get number of vertices currently in refraction
 * for resulting in extreme angles.
 * \return Number of vertices refracted for resulting in extreme angles. 
 */

int Refractory::getNumVertRefractedAng (void) const
{
  return ang.size();
}

/** Get an iterator pointing to the first in the collection
 * of vertices refracted for resulting in extreme angles.
 * \return Iterator pointing to the first refracted vertex.
 */

vp_cit Refractory::beginAng (void) const
{
  return ang.begin();
}

/** Get an iterator pointing to one past the last in the collection
 * of vertices refracted for resulting in extreme angles.
 * \return Iterator pointing to one past the last refracted vertex.
 */

vp_cit Refractory::endAng (void) const
{
  return ang.end();
}

/** Get number of vertices currently in refraction
 * for having unacceptably small virtual displacements.
 * \return Number of vertices refracted for resulting
 * in small virtual displacements.. 
 */

int Refractory::getNumVertRefractedSmallDisp (void) const
{
  return small_disp.size();
}

/** Get an iterator pointing to the first in the collection
 * of vertices refracted for having an unacceptably small virtual displacement.
 * \return Iterator pointing to the first refracted vertex.
 */

vp_cit Refractory::beginSmallDisp (void) const
{
  return small_disp.begin();
}

/** Get an iterator pointing to one past the last in the collection
 * of vertices refracted for having an unacceptably small virtual displacement.
 * \return Iterator pointing to one past the last refracted vertex.
 */

vp_cit Refractory::endSmallDisp (void) const
{
  return small_disp.end();
}

