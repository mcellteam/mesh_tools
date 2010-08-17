// Author: Justin Kinney
// Date: Sep 2008

#include "state.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "log.h"
#include "edge.h"
#include "controls.h"
#include "container.h"
#include "virtual_disp.h"
#include "gain_schedule.h"
#include "octree_agent_face.h"
#include "octree_visitor_update.h"
#include "intersecting_faces.h"

using std::cout;
using std::endl;

// DEBUG
//const int TARGET_VERTEX_INDEX = 58340;
//const std::string TARGET_VERTEX_NAME("d000");
// DEBUG
  
State * State::only_one = NULL;

State & State::instance(void)
{
  // Not thread-safe.
  // -- lock mutex
  if (only_one == NULL)
    only_one = new State();
  // -- unlock mutex
  return *only_one;
}

State::State (void)
  :ae(),ea(),vd2(0)
{
  ae.reserve(Controls::instance().get_vector_reserve());
};

/** Collect vertices whose closest point may now lie
 *  anywhere.
 * \param[in] v The current vertex.
 * \param[in] lower Lower corner of desired search region.
 * \param[in] upper Upper corner of desired search region.
 * \return Vertex pointers requiring full closest point search.
 */

v_set State::getVertsForFullClosestPtSearch (Vertex const * const v,
                                                vector3 const & lower,
                                                vector3 const & upper)
{
  // make a visitor
  Octree_Visitor_Face visitor(Vector3r(lower.p[0],lower.p[1],lower.p[2]),
                              Vector3r(upper.p[0],upper.p[1],upper.p[2]));
  // execute visitor
  Container::instance().octree->visit( visitor );
  //visitor.sort();
  v_set fs;
    // for each face in box
    for (fp_cit j=visitor.mybegin();j!=visitor.myend();++j)
    {
      (*j)->clearFlag();
      for (int k=0;k<3;++k)
      {
        Vertex * vv = (*j)->getVertex(k);
        // DEBUG
//        if (vv->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
//        {
//          cout << "\nState::getVertsForFullClosestPtSearch: "
//                << "Processing target vertex.\n";
//          if (vv->getClosestFace()==NULL) 
//          {
//            cout << "\nState::getVertsForFullClosestPtSearch: "
//                  << "target vertex closest face == NULL.\n";
//          }
//          else
//          {
//            cout << "\nState::getVertsForFullClosestPtSearch: "
//                  << "target vertex closest face != NULL.\n";
//          }
//        }
//        // DEBUG
        // if face vertex has a closest point and
        // closest point lies on an adjacent face of active vertex, v
        // then add vertex to affected vertex set
        if (vv->getClosestFace()!=NULL && 
            v->faceIsAdjacent(vv->getClosestFace())==true )
        {
          fs.insert(vv);
          // DEBUG
//          if (vv->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
//          {
//            cout << "\nState::getVertsForFullClosestPtSearch: "
//                  << "Target vertex added to full set.\n";
//          }
          // DEBUG
        }
//        else
//        {
//          // DEBUG
//          if (vv->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
//          {
//            cout << "\nState::getVertsForFullClosestPtSearch: "
//                  << "Target vertex NOT added to full set.\n";
//          }
//          // DEBUG
//        }
      }
    }
  // adjacent vertices to current vertex need full search for closest point
  // since vertx normal will have changed
  // for each adjacent face of vertex attempting move
  for (fp_cit k=v->begin();k!=v->end();++k)
  {
    for (int i=0;i<3;++i)
    {
      fs.insert((*k)->getVertex(i));
    }
  } 
  return fs;
}

/** Collect vertices whose closest point may now lie
 *  on an adjacent face of current vertex.
 *
 * \param[in] lower Lower corner of desired search region.
 * \param[in] upper Upper corner of desired search region.
 * \return Vertex pointers requiring partial closest point search.
 */

v_set State::getVertsForPartialClosestPtSearch (vector3 const & lower,
                                                vector3 const & upper)
{
  v_set ps;
  // make a visitor
  Octree_Visitor_Face visitor(Vector3r(lower.p[0],lower.p[1],lower.p[2]),
                              Vector3r(upper.p[0],upper.p[1],upper.p[2]));
  // execute visitor
  Container::instance().octree->visit( visitor );
  // for each face in box
  for (fp_it j=visitor.mybegin();j!=visitor.myend();++j)
  {
    (*j)->clearFlag();
    // for each vertex of face
    for (int k=0;k<3;++k)
    {
      // store face vertices
      ps.insert((*j)->getVertex(k));
      // OR
      // do not store
      // decide with some degree of accuracy
      // whether vertex has been searched
      // if it has been searched, then continue
      // else perform search
      // OR
      // do store, but first
      // check if frozen
      // check if in full search
      // then decide with some degree of accuracy
      // whether vertex has been searched
      // if it has been searched, then continue
      // else store
    }
  }
  return ps;
}

/** Record all edge angles from adjacent faces to input vertex.
 * \param[in] v Vertex of interest.
 */

void State::recordVertAdjFaceEdgeAngles (Vertex const * const v)
{
  ea.clear();
  // NOTE THIS ALGORITHM LIKELY TO CHECK EDGES TWICE.
  // for each adjacent face of current vertex
  for (fp_cit i=v->begin();i!=v->end();++i)
  {
    ea[(*i)->getEdge(0)] = (*i)->getEdge(0)->getAngle();
    ea[(*i)->getEdge(1)] = (*i)->getEdge(1)->getAngle();
    ea[(*i)->getEdge(2)] = (*i)->getEdge(2)->getAngle();
  }
}

/** Search for new closest point to input vertices.
 *
 * \param[in] v The current vertex.
 * \param[in] fs Hashed set of vertex pointers requiring full closest point search.
 * \param[in] ps Vector of vertex pointers requiring partial closest point search.
 * \return Collection of vertex counts quantifying how often
 * the closest face to vertex changed.
 */

Search_Stats State::updateClosestFaceToVertices (Vertex * const v,
                                                 v_set & fs,
                                                 v_set & ps)
{
  Search_Stats ss = {0,0,0,0};
  Container & c(Container::instance());
  // for each collected vertex requiring full closest point search
  for (vs_it i=fs.begin();i!=fs.end();++i)
  {
    // DEBUG
//    if ((*i)->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
//    {
//      cout << "\nState::updateClosestFaceToVertices: "
//            << "Processing target vertex.\n"; 
//      (*i)->print(cout);
//      //Container::instance().checkClosestFace(1,"FULL");
//    }
    // DEBUG
    // if vertex is not frozen
    if (c.vertexIsFrozen(*i)==false)
    {
      Face * new_face = NULL;
      vector3 p;
      double sqd = 1E30;
      //int neighbor_count = 0;
      bool closest_pt_found = c.findClosestPtToVertex(*i,p,sqd,new_face);
      if (closest_pt_found==true)
      {
        // update statistics
        Face const * old_face = (*i)->getClosestFace();
        // if old and new face are different
        if (old_face!=new_face)
        {
          ss.fs_changed++;
//          // check if vertex of interest adjacent to moved vertex?
//          if ((*i)->faceIsAdjacent(new_face)==true) ss.fs_changed_adj++;
        }
        // update closest face
        (*i)->setFace(new_face);
        //(*i)->setNeighborCount(neighbor_count);
        // DEBUG
//        if ((*i)->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
//        {
//          cout << "\nTARGET VERTEX IN FULL SEARCH\n"; 
//          (*i)->print(cout);
//          cout << "OLD FACE\n"; 
//          old_face->print(cout);
//          cout << "NEW FACE\n"; 
//          new_face->print(cout);
//          cout << "closest_point ["
//                << p.p[0] << " "
//                << p.p[1] << " "
//                << p.p[2] << "], sqd =  "
//                << sqd << endl;
//          cout << "AFTER update face.\n";
//          //Container::instance().checkClosestFace(1,"FULL");
//        }
        // DEBUG
      }
      else
      {
        // update statistics
        Face const * old_face = (*i)->getClosestFace();
        // if old and new face aredifferent
        if (old_face!=NULL) ss.fs_changed++;
        // update closest face
        (*i)->setFace(NULL);
      }
      double gain;
      if (*i==v)
      {
        gain = Gain_Schedule::instance().getGain();
      }
      else
      {
        //gain = Gain_Schedule::instance().getRefGain();
        gain = Gain_Schedule::instance().getMaxGain();
      }
      Virtual_Disp::instance().updateVirtualDisp(*i,closest_pt_found,gain);
    }
  }
  // for each collected vertex not requiring full search
  for (vs_it i=ps.begin();i!=ps.end();++i)
  {
    // if vertex not already processed, i.e. found in fs hashset and not frozen
    if ((binary_search(fs.begin(),fs.end(),*i)==false) && 
        (c.vertexIsFrozen(*i)==false))
    {
      ///// check all adjacent faces of moved vertex to see /////
      ///// if this vertex's closest point now lies on adjacent face /////
      // collect adjacent faces of moved vertex
      // to check as new closest face to this vertex
      // but first make sure each adjacent face
      // is NOT an adjacent face of this vertex
      // THIS IS NOT NECESARY since ps vertices
      // are not fs vertices which included
      // adjacent vertices to moved vertex.
      //
      Face * new_face;      // new closest face to this vertex
      vector3 p;            // new closest point to this vertex
      double old_sqd_dist;  // square of extracellular width for this vertex
      bool new_closest_pt_found = false;
      // record the squared extracellular width to use as measuring stick
      if ((*i)->getClosestFace()!=NULL)
      {
        old_sqd_dist = (*i)->getSqSepDist();
        new_face = (*i)->getClosestFace(); // new closest face to this vertex
      }
      else
      {
        old_sqd_dist = 1E30;
        new_face = NULL;
      }
      // check moved vertex adj faces for new closet face to this vertex
      int dummy;
      double search_radius = Controls::instance().get_update_region_size();
      if (c.findClosestPtToVertexAmongFaces(*i,v->begin(),v->end(),search_radius,
                                            p,old_sqd_dist,new_face,dummy)==true)
      {
        Face const * old_face = (*i)->getClosestFace();
        new_closest_pt_found = true;
        // update statistics
        // if old and new face are different
        if (old_face!=new_face)
        {
          ss.ps_changed++;
          // check if new face is adjacent to vertex
          if ((*i)->faceIsAdjacent(new_face)==true) ss.ps_changed_adj++;
        }
        // update closest face
        (*i)->setFace(new_face);
        // DEBUG
//        if ((*i)->isMatch(TARGET_VERTEX_INDEX,TARGET_VERTEX_NAME)==true)
//        {
//          cout << "\nTARGET VERTEX IN PARTIAL SEARCH\n"; 
//          (*i)->print(cout);
//          cout << "OLD FACE\n"; 
//          if (old_face!=NULL)
//            old_face->print(cout);
//          else
//            cout << "NULL\n";
//          cout << "NEW FACE\n"; 
//          if (new_face!=NULL)
//            new_face->print(cout);
//          else
//            cout << "NULL\n";
//          cout << "closest_point ["
//                << p.p[0] << " "
//                << p.p[1] << " "
//                << p.p[2] << "], sqd =  "
//                << old_sqd_dist << endl;
//          cout << "AFTER update face.\n";
//          //Container::instance().checkClosestFace(1,"PARTIAL");
//        }
        // DEBUG
      }
      // if new closest point was found
      if (new_closest_pt_found==true)
      {
        // update sets with vertex squared virtual displacement
        //Virtual_Disp::instance().setVirtualDisp(*i,(*i)->getSqVirtualDisp(Gain_Schedule::instance().getRefGain()));
        Virtual_Disp::instance().setVirtualDisp(*i,(*i)->getSqVirtualDisp(Gain_Schedule::instance().getGain()));
      }
    }
  }
  return ss;
}

/** Update the location of adjacent faces to moved vertex in octree.
 * \param[in] v Last moved vertex.
 */

void State::updateAdjacentFacesInTree (Vertex const * const v)
{
  Container & c(Container::instance());
  // make an agent (probably a singleton)
  Octree_Agent_Face agent;
  // for each adjacent face of input vertex
  for (fp_cit i=v->begin();i!=v->end();++i)
  {
    // remove and reinsert face in octree
    c.octree->remove( *(*i), agent );
    c.octree->insert( *(*i), agent );
  }
}

/** The virtual displacement of all vertices on adjacent
 * faces of moved vertex will have changed, so update
 * the appropriate maps.
 *
 * \param[in] v The current vertex.
 */

void State::updateVertexVD (Vertex * const v)
{
  // collect adjacent edges to current vertex
  vec_ep e;
  v->getAdjacentEdges(e);
  // collect edges from current vertex adjacent faces
  ae.clear();
  for (fp_cit i=v->begin();i!=v->end();++i)
  {
    // store face edges
    ae.push_back((*i)->getEdge(0));
    ae.push_back((*i)->getEdge(1));
    ae.push_back((*i)->getEdge(2));
  }
  // sort and keep unique 
  sort(ae.begin(),ae.end());
  ae.assign(ae.begin(),unique(ae.begin(),ae.end()));
  // set for storing nearby vertices to update
  v_set nearby;
  // for each collected edge from current vertex adjacent faces
  for (ep_it i=ae.begin();i!=ae.end();++i)
  {
    // if edge is not adjacent to current vertex
    if ( find(e.begin(),e.end(),*i)==e.end())
    {
      // insert all four edge vertices (v1,v2,o1,o2) into set
      nearby.insert((*i)->getO1());
      nearby.insert((*i)->getO2());
    }
  }
  // update collected vertices excluding current vertex
  // since it has already been updated
  // for each collected vertex
  for (vs_it i=nearby.begin();i!=nearby.end();++i)
  {
    // if vertex has a closest point, then update sets
    if ((*i)->getClosestFace()!=NULL && *i!=v)
    {
      //double gain = Gain_Schedule::instance().getRefGain();
      double gain = Gain_Schedule::instance().getGain();
      double sq_vd = (*i)->getSqVirtualDisp(gain);
      Virtual_Disp::instance().setVirtualDisp(*i,sq_vd);
    }
  }
}

/** Attempt move by input vertex.
 * \param[in] v Vertex to move.
 * \param[in] new_pos New location of vertex after move.
 * \param[out] face_intersection True if interecting faces were detected
 * and move was aborted.
 * \param[out] extreme_angle True if extreme angles (near 0 or 2*pi)
 * were detected and move was aborted.
 * \param[out] outside_octree True if vertex move breached octree boundary
 * and move was aborted.
 * \return True if move successful; false otherwise and vertex does not move.
 */

bool State::assignNewVertexCoords (Vertex * const v,
                                   vector3 const * const new_pos,
                                   vector3 const & old_pos,
                                   bool & face_intersection,
                                   bool & extreme_angle,
                                   bool & outside_octree)
{
  Intersecting_Faces & i_f(Intersecting_Faces::instance());
  // store current vertex position
  //double pO[3]={*v->getCoord(0),*v->getCoord(1),*v->getCoord(2)};
  // store adjacent face bounding boxes
  Vector3r adjacent_face_lower[v->getNumAdjFaces()];
  Vector3r adjacent_face_upper[v->getNumAdjFaces()];
  v->recordAdjFaceBoundingBoxes(&adjacent_face_lower[0],&adjacent_face_upper[0]);
  // define local region
  vector3 lower,upper;
  v->defineLocalRegion(lower,upper);
  // DEBUG
//  cout << "\nState::assignNewVertexCoords: "
//        << "local region ["
//        << lower.p[0] << " "
//        << upper.p[0] << " "
//        << lower.p[1] << " "
//        << upper.p[1] << " "
//        << lower.p[2] << " "
//        << upper.p[2] << "]\n";
  // DEBUG
  // collect vertices in local region
  // whose closest point is on adjacent face of current vertex
  v_set fs = getVertsForFullClosestPtSearch(v,lower,upper);
  // collect edge angles
  recordVertAdjFaceEdgeAngles(v);
  // set current position to holding position
  v->setNewPos(new_pos);
  // update adjacent face normals
  v->updateAdjacentFaceNormals();
  updateVertexNormals(fs);
  // update octree
  updateAdjacentFacesInOctree(v,&adjacent_face_lower[0],&adjacent_face_upper[0]);
  // if no faces intersect and no edges have small angles
  face_intersection = (i_f.vertAdjFacesHaveNewInt(v) && 
                       Controls::instance().get_intersection_weight()!=100.0);
  extreme_angle = State::instance().extremeAngleFound();
  outside_octree = Container::instance().vertexOutsideOctreeBounds(new_pos);
  if (face_intersection==false && extreme_angle==false && outside_octree==false)
  {
    // process intersecting faces, if any
    if (i_f.getCountOfIntFaces(false)>0)
    {
      // store vertices for which niceness may have changed
      // i.e. vertices of faces that were intersected
      // but now aren't, and vice versa
      hashset_v ni = i_f.getNiceCheckSet(v);
      // detect niceness changes and collect changed vertices
      i_f.getNiceSet(fs,ni);
    }
    // add to previous collection of vertices the collection 
    // of all vertices in the same local region as before
    v_set ps = getVertsForPartialClosestPtSearch(lower,upper);
    // update closest point and global energy 
    // for affected vertices collected before
    assert(Intersecting_Faces::instance().intFacesAreSymmetric());
    Search_Stats ss = updateClosestFaceToVertices(v,fs,ps);
    // update sets with squared virtual displacement of nearby vertices
    updateVertexVD(v);
    // update global energy due to this vertex
    // recalcEnergyAfterVertMove(v);
    Log::instance().updateClosestPtSearchStats(fs.size(),ps.size(),ss);
    return true;
  }
  else
  {
    // store current bounding boxes of adjacent faces
    Vector3r new_adjacent_face_lower[v->getNumAdjFaces()];
    Vector3r new_adjacent_face_upper[v->getNumAdjFaces()];
    v->recordAdjFaceBoundingBoxes(&new_adjacent_face_lower[0],
                                  &new_adjacent_face_upper[0]);
    // move vertex back
    v->setPos(old_pos.p[0],old_pos.p[1],old_pos.p[2]);
    // update adjacent face normals
    v->updateAdjacentFaceNormals();
    updateVertexNormals(fs);
    // update octree
    updateOctree(v,&new_adjacent_face_lower[0],&new_adjacent_face_upper[0],
                 &adjacent_face_lower[0],&adjacent_face_upper[0]);
    return false;
  }
}

/** Require vertex move to not worsen the condition of extreme angles.
 * \return True if vertex move is an extreme angle violation;
 * false if move passed the extreme angle test.
 */

bool State::extremeAngleFound (void) const
{
  Controls & cs(Controls::instance());
  // for each element in hashtable (edge*->double)
  // Note sequence of elements in hash_map is not necessarily repeatable.
  for (edhm_cit j=ea.begin();j!=ea.end();++j)
  {
    double new_angle=(*j).first->getAngle();
    // if new angle is extreme
    if (fabs(new_angle)<cs.get_edge_angle_threshold() ||
        fabs(2*cs.get_pi()-new_angle)<cs.get_edge_angle_threshold())
    {
      // if old angle is not extreme or angle change is wrong
      // then record violation
      if (   !(fabs((*j).second)<cs.get_edge_angle_threshold()  ||
               fabs(2*cs.get_pi()-(*j).second)<cs.get_edge_angle_threshold()) ||
             angleChangeIsWrong((*j).second,new_angle) )
      {
        return true;
      }
    }
  }
  // no adjacent edges of current vertex violate edge angle threshold
  return false;
}

/** Determine if extreme angle improved.
 * \param[in] old_angle Extreme angle before the last vertex move.
 * \param[in] new_angle Extreme angle after the last vertex move.
 * \return False if the change improved the angle; true otherwise.
 */

bool State::angleChangeIsWrong (float const & old_angle,float const & new_angle) const
{
  if (old_angle < Controls::instance().get_pi())
  {
    // angle should increase towards PI
    // angle increase is correct
    // if angle increases, return false
    // if angle decreases, return true
    if ( new_angle>old_angle ) 	{return false;}
    else 			{return true;} 
  }
  else
  {
    // assume old_angle > PI, i.e. not exactly PI
    // angle should decrease towards PI
    // angle decrease is correct
    // if angle decreases, return false
    // if angle increases, return true
    if ( new_angle<old_angle ) 	{return false;}
    else 			{return true;} 
  }
}

/** Update vertex normals of vertices around moved vertex.
 *
 * \param[in] fs Hashed set of vertex pointers requiring full closest point search.
 */

void State::updateVertexNormals (v_set & fs)
{
  // for each collected vertex requiring full closest point search
  for (vs_it i=fs.begin();i!=fs.end();++i)
  {
    (*i)->setNormal();
  }
}

/** Update the octree with the new position of input adjacent faces to moved vertex.
 * \param[in] v Moved vertex.
 * \param[in] adjacent_face_lower Pointer to lower corner of pre-move bounding box of 
 * first adjacent face to moved vertex.
 * \param[in] adjacent_face_upper Pointer to upper corner of pre-move bounding box of 
 * first adjacent face to moved vertex.
 */

void State::updateAdjacentFacesInOctree (Vertex * const v,
                                         Vector3r * const adjacent_face_lower,
                                         Vector3r * const adjacent_face_upper)
{
  // store current bounding boxes of adjacent faces
  Vector3r new_adjacent_face_lower[v->getNumAdjFaces()];
  Vector3r new_adjacent_face_upper[v->getNumAdjFaces()];
  v->recordAdjFaceBoundingBoxes(&new_adjacent_face_lower[0],
                                &new_adjacent_face_upper[0]);
  // update tree
  updateOctree(v,adjacent_face_lower,adjacent_face_upper,
                 &new_adjacent_face_lower[0],&new_adjacent_face_upper[0]);
}

/** Update the location of moved vertex adjacent faces in the octree.
 * \param[in] v Moved vertex.
 * \param[in] old_adjacent_face_lower Pointer to lower corner of old bounding box of 
 * first adjacent face to moved vertex.
 * \param[in] old_adjacent_face_upper Pointer to upper corner of old bounding box of 
 * first adjacent face to moved vertex.
 * \param[in] new_adjacent_face_lower Pointer to lower corner of new bounding box of 
 * first adjacent face to moved vertex.
 * \param[in] new_adjacent_face_upper Pointer to upper corner of new bounding box of 
 * first adjacent face to moved vertex.
 */

void State::updateOctree (Vertex * const v,
                          Vector3r * const old_adjacent_face_lower,
                          Vector3r * const old_adjacent_face_upper,
                          Vector3r * const new_adjacent_face_lower,
                          Vector3r * const new_adjacent_face_upper)
{
  int j=0;
  // for each adjacent face
  for (fp_cit i=v->begin();i!=v->end();++i)
  {
    Octree_Visitor_Update visitor(*i,old_adjacent_face_lower[j],
                                     old_adjacent_face_upper[j],
                                     new_adjacent_face_lower[j],
                                     new_adjacent_face_upper[j]);
    // execute visitor
    Container::instance().octree->visit( visitor );
    j++;
  }
}
