// Author: Justin Kinney
// Date: Sep 2008

#include "intersecting_faces.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "nice.h"
#include "controls.h"
#include "opttritri.h"
#include "container.h"

using std::cout;
using std::endl;

Intersecting_Faces * Intersecting_Faces::only_one = NULL;

Intersecting_Faces & Intersecting_Faces::instance(void)
{
  // Not thread-safe.
  // -- lock mutex
  if (only_one == NULL)
    only_one = new Intersecting_Faces();
  // -- unlock mutex
  return *only_one;
}

Intersecting_Faces::Intersecting_Faces (void)
:intf()
{
}

/** Return the stored collection of faces that intersect the input Face.
 * \param[in] face Intersected face of interest.
 * \return Pointer to collection of interecting faces of input face.
 */

vec_fp * Intersecting_Faces::getIntersectingFacesRHS (Face const * const face)
{
  htff_it i = intf.find(const_cast<Face*>(face));
  // ensure that face of interest is found
  // in the table of intersected faces
  if (i==intf.end())
  {
    cout << "Intersecting_Faces::getIntersectingFacesRHS: Error: "
          << "Face was not intersected LHS.\n";
    face->print(cout);
    assert(i!=intf.end());
    exit(1);
  }
  return &((*i).second);
}

/** Calculate and return the number of intersecting face pairs in model
 * or only the number of intersecting pairs from faces of the same object.
 * \param[in] detect_self If true then only count intersections
 * between faces of the same object. If false, then count all intersecting face pairs.
 * \return Count of pairs of intersecting faces.
 */

int Intersecting_Faces::getCountOfIntFaces (bool detect_self)
{
  int count = 0;
  // for each element of hash_map
  // NOTE sequence of elements in hash_map is not necessarily repeatable.
  for (htff_it i = intf.begin();i!=intf.end();++i)
  {
    // if this face has intersecting faces
    assert((*i).second.empty()==false);
    // for each intersecting face
    for (fp_it j=(*i).second.begin();j!=(*i).second.end();++j)
    {
      // assert that intersecting face pair (lhs,rhs) and (rhs,lhs) are in map
      assert(faceIsIntersectedLHS(*j)==true);
      if (faceIntersectsFace(*j,(*i).first)==false)
      {
        (*j)->print(cout);
        (*i).first->print(cout);
        assert(faceIntersectsFace(*j,(*i).first)==true);
      }
      if (detect_self==true)
      {
        if ((*i).first->getVertex(0)->getObject()->getName() ==
            (*j)->getVertex(0)->getObject()->getName())
        {
          count++;
        }
      }
      else
      {
        count++;
      }
    }
  }
  assert((count%2) == 0);
  return count/2;
}

/** Determine if an intersection between two faces is recorded in class.
 * \param[in] lhs Face pointer to use as key of intersecting faces table. 
 * \param[in] rhs Putative intersecting face of lhs face. 
 * \return True if rhs is found in collection
 * of intersecting faces for lhs
 * (with note to the directionality of check); otherwise false.
 */

bool Intersecting_Faces::faceIntersectsFace (Face const * const lhs,
                                             Face const * const rhs)
{
  // get this face's intersecting faces vector
  //cout << "AAA";cout.flush();
  vec_fp const * ifv=getIntersectingFacesRHS(lhs);
  // look for face needle in face haystack's intersecting faces vector
  fp_cit i=find((*ifv).begin(),(*ifv).end(),rhs);
  // return true if found
  // return false if no found
  return (i!=(*ifv).end());
}

/** Verify symmetry of intersecting faces. In other words, for every
 * LHS->RHS pair look for matching RHS->LHS pair.
 * \return True if all intersecting face pairs are symmetric; otherwise false.
 */

bool Intersecting_Faces::intFacesAreSymmetric (void)
{
  bool symmetric = true;
  for (htff_it i = intf.begin();i!=intf.end();i++)
  {
    Face * lhs = (*i).first;
    for (fp_it j = (*i).second.begin();j!=(*i).second.end();j++)
    {
      Face * rhs = *j;
      if (faceIntersectsFace(rhs,lhs)==false)
      {
        cout << "Intersecting_Faces::intFacesAreSymmetric: "
              << "Error. Following faces not found as LHS->RHS pair.\n";
        rhs->print(cout);
        lhs->print(cout);
        symmetric = false;
      }
    }
  }
  return symmetric;
}

/** Strongly determine if input face is recorded in class as being intersected.
 * \param[in] face Face of interest.
 * \return True if face is used as key in intersecting faces container
 * and has faces stored as intersecting.; otherwise false.
 */

bool Intersecting_Faces::faceIsIntersectedRHS (Face const * const face) const
{
  //return intf[const_cast<Face*>(f)].empty()==false;
  htff_cit i = intf.find(const_cast<Face*>(face));
  if (i==intf.end())
  {
    cout << "Intersecting_Faces::faceIsIntersectedRHS: Error"
          << "Face was not intersected LHS.\n";
    assert(i!=intf.end());
    exit(1);
  }
  else
  {
    return ((*i).second).empty()==false;
  }
}

/** Search for new face intersections with input face.
 * \param[in] face Face of interest.
 * \param[out] int_faces Collection of faces that intersect face.
 * \return True if new face intersections are detected; false otherwise.
 */

bool Intersecting_Faces::detectNewFaceInts (Face * const face,vec_fp & int_faces)
{
  bool int_found = false;
  // make a visitor
  Vector3r lower,upper;
  face->getBoundingBox(lower,upper);
  Octree_Visitor_Face visitor(lower,upper);
  // execute visitor
  Container::instance().octree->visit( visitor );
  // for each unique face
  for (fp_cit i=visitor.mybegin();i!=visitor.myend();++i)
  {
    (*i)->clearFlag();
    if (int_found==true) continue;
    // if unique face is NOT same as input face
    if (*i!=face)
    {
      // if faces intersect
      if (checkFaceFaceInts(face,*i)==true)
      { 
        int_faces.push_back(*i);
        if (Controls::instance().get_strict_face_intersection_prevention()==true)
        {
          // if this face previously had no intersections, then reject move
          if (faceIsIntersectedLHS(face)==false) int_found=true;
          // else face was previously intersected
          else
          {
            // grab previous intersecting faces of this face
  //cout << "BBB";cout.flush();
            vec_fp *fv=getIntersectingFacesRHS(face);
            bool same=false;
            // for each previous intersecting face
            for (fp_it p=(*fv).begin();p!=(*fv).end();++p)
            {
              // if current intersecting face is of same object
              //  as previous intersecting face
              Object const * const o1 = (*p)->getVertex(0)->getObject();
              Object const * const o2 = (*i)->getVertex(0)->getObject();
              if (o1 == o2) same = true;
            }
            // if intersecting objects have changed from before to now
            // then reject move
            if (same==false) int_found=true;
          }
        }
        else
        {
          // if current face if of same object
          // as intersecting face, then reject move
          Object const * const o1 = face->getVertex(0)->getObject();
          Object const * const o2 = (*i)->getVertex(0)->getObject();
          if (o1==o2) int_found=true;
        }
      }
    }
  }
  //Container::instance().checkFaces("Intersecting_Faces::detectNewFaceInts");
  return int_found;
}

/** Find and record all intersecting faces of input face. 
 *
 * \param[in] face Face of interest.
 * \return True if current face is intersected, false otherwise.
 */

bool Intersecting_Faces::findAndRecordNewFaceInt (Face * const face)
{
  bool flag = false;
  // make a visitor
  Vector3r lower,upper;
  face->getBoundingBox(lower,upper);
  Octree_Visitor_Face visitor(lower,upper);
  // execute visitor
  Container::instance().octree->visit( visitor );
  //visitor.sort();
  // for each unique face
  for (fp_it i=visitor.mybegin();i!=visitor.myend();++i)
  {
    (*i)->clearFlag();
    // if unique face is NOT same as input face
    if (*i != face)
    {
      // if faces intersect
      if (checkFaceFaceInts(face,*i)==true)
      {
        if (Controls::instance().get_disable_messages()==false)
        {
          cout << "\nIntersecting_Faces::findAndRecordNewFaceInt :"
                << "face intersection found.\n";
          face->print(cout);
          (*i)->print(cout);
          cout << endl;
        }
        // DEBUG
        //if ( ( face->isMatch(31843,"Pre") && (*i)->isMatch(31851,"Pre") ) ||
        //     ( face->isMatch(31851,"Pre") && (*i)->isMatch(31843,"Pre") ) )
//        if ( ( face->isMatch(204423,"Pre") && (*i)->isMatch(204396,"Pre") ) ||
//             ( face->isMatch(204396,"Pre") && (*i)->isMatch(204423,"Pre") ) )
//        {
//          cout << "\n\nTarget faces intersect.\n";
//          exit(1);
//        }
        // DEBUG
        // save intersecting face* to this face's intersecting face vector
        addFaceToFace(face,*i);
        // return
        flag = true;
      }
    }
  }
  //Container::instance().checkFaces("Intersecting_Faces::findAndRecordNewFaceInt");
  return flag;
}

/** Check if any adjacent face of current vertex is intersected. 
 *
 * \param[in] v The current vertex.
 * \return True if any adjacent face of vertex is intersected, false otherwise.
 */

bool Intersecting_Faces::vertAdjFacesHaveNewInt (Vertex const * const v)
{
  // for each adjacent face of current vertex
  for (fp_cit i=v->begin();i!=v->end();++i)
  {
    vec_fp dummy;
    //  if face intersects any other face, then return true
    if (detectNewFaceInts(*i,dummy)==true){return true;}
  }
  // no adjacent faces of current vertex intersect any other face
  return false;
}

/** Record an intersection between two faces. 
 *
 * \param[in] lhs Face pointer to use as key of intersecting faces table. 
 * \param[in] rhs Intersecting face of lhs face. 
 */

void Intersecting_Faces::addFaceToFace (Face * const lhs,
                                        Face const * const rhs)
{
  // if this face not in intf
  if (faceIsIntersectedLHS(lhs)==false)
  {
    // create new face vector in table
    vec_fp nv;
    intf[lhs]=nv;
  }
  // face f not already in vector
  if (faceIntersectsFace(lhs,rhs)==false)
  {
    // add face f to this face's intersecting faces vector
    intf[lhs].push_back(const_cast<Face*>(rhs));
  }
}

/** Remove records of intersections for input face
 * (with note that search is not exhaustive).
 *
 * \param[in] face Face of interest.
 */

void Intersecting_Faces::setFaceNotIntersectedLHS (Face const * const face)
{
  // if this face is in intf table 
  if (faceIsIntersectedLHS(face)==true)
  {
    // get this face's intersecting faces vector
  //cout << "CCC";cout.flush();
    vec_fp *ifv=getIntersectingFacesRHS(face);
    // if intersecting faces vector is not empty
    if ((*ifv).empty()==false)
    {
      // for each face in intersection vector
      for (fp_it i=(*ifv).begin();i!=(*ifv).end();++i)
      {
        // remove this face from intersection vector of that face
        removeFaceFromFaceInt(*i,face);
      }
      (*ifv).clear();
    }
    // remove element from table
    intf.erase(const_cast<Face*>(face));
  }
}

/** Remove record of intersection between input faces.
 *
 * \param[in] lhs Face pointer to use as key of intersecting faces table. 
 * \param[in] rhs Intersecting face of lhs face that is to be removed. 
 */

void Intersecting_Faces::removeFaceFromFaceInt (Face * const lhs,
                                                Face const * const rhs)
{
  // if this face has intersecting faces in hashmap intf
  if (faceIsIntersectedLHS(lhs)==true)
  {
    // get this face's intersecting faces vector
  //cout << "DDD";cout.flush();
    vec_fp *ifv=getIntersectingFacesRHS(lhs);
    // look for face rhs in face lhs's intersecting faces vector
    fp_it i=find((*ifv).begin(),(*ifv).end(),rhs);
    // if found
    if (i!=(*ifv).end())
    {
      // remove face f from this face's intersecting faces vector
      (*ifv).erase(remove((*ifv).begin(),(*ifv).end(),rhs),(*ifv).end());
      // if intersecting face vector is now empty
      if (faceIsIntersectedRHS(lhs)==false)
      {
        // then remove intersecting face vector from hashtable
        setFaceNotIntersectedLHS(lhs);
      }
    }
  }
}

/** Remove intersections of adjacent faces to current vertex from
 * intersection list, and add all vertices of intersecting faces
 * to set of vertices whose niceness may have changed.
 *
 * \param[in] v The current vertex.
 * \param[in] vertices Collection of vertices whose niceness may have changed.
 */

void Intersecting_Faces::removeOldIntersections (Vertex const * const v,
                                                 hashset_v & vertices)
{
  // collect the vertices of all intersecting faces
  // of all adjacent faces to current vertex
  // into vset, i.e. hash_set of Vertex*s.
  //
  // for every adjacent face
  for (fp_cit k=v->begin();k!=v->end();++k)
  {
    // if adjacent face has intersecting faces
    if (faceIsIntersectedLHS(*k)==true)
    {
      // for each intersecting face of adjacent face
  //cout << "EEE";cout.flush();
      vec_fp *fv=getIntersectingFacesRHS(*k);
      for (fp_it p=(*fv).begin();p!=(*fv).end();++p)
      {
        // add intersecting face vertices to vset2
        vertices.insert((*p)->getVertex(0));
        vertices.insert((*p)->getVertex(1));
        vertices.insert((*p)->getVertex(2));
        // remove adjacent face from intersecting face's vector
        removeFaceFromFaceInt(*p,*k);
      }
    }
    // clear adjacent face's intersecting face vector
    setFaceNotIntersectedLHS(*k);
  }
}

/** Check for face intersections of each adjacent face of current vertex,
 * and add all vertices of intersecting faces
 * to set of vertices whose niceness may have changed.
 *
 * \param[in] v The current vertex.
 * \param[in] vertices Collection of vertices whose niceness may have changed.
 */

void Intersecting_Faces::updateNewIntersections (Vertex const * const v,
                                                 hashset_v & vertices)
{
  vec_fp int_faces;
  // for every adjacent face
  for (fp_cit k=v->begin();k!=v->end();++k)
  {
    int_faces.clear();
    // if adjacent face is currently intersected
    if (detectNewFaceInts(*k,int_faces))
    {
      // if intersected face was previously not intersected
      if (faceIsIntersectedLHS(*k)==false)
      {
        // create new face vector in table
        vec_fp nv;
        intf[*k]=nv;
      }
      // for each intersecting face of adjacent face
      //cout << "FFF";cout.flush();
      // PROBLEM OCCURS HERE
      // Suspect problem is intersecting faces were
      // never added to hash table
      //vec_fp *fv=getIntersectingFacesRHS(*k);
      //for (fp_it p=(*fv).begin();p!=(*fv).end();++p)
      for (fp_it p=int_faces.begin();p!=int_faces.end();++p)
      {
        // add intersecting face vertices to vset2
        vertices.insert((*p)->getVertex(0));
        vertices.insert((*p)->getVertex(1));
        vertices.insert((*p)->getVertex(2));
        // if intersected face was previously not intersected
        if (faceIsIntersectedLHS(*p)==false)
        {
          // create new face vector in table
          vec_fp nv;
          intf[*p]=nv;
        }
        // add adjacent face* to intersected face
        addFaceToFace(*p,*k);
        addFaceToFace(*k,*p);
      }
    }
  }
}

/** Collect vertices whose niceness may have changed after vertex move
 * due to face intersections.
 *
 * \param[in] v The current vertex.
 * \return Collection of vertices whose niceness may have changed.
 */

hashset_v Intersecting_Faces::getNiceCheckSet (Vertex const * const v)
{
  hashset_v ni;
  // add current vertex
  ni.insert(const_cast<Vertex*>(v));
  removeOldIntersections(v,ni);
  updateNewIntersections(v,ni);
  return ni;
}

/** Check niceness of each collected vertex, and, if niceness changed,
 * then add vertex to set of verticess requiring a full closest point search. 
 *
 * \param[out] full_search_faces The collection of vertices
 * requiring full search.
 * \param[in] vertices Collection of vertices whose niceness may have changed.
 */

void Intersecting_Faces::getNiceSet (v_set & full_search_faces,
                                     hashset_v & vertices)
{
  // for each vertex whose niceness may have changed
  for (hv_it i=vertices.begin();i!=vertices.end();++i)
  {
    // update niceness
    if (Nice::instance().updateVertexNiceness(*i))
    {
      // if niceness changed
      // i.e. was nice and now not,
      // or was not nice and now is
      // then add vertex to set 
      // requiring full search for closest point
      full_search_faces.insert(*i);
    }
  }
}

/** Find and record all face intersections.
*/

void Intersecting_Faces::findAllFaceIntersections (void)
{
  cout << "Find all face intersections....................";
  cout.flush();
  Container & c(Container::instance());
  int k=0;
  double goal = 0.2;
  double inc = 0.2;
  double a = 1.0/c.o.size();
  if (goal<a)
  {
    goal=a;
    inc = a;
  }
  cout << "0%..";
  cout.flush();
  // clear table
  intf.clear();
  // for each object
  for (o_it i=c.o.begin();i!=c.o.end();++i)
  {
    // for each face in object
    for (f_it j=i->f.begin();j!=i->f.end();++j)
    {
      findAndRecordNewFaceInt(&(*j));
    }
    // track progress
    double progress = static_cast<double>(k++)*a;
    if (progress>goal || !distinguishable(progress,goal))
    {
      cout << static_cast<int>(goal*100) << "%..";
      cout.flush();
      goal+=inc;
    }
  }
  //cout << "100%..complete.\n";
  cout << "complete.\n";
  cout.flush();
}

/** Calculate total force vector on current face
 * due to intersection with other faces.
 *
 * \param[in] face The current face.
 * \param[out] total_force Intersection force is added to total_force.
 */

void Intersecting_Faces::getFaceIntersectionForce (Face * const face,
                                                   vector3 & total_force)
{
  vector3 d;
  // if this face has intersecting faces
  assert(faceIsIntersectedLHS(face)==true && faceIsIntersectedRHS(face)==true);
  // get intersecting faces
  //cout << "GGG";cout.flush();
  vec_fp* ptr = getIntersectingFacesRHS(face);
  // for each intersecting face
  for (fp_it i=ptr->begin();i!=ptr->end();++i)
  {
    getFaceFaceIntForce(face,*i,d);
  }
  for (int i=0;i<3;++i) {total_force.p[0]+=d.p[0];}
}

/** Calculate force vector on lhs face due to intersection with rhs face.
 *
 * \param[in] lhs The current face.
 * \param[in] rhs The other face.
 * \param[out] total_force Intersection force is added to total_force.
 */

void Intersecting_Faces::getFaceFaceIntForce (Face * const lhs,
                                              Face * const rhs,
                                              vector3 & total_force)
{
  assert(faceIntersectsFace(lhs,rhs) && faceIntersectsFace(rhs,lhs));
  // get current face unit normal
  vector3 const * n1 = lhs->getNormal();
  double m1=1.0/sqrt(n1->dot(*n1));
  vector3 nn1((*n1).p[0]*m1,(*n1).p[1]*m1,(*n1).p[2]*m1);
  // get intersecting face unit normal
  vector3 const * n2 = rhs->getNormal();
  double m2=1.0/sqrt(n2->dot(*n2));
  vector3 nn2((*n2).p[0]*m2,(*n2).p[1]*m2,(*n2).p[2]*m2);
  // compute resultant vector
  vector3 R(nn1+nn2);
  //
  vector3 T(nn1.cross(nn2));
  vector3 P(T.cross(R));
  // compute force
  P *= Controls::instance().get_intersection_weight()/sqrt(P.dot(P));
  total_force += P;
}

/** For a pair of faces, each edge of each face is checked
 * for intersection with other face.
 *
 * \param[in] cf The current face.
 * \param[in] of The other face.
 * \return True if faces intersect, otherwise false.
 */

bool Intersecting_Faces::checkFaceFaceInts (Face const * const cf,
                                            Face const * const of) const
{
  // get number of unique vertices between current and other face
  // single_shared_vert[0]==cf index of single shared vertex
  // single_shared_vert[1]==of index of single shared vertex
  // where index == 0,1,2
  int single_shared_vert[2]={-1,-1};
  int num_unique = getNumUniqueVerts(cf,of,single_shared_vert);
  // DEBUG
//  bool myflag = false;
//  //if ( ( cf->isMatch(31843,"Pre") && of->isMatch(31851,"Pre") ) ||
//  //     ( cf->isMatch(31851,"Pre") && of->isMatch(31843,"Pre") ) )
//  //if ( ( cf->isMatch(201480,"Pre") && of->isMatch(231208,"Pre") ) ||
//  //     ( cf->isMatch(231208,"Pre") && of->isMatch(201480,"Pre") ) )
//  if ( ( cf->isMatch(204396,"Pre") && of->isMatch(204423,"Pre") ) ||
//       ( cf->isMatch(204423,"Pre") && of->isMatch(204396,"Pre") ) )
//  {
//    myflag = true;
//    cout << "\n\nIntersecting_Faces::checkFaceFaceInts: "
//          << "num_unique = " << num_unique << endl;
//  }
//  // DEBUG
  // if polygons are identical
  if (num_unique==3) return true;
  // if shared edge
  else if (num_unique==4)
  {
    // if not parallel, i.e. NOT coplanar
    if (facesParallel(cf,of)==false) return false;
    // faces are parallel, i.e. coplanar
    // check polygon edges intersection
    return checkEdgeEdgeIntersection(cf,of);
  }
  // if single shared vertex 
  else if (num_unique==5)
  {
    // if parallel, i.e. coplanar
    if (facesParallel(cf,of)==true)
    {
      // check polygon edges intersection
      return checkEdgeEdgeIntersection(cf,of);
    }
    // faces are NOT parallel, i.e. NOT coplanar
    vector3 const * cv0, * cv1, * cv2;
    vector3 const * ov0, * ov1, * ov2;
    cf->getVertexCoord(cv0,cv1,cv2);
    of->getVertexCoord(ov0,ov1,ov2);
    vector3 const * target1,* target2;
    if      (single_shared_vert[0]==0) { target1 = cv1; target2 = cv2; }
    else if (single_shared_vert[0]==1) { target1 = cv0; target2 = cv2; }
    else                               { target1 = cv0; target2 = cv1; }
    // check intersection of current face edge with other face
    result r;
    if (intersect_triangle3(target1,target2,of->getNormal(),ov0,ov1,ov2,r)==true) return true;
    // check intersection of other face edge with current face
    if      (single_shared_vert[1]==0) { target1 = ov1; target2 = ov2; }
    else if (single_shared_vert[1]==1) { target1 = ov0; target2 = ov2; }
    else                               { target1 = ov0; target2 = ov1; }
    if (intersect_triangle3(target1,target2, cf->getNormal(),cv0,cv1,cv2,r)==true) return true;
    // no face intersection
    return false;
  }
  // no shared vertices
  else
  {
    // if coplanar
    if (facesCoplanar(cf,of)==true)
    {
      // check polygon edges intersection
      return checkEdgeEdgeIntersection(cf,of);
    }
    // faces are NOT parallel, i.e. NOT coplanar
    // get face vertex coordinates
    vector3 const * cv0, * cv1, * cv2;
    vector3 const * ov0, * ov1, * ov2;
    cf->getVertexCoord(cv0,cv1,cv2);
    of->getVertexCoord(ov0,ov1,ov2);
    return NoDivTriTriIsect(cv0,cv1,cv2,ov0,ov1,ov2);
  }
}

/** Count number of unique vertices for pair of faces,
 * and, if exactly one shared vertex, then identify
 * shared vertex in each face.
 *
 * \param[in] cf The current face.
 * \param[in] of The other face.
 * \param[out] single_shared If exactly one shared vertex,
 * then save it's identity in each face.
 * \return Number of unique vertices for pair of faces.
 */

int Intersecting_Faces::getNumUniqueVerts (Face const * const cf,
                                           Face const * const of,
                                           int * const single_shared) const
{
  int num_unique = 0;
  // for each current face vertex
  for (int j=0;j<3;++j)
  {
    // for each other face vertex
    for (int k=0;k<3;++k)
    {
      if (cf->getVertex(j)!=of->getVertex(k))
      {
        num_unique++;
      }
      else
      {
        single_shared[0]=j;
        single_shared[1]=k;
      }
    }
  }
  return num_unique-3;
}

/** Check if pair of faces are parallel.
 *
 * \param[in] cf The current face.
 * \param[in] of The other face.
 * \return True if faces are parallel, otherwise false.
 */

bool Intersecting_Faces::facesParallel (Face const * const cf,
                                        Face const * const of) const
{
  // get face normals
  vector3 const * cn = cf->getNormal();
  vector3 const * on = of->getNormal();
  // are current face and other face parallel
  // i.e. is angle between normals equal to zero?
  // i.e. is the square of the cosine of the angle equal to 1?
  if (!distinguishable(cn->dot(*on)*cn->dot(*on),cn->dot(*cn)*on->dot(*on),1e-5)) return true;
  else return false;
}

/** Check if pair of faces are coplanar.
 *
 * \param[in] cf The current face.
 * \param[in] of The other face.
 * \return True if faces are coplanar, otherwise false.
 */

bool Intersecting_Faces::facesCoplanar (Face const * const cf,
                                        Face const * const of) const
{
  if (facesParallel(cf,of)==true)
  {
    // get one face normal
    vector3 const * cn = cf->getNormal();
    // get vector between any two vertices of faces
    Vertex const * va = cf->getVertex(0);
    int i = 0;
    while (of->getVertex(i)->isMatch(va->getIndex(),va->getObject()->getName())==true)
    {
      i++;
      assert(i<3);
    }
    Vertex const * vb = of->getVertex(i);
    vector3 const * pa = va->getPos();
    vector3 const * pb = vb->getPos();
    vector3 diff = *pa-*pb;
    // DEBUG
    //if (flag)
    //{
    //  cout << "Intersecting_Faces::facesCoplanar: "
    //        << "pa = ";
    //  pa->print(cout);
    //  cout << "\nIntersecting_Faces::facesCoplanar: "
    //        << "pb = ";
    //  pb->print(cout);
    //  cout << "\nIntersecting_Faces::facesCoplanar: "
    //        << "diff = ";
    //  diff.print(cout);
    //  cout << "\n";
    //  cout << "Intersecting_Faces::facesCoplanar: "
    //        << "dot = " << fabs(cn->dot(diff)) << "\n";
    //}
    // DEBUG
    // are current face and other face coplanar
    // i.e. is angle between one face normal and
    //      vector between any two vertices of faces equal to zero?
    // i.e. is the square of the cosine of the angle equal to 1?
    if (fabs(cn->dot(diff)) < Controls::instance().get_epsilon()) return true;
    else return false;
  }
  else
  {
    return false;
  }
}

/** For a pair of faces, each edge of each face is checked
 * for intersection with an edge of the other face.
 *
 * \param[in] cf The current face.
 * \param[in] of The other face.
 * \return True if faces intersect, otherwise false.
 */

bool Intersecting_Faces::checkEdgeEdgeIntersection (Face const * const cf,
                                                    Face const * const of) const
{
  Controls & cs(Controls::instance());
  // algorithm from 
  // David Eberly
  // Geometric Tools, LLC
  // http://www.geometrictools.com/
  int pairs[3][2] = {{0,1},{1,2},{2,0}};
  // cpvc = current_polygon_vertex_coordinates
  // opvc = other_polygon_vertex_coordinates
  // cv   = current_vertex
  // ov   = other_vertex
  Vertex *cv[2],*ov[2];
  // for each current face edge
  for (int i=0;i<3;++i)
  {
    // for each other face edge
    for (int j=0;j<3;++j)
    {
      cv[0] = cf->getVertex(pairs[i][0]);
      cv[1] = cf->getVertex(pairs[i][1]);
      ov[0] = of->getVertex(pairs[j][0]);
      ov[1] = of->getVertex(pairs[j][0]);
      // if the edges do not share a vertex
      if (cv[0]!=ov[0]&&cv[0]!=ov[1]&&cv[1]!=ov[0]&&cv[1]!=ov[1])
      {
        vector3 del(*(ov[0]->getPos())-*(cv[0]->getPos()));
        vector3 e0(*(cv[1]->getPos())-*(cv[0]->getPos()));
        vector3 e1(*(ov[1]->getPos())-*(ov[0]->getPos()));
        vector3 e0p(e0.p[1],-e0.p[0],0);
        vector3 e1p(e1.p[1],-e1.p[0],0);
        // and the edges are not parallel
        if (e0.dot(e1p)!=0.0)
        {
          // compute scalars
          double t = del.dot(e0p);
          double s = del.dot(e1p);
          if ( (t > cs.get_my_double_epsilon()) && 
               (t < (1.0-cs.get_my_double_epsilon())) && 
               (s > cs.get_my_double_epsilon()) && 
               (s < (1.0-cs.get_my_double_epsilon())) )
          {
            return true;
          }
        }
      }
    }
  }
  return false;
}

