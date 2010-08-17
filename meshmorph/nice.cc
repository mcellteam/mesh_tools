// Author: Justin Kinney
// Date: Sep 2008

#include "nice.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "controls.h"
#include "container.h"
#include "opttritri.h"

using std::cout;
using std::endl;

Nice * Nice::only_one = NULL;

Nice & Nice::instance()
{
  // Not thread-safe.
  // -- lock mutex
  if (only_one == NULL)
    only_one = new Nice();
  // -- unlock mutex
  return *only_one;
}

Nice::Nice(void)
  :nonnice()
{
}

/** Check if vertex is recorded as being nice.
 * \param[in] v Vertex of interest.
 * \return True if vertex is nice; false otherwise.
 */

bool Nice::vertexIsNice (Vertex const * const v) const
{
  return nonnice.find(const_cast<Vertex*>(v))==nonnice.end();	
}

/** Set nice value of vertex with code described in nice.cc.
 * \param[in] newval Set niceval of input vertex to new value.
 * \param[in] v Vertex of interest.
 */

void Nice::setVertexNiceVal (int const & newval,Vertex * const v)
{
  // newval == 0, nice
  // newval == 1, nonnice in different object
  // newval == 2, nonnice in same object
  if (newval==0)
  {
    vhm_it i = nonnice.find(v);
    if (i!=nonnice.end())
    {
      nonnice.erase(i);
    }
  }
  else
  {
      nonnice[v]=newval;
  }
}

/** Get nice value of vertex with code described in nice.cc.
 * \param[in] v Vertex of interest.
 * \return Nice value of vertex.
 */

int Nice::getVertexNiceVal (Vertex const * const v) const
{
  if (vertexIsNice(const_cast<Vertex*>(v)))
  {
    return 0;
  }
  else
  {
     vhm_cit i = nonnice.find(const_cast<Vertex*>(v));	
     return (*i).second;
  }
}

/** Identify a location near the vertex of interest
 * guaranteed to be outside of the parent object of vertex.
 * \param[in] v Vertex of interest.
 * \param[out] extra_obj_pt Location outside of vertex parent object,
 * if successfully found.
 * \param[in] adjacent_face Iterator pointing to vertex adjacent face
 * to use as parent of extracellular location.
 * \return 1 if extra-object location identified; 0 otherwise.
 */

bool Nice::findExtraPoint (Vertex const * const v,
                          vector3 & extra_obj_pt,
                          fp_cit const & adjacent_face) const
{
  // calculate ray from adjacent face centroid
  // along face normal small distance
  vector3 face_centroid;
  getVertAdjFaceRay(face_centroid,extra_obj_pt,adjacent_face);
  // find all intersected faces along ray
  vec_fp crossed_faces,edge_faces;
  // find intersected faces along ray
  face_grp fg = findIntFacesAlongRay(v,face_centroid,extra_obj_pt,false);
  sort(fg.crossed_faces.begin(),fg.crossed_faces.end());
  // prohibit edge intersections, since they are tricky and avoidable
  if (fg.edge_faces.empty()==false) return false;
  // from intersected faces, determine 
  vec_op odd_objects;
  getPenetratedObjs(fg.crossed_faces,fg.edge_faces,odd_objects);
  // early exit if no odd objects
  if (odd_objects.empty()==true) return true;
  // odd objects found
  std::pair<op_it,op_it> i;
  i=equal_range(odd_objects.begin(),odd_objects.end(),v->getObject());
  if (i.first!=i.second)
  {
    // vertex parent object was one of odd meshes
    return false;
  }
  // vertex parent object was NOT one of odd meshes
  return true;
}

/** Find and record all nonnice vertices in model.
 */

void Nice::findNonniceVertices (void)
{
  cout << "Find nice vertices.............................";
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
  printf("0%%..");
  fflush(stdout);
  // for each object in container
  for (o_it i=c.o.begin();i!=c.o.end();++i)
  {
    // for each vertex in object
    for (v_it j=i->v.begin();j!=i->v.end();++j)
    {
      updateVertexNiceness(&(*j));
    }
    // track progress
    double progress = static_cast<double>(k++)*a;
    if (progress>=goal)
    {
      printf("%d%%..",static_cast<int>(goal*100));
      fflush(stdout);
      goal+=inc;
    }
  }
//  printf("100%%..");
//  fflush(stdout);
  cout << "complete.\n";
  cout.flush();
}

/** Evaluate niceness of vertex with ray tracing
 * and record and return result.
 * \param[in] v Vertex of interest.
 * \return True if vertex niceness changed since last check;
 * false otherwise.
 */

bool Nice::updateVertexNiceness (Vertex * const v)
{
  vec_op all_crossed_objects;
  // collect objects inside which vertex lies
  getCrossedObjects(v,all_crossed_objects);
  // update niceness of vertex based on all_crossed_objects
  return setVertexNiceness(v,all_crossed_objects);
}

/** Identify a location guaranteed to be outside 
 * the parent object of the vertex of interest
 * and identify object crossings
 * between point and current vertex.
 * \param[in] v Vertex of interest.
 * \param[out] extra_obj_pt Location guaranteed to be
 * outside vertex parent object.
 * \param[out] crossed_objects Collection of crossed objects
 * between vertex and extra-object location.
 * \return True if extra-object location found; false otherwise.
 */

bool Nice::getPointOutsideObject(Vertex const * const v,
                                 vector3 & extra_obj_pt,
                                 vec_op & crossed_objects) const
{
  // find location (extra_obj_pt) outside of current object
  bool edge_intersected = true;
  fp_cit adjacent_face= v->begin();
  do {
    // if all adjacent faces of vertex failed
    if (adjacent_face==v->end())
    {
      cout << "\n\nNice::collectCrossed: "
            << "Intersected edge along path from vertex "
            << "to extra-object location on every adjacent face.\n"
            << "Warning! Setting undetermined vertex to 'nice'.\n";
      v->print(std::cout);
      crossed_objects.clear();
      return true;
    }
    // while searching for location outisde of current object
    while (findExtraPoint(v,extra_obj_pt,adjacent_face)==false)
    {
      adjacent_face++;
      if (adjacent_face==v->end())
      {
        cout << "\n\nNice::collectCrossed: "
              << "Failed to find valid extracellular point on any adjacent face.\n"
              << "Warning! Setting undetermined vertex to 'nice'.\n";
        v->print(std::cout);
        crossed_objects.clear();
        return false;
      }
    }
    // grab intersected objects between vertex and extra-object point
    edge_intersected = getCrossedObjFromVertToExtra(v,extra_obj_pt,crossed_objects);
    // if a face edge was intersected, then search again with different ray
    if (edge_intersected==true) adjacent_face++;
  }
  while (edge_intersected==true);
  return true;
}

/** Collect and return objects inside which vertex lies
 * as determined by ray tracing.
 * \param[in] v Vertex of interest.
 * \param[out] all_crossed_objects Collection of objects inside which vertex lies.
 */

void Nice::getCrossedObjects (Vertex const * const v,
                              vec_op & all_crossed_objects) const
{
  // find location, extra_object_pt, outside of vertex parent object
  // and grab intersected objects between vertex and extra-object point
  // and return as crossed_objects
  vector3 extra_object_pt;
  vec_op crossed_objects;
  bool a = getPointOutsideObject(v,extra_object_pt,crossed_objects);
  // if no extra object location found then return
  if (a==false) return;
  // get 6 rays along principal axes from extra-object location
  // sorted from shortest to longest on index 0 to 5.
  double rays[6][3];
  getRaysToWorldLimit(extra_object_pt,rays); // returns ray
  // do while an axis needs to be found and an axis is unchecked.
  bool edge_intersected = true;
  vector3 world_limit;
  int select_axis = 0;
  do {
    world_limit.p[0]=rays[select_axis][0];
    world_limit.p[1]=rays[select_axis][1];
    world_limit.p[2]=rays[select_axis][2];
    edge_intersected = getCrossedObjFromExtraToLimit(extra_object_pt,
                                                     world_limit,
                                                     all_crossed_objects);
    if (edge_intersected==true) select_axis++;
  } while (edge_intersected==true && select_axis<6);
  // error checking
  if (select_axis==6)
  {
    cout << "\n\nNice::collectCrossed: ERROR: All six principal coordinates"
          << " resulted in face edge collisions.\n";
    v->print(std::cout);
    cout << endl;
    assert(select_axis!=6);
    exit(1);
  }
  // merge object collections from vertex to extra-object point
  // to object collection from extra-object location to world limit
  sort(all_crossed_objects.begin(),all_crossed_objects.end());
  for (op_it i=crossed_objects.begin();i!=crossed_objects.end();++i)
  {
    std::pair<op_it,op_it> j;
    j=equal_range(all_crossed_objects.begin(),all_crossed_objects.end(),*i);
    // if object in crossed_objects is found in all_crossed_objects
    // then remove from all_crossed_objects
    if (j.first!=j.second){all_crossed_objects.erase(j.first);}
    // else add it
    else
    {
      all_crossed_objects.push_back(*i);
    }
  }
}

/** Record new nice value based on collection of objects
 * crossed in ray trace from vertex to world limits.
 * \param[in] v Vertex of interest.
 * \param[in] crossed_objects Collection of objects crossed in ray trace
 * from vertex to world limits.
 * \return True if vertex niceness changed since last check; false otherwise.
 */

bool Nice::setVertexNiceness (Vertex * const v,vec_op & crossed_objects)
{
  // new_nice==0 means nice
  // new_nice==1 means nonnice in different object
  // new_nice==2 means nonnice in same object
  int old_nice = getVertexNiceVal(v);
  int new_nice = 0;
  // if crossed_objects is not empty
  if (crossed_objects.empty()==false)
  {
    // vertex is not nice
    new_nice++;
    // if vertex is inside self object
    std::pair<op_it,op_it> i;
    i=equal_range(crossed_objects.begin(),crossed_objects.end(),v->getObject());
    if (i.first!=i.second){new_nice++;}
  }
  setVertexNiceVal(new_nice,v);
  // if vertex niceness changes then return true, else return false
  return old_nice!=new_nice;
}

/* Calculate rays along each principal axis
 * from location to world limits.
 * \param[in] point Location.
 * \param[out] rays End points of rays to world limits
 * sorted from shortest to longest path. 
 */

void Nice::getRaysToWorldLimit (vector3 const & point,
                                double rays[6][3]) const
{
  Container  & c(Container::instance());
  // map: (distance to world bounds)double->(direction code)int
  map_di dist;
  dist.insert(std::make_pair(fabs(point.p[0]-c.getWorld(0)),0));
  dist.insert(std::make_pair(fabs(point.p[0]-c.getWorld(1)),1));
  dist.insert(std::make_pair(fabs(point.p[1]-c.getWorld(2)),2));
  dist.insert(std::make_pair(fabs(point.p[1]-c.getWorld(3)),3));
  dist.insert(std::make_pair(fabs(point.p[2]-c.getWorld(4)),4));
  dist.insert(std::make_pair(fabs(point.p[2]-c.getWorld(5)),5));
  int j=0;
  // for each map element
  for (di_it i=dist.begin();i!=dist.end();++i)
  {
    rays[j][0] = point.p[0];
    rays[j][1] = point.p[1];
    rays[j][2] = point.p[2];
    int index=(*i).second;
    if		(index==0)
    {
      rays[j][0] = point.p[0]-2*(c.getWorld(1)-c.getWorld(0));	// end x
    }
    else if (index==1)
    {
      rays[j][0] = point.p[0]+2*(c.getWorld(1)-c.getWorld(0));	// end x
    }
    else if (index==2)
    {
      rays[j][1] = point.p[1]-2*(c.getWorld(3)-c.getWorld(2));	// end y
    }
    else if (index==3)
    {
      rays[j][1] = point.p[1]+2*(c.getWorld(3)-c.getWorld(2));	// end y
    }
    else if (index==4)
    {
      rays[j][2] = point.p[2]-2*(c.getWorld(5)-c.getWorld(4));	// end z
    }
    else if (index==5)
    {
      rays[j][2] = point.p[2]+2*(c.getWorld(5)-c.getWorld(4));	// end z
    }
    j++;
  }
}	

/** Identify and collect objects crossed during ray tracing
 * from input vertex to extra-object location.
 * \param[in] v Vertex of interest.
 * \param[in] extra_obj_pt Location guaranteed to be
 * outside vertex parent object.
 * \param[out] crossed_objects Collection of object crossed
 * during ray trace from vertex to extra-object point.
 * \return True if any edges were intersected by ray; false otherwise.
 */

bool Nice::getCrossedObjFromVertToExtra (Vertex const * const v,
                                         vector3 const & extra_obj_pt,
                                         vec_op & crossed_objects) const
{
  vector3 origin;
  origin.p[0] = *v->getCoord(0);
  origin.p[1] = *v->getCoord(1);
  origin.p[2] = *v->getCoord(2);
  // find and return crossed objects between origin point and end point
  vec_fp crossed_faces,edge_faces;
  face_grp fg = findIntFacesAlongRay(v,origin,extra_obj_pt,false);
  // prohibit edge intersection, since it is tricky and avoidable
  if (fg.edge_faces.empty()==false) return true;
  // remove current vertex adjacent faces from crossed_faces
  if (fg.crossed_faces.empty()==false)
  {
    sort(fg.crossed_faces.begin(),fg.crossed_faces.end());
    // for each adjacent face
    for (fp_cit i=v->begin();i!=v->end();++i)
    {
      std::pair<fp_it,fp_it> j;	
      j=equal_range(fg.crossed_faces.begin(),fg.crossed_faces.end(),*i);
      // if adjacent face is found in crossed_faces
      // then remove from crossed_faces
      if (j.first!=j.second){fg.crossed_faces.erase(j.first);}
    }
  }
  // identify crossed objects
  getPenetratedObjs(fg.crossed_faces,fg.edge_faces,crossed_objects);
  return false;
}

/** Find and return objects crossed during ray trace 
 * from extra-object point to world limit.
 * \param[in] extra_obj_pt Location guaranteed to be
 * outside of object of interest.
 * \param[in] world_limit Ray will be traced to the world limit.
 * \param[out] crossed_objects Collection of objects crossed during ray trace.
 * \return True if any edges were intersected by ray; false otherwise.
 */

bool Nice::getCrossedObjFromExtraToLimit (vector3 const & extra_obj_pt,
                                          vector3 const & world_limit,
                                          vec_op & crossed_objects) const
{
  vec_fp crossed_faces,edge_faces;
  // keep intersected faces
  face_grp fg = findIntFacesAlongRay(NULL,extra_obj_pt,world_limit,false);
  // prohibit edge intersection, since it is tricky and avoidable
  if (fg.edge_faces.empty()==false) return true;
  // get objects penetrated by ray
  // with endpoints in opposite opposite world partitions,
  // i.e. one inside and one outside.
  getPenetratedObjs(fg.crossed_faces,fg.edge_faces,crossed_objects);
  return false;
}

/** Calculate ray from adjacent face centroid
 * along face normal small distance.
 * \param[out] centroid One endpoint of ray.
 * \param[out] extra_obj_pt Other endpoint of ray.
 * \param[in] f Iterator pointing to adjacent face of vertex of interest.
 */

void Nice::getVertAdjFaceRay (vector3 & centroid,
                              vector3 & extra_obj_pt,
                              fp_cit const & f) const
                              //int const & index) const
{
  // compute centroid of first adjacent face
  //Face * f = v->getAdjacentFace(index);
  for (int i=0;i<3;++i)
  {
    centroid.p[i] = *(*f)->getVertex(0)->getCoord(i)+
                    *(*f)->getVertex(1)->getCoord(i)+
                    *(*f)->getVertex(2)->getCoord(i);
  }
  centroid *= 1.0/3.0;
  // get normal info
  vector3 const * n = (*f)->getNormal();
  double L = Controls::instance().get_epsilon()/sqrt(n->dot(*n));
  extra_obj_pt = centroid+(*n)*L;
}

/** Check a collection of faces for interection with input line segment.
 * \param[in] v Vertex of interest.
 * \param[in] origin One end of line segment.
 * \param[in] end Other end of line segment.
 * \return Collection of faces that interect input line segment
 * along an edge or strictly interior of face.
 */

face_grp Nice::findIntFacesAlongRay (Vertex const * const v,
                                     vector3 const & origin,
                                     vector3 const & end,
                                     bool gate) const
{
  face_grp fg;
  // DEBUG
  vector3 lower(0,0,0);
  vector3 upper(0,0,0);
  for (int i=0;i<3;i++)
  {
    if (origin.p[i]<end.p[i])
    {
      lower.p[i] = origin.p[i];
      upper.p[i] = end.p[i];
    }
    else
    {
      lower.p[i] = end.p[i];
      upper.p[i] = origin.p[i];
    }
  }
  Octree_Visitor_Face visitor(Vector3r(lower.p[0],lower.p[1],lower.p[2]),
                              Vector3r(   upper.p[0],   upper.p[1],   upper.p[2]));
  // DEBUG
  // make a visitor
  //Octree_Visitor_Face visitor(Vector3r(origin.p[0],origin.p[1],origin.p[2]),
  //                            Vector3r(   end.p[0],   end.p[1],   end.p[2]));
  // execute visitor
  Container::instance().octree->visit( visitor );
  // DEBUG
  if (gate)
  {
    cout << "\nNice::findIntFacesAlongRay : ";
    v->print(cout);
    cout << "Nice::findIntFacesAlongRay : "
        << "ray origin = [" << origin.p[0]
        << " " << origin.p[1]
        << " " << origin.p[2] << "]" << endl;
    cout << "Nice::findIntFacesAlongRay : "
        << "ray end = [" << end.p[0]
        << " " << end.p[1]
        << " " << end.p[2] << "]" << endl;
    cout << "Nice::findIntFacesAlongRay : "
          << "number of faces to check = " << visitor.myend()-visitor.mybegin() << endl;
  }
  // DEBUG
  // for each face assumed to contain only unique elements
  for (fp_it j=visitor.mybegin();j!=visitor.myend();++j)
  {
    // DEBUG
    if (gate)
    {
      cout << "Nice::findIntFacesAlongRay :\n";
      (*j)->print(cout);
    }
    // DEBUG
    (*j)->clearFlag();
    // skip adjacent faces since vertex is coincident
    // with the face which would result in edge intersection
    if (v!=NULL)
    {
      if (v->faceIsAdjacent(*j)==true) continue;
    }
    // check for intersection
    result r;
    intersect_triangle3(&origin,&end,
                        (*j)->getNormal(),
                        (*j)->getVertex(0)->getPos(),
                        (*j)->getVertex(1)->getPos(),
                        (*j)->getVertex(2)->getPos(),r);
    // does point intersect polygon
    if (r.line_flag==true)
    {
      if (r.poly_edge_flag==true)
      {
        // add polygon_index to edge array
        fg.edge_faces.push_back(*j);
      }
      else if (r.poly_flag==true)
      {
        // add polygon_index to crossed array
        fg.crossed_faces.push_back(*j);
      }
    }
  }
  return fg;
}

/** Collect the objects for which ray not only intersected
 * but certainly changed sides of object,
 * i.e. inside to outside or outside to inside.
 * \param[in] crossed_faces Faces intersected by ray
 * strictly in face interior.
 * \param[in] crossed_faces_on_edge Faces intercted by ray on edge.
 * \param[out] odd_objects Objects for which ray changed sides.
 */

void Nice::getPenetratedObjs (vec_fp & crossed_faces,
                              vec_fp & crossed_faces_on_edge,
                              vec_op & odd_objects) const
{
  vec_op ol;
  // if any faces were intersected on edge by ray
  if (crossed_faces_on_edge.empty()==false)
  {
    vec_op objs,tmp_objs;
    // for each face crossed on edge
    for (fp_it i=crossed_faces_on_edge.begin();
              i!=crossed_faces_on_edge.end();++i)
    {
      // record face parent object
      objs.push_back(const_cast<Object*>((*i)->getVertex(0)->getObject()));
    }
    // coalesce pairs of edge hits into single object hits
    getObjectsFromEdgeHits(objs,tmp_objs);
    // keep odd objects 
    findOddObjects(tmp_objs,ol);
  }

  // for each face crossed interior
  for (fp_it i=crossed_faces.begin();i!=crossed_faces.end();++i)
  {
    // add face parent object to odd objects from last step
    ol.push_back(const_cast<Object*>((*i)->getVertex(0)->getObject()));
  }
  findOddObjects(ol,odd_objects);
}

/* Filter input object collection by keeping one object if odd number
 * of object name found and zero if an even count is made.
 * \param[in] input_objects Collection of object crossings by ray.
 * \param[out] output_objects Keep only objects crossed an odd number of times.
 */

void Nice::findOddObjects (vec_op & input_objects,
                           vec_op & output_objects) const
{
  // sort object list
  sort(input_objects.begin(),input_objects.end());

  output_objects.clear();
  int count=input_objects.size();
  int k=0;
  // while more Object*s to process
  while (k<count)
  {
    // if not on last Object*
    if (k+1!=count)
    {
      // skip identical pairs of object indices
      if (input_objects[k]==input_objects[k+1]) {k++;k++;}
      // odd object
      else {
        output_objects.push_back( input_objects[k]);
        k++;
      }
    }
    else
    { // add remaining object to odd object list
      output_objects.push_back(input_objects[k]);
      k++;
    }
  }
}

/** Coalesce pairs of same-object edge hits into single object hits.
 * \param[in] edge_hits Collection of objects intersected on edge by ray.
 * \param[out] objs Collection of intersected objects inferred from edge hits.
 */

void Nice::getObjectsFromEdgeHits (vec_op & edge_hits,vec_op & objs) const
{
  objs.clear();
  int count = edge_hits.size();
  // Expecting each edge hit to register twice in record,
  // once for each facec adjacent to edge.
  // Assert an even number of edge hits,
  assert((count%2)==0);
  // sort edge hit list by parent object
  sort(edge_hits.begin(),edge_hits.end());
  // for each Object
  int k=0;
  for (op_it i=edge_hits.begin();i!=edge_hits.end();++i)
  {
    // skip every other object in container
    if ((++k%2)==0) continue;
    // for all but last object
    if ((i+1)!=edge_hits.end())
    {
      // if this element is identical to next
      if (*i==*(i+1))
      {
        // record object
        objs.push_back(*i);
      }
      else
      {
        // error
        cout << "\n\nNice::getObjectsFromEdgeHits:  "
             << "Error. Asymmetric edge intersection. "
             << "Only one adjacent face parent object reported.\n";
        assert(*i==*(i+1));
        exit(1);
      }
    }
  }
}

/** Count the number of nonnice vertices in model.
 * \param[in] detect_self If true then only count vertices
 * inside their parent object. 
 *\return Number of nonnice vertices.
 */

int Nice::getNonniceCount (bool detect_self)
{
  // if detect_self == true
  // then only count self intersections
  // else count all intersections
  int count = 0;
  // for each element of nice
  for (vhm_it j=nonnice.begin();j!=nonnice.end();++j)
  {
    // if zero then should have been removed
    assert( ((*j).second>0) && ((*j).second<3) );
    if (detect_self==true)
    {
      if ((*j).second==2)
      {
        count++;
      }
    }
    else
    {
      count++;
    }
  }
  return count;
}


