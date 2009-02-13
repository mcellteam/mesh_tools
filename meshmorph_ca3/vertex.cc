// Author: Justin Kinney
// Date: Sep 2008

#include "vertex.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "edge.h"
#include "nice.h"
#include "controls.h"
#include "container.h"
#include "intersecting_faces.h"

using std::cout;
using std::endl;

Vertex::Vertex (const Vertex& rhs)
:index(rhs.index),cl(rhs.cl),o(rhs.o),
  f(rhs.f),n(rhs.n),p(rhs.p),r(0)
{
}

Vertex& Vertex::operator = (const Vertex& rhs)
{
  cout << "Assignment operator prohibited on instances of vertex class.\n";
  cout << "Vertex " << rhs.index << endl;
  exit(1);
}

/** Create new instance of Vertex class. 
 *
 * \param[in] line Single line from input mesh file.
 * \param[in] q Pointer to parent Object.
 */

Vertex::Vertex (int const & in,double const & x, double const & y,double const & z,Object * const q)
//Vertex::Vertex (int in,double x, double y,double z,Object * const q)
:index(in),cl(NULL),o(q),f(),n(),p(x,y,z),r(0)
{
}

///** Create new instance of Vertex class. 
// *
// * \param[in] line Single line from input mesh file.
// * \param[in] q Pointer to parent Object.
// */
//
//Vertex::Vertex (char const * line,Object * const q)
//:index(0),cl(NULL),o(q),f(),n(),p()
//{
//  int i;
//  char val[80];
//  char * eptr;
//  char const * cp = line;
//
//  // get past 'Vertex'
//  while (strchr("Vertx",*line)!=NULL) {line++;}
//
//  // grab vertex index
//  while (strchr(" \t,",*line)!=NULL) {line++;}
//  i=0;
//  while (strchr("0123456789+-eE.",*line)!=NULL)
//  {
//    val[i++] = *line++;
//  }
//  val[i]=0;
//  index = (int) strtod(val,&eptr);
//  if (val==eptr)
//  {
//    index=0;
//    printf("Error in reading vertex index\n");
//    return;
//  }
//
//  // grab x coord
//  while (strchr(" \t,",*line)!=NULL) {line++;}
//  i=0;
//  while (strchr("0123456789+-eE.",*line)!=NULL)
//  {
//    val[i++] = *line++;
//  }
//  val[i]=0;
//  p.p[0]=strtod(val,&eptr);
//  if (val==eptr)
//  {
//    printf("Error in reading vertex\n");
//    printf("Error in reading vertex: string %s\n",cp);
//    return;
//  }
//
//  // grab y coord
//  while (strchr(" \t,",*line)!=NULL) line++;
//  i=0;
//  while (strchr("0123456789+-eE.",*line))
//  {
//    val[i++] = *line++;
//  }
//  val[i]=0;
//  p.p[1] = strtod(val,&eptr);
//  if (val==eptr)
//  {
//    printf("Error in reading vertex\n");
//    printf("Error in reading vertex: string %s\n",cp);
//    return;
//  }
//
//  // grab z coord
//  while (strchr(" \t,",*line)!=NULL) line++;
//  i=0;
//  while (strchr("0123456789+-eE.",*line))
//  {
//    val[i++] = *line++;
//  }
//  val[i]=0;
//  p.p[2] = strtod(val,&eptr);
//  if (val==eptr)
//  {
//    printf("Error in reading vertex\n");
//    printf("Error in reading vertex: string %s\n",cp);
//    return;
//  }
//}

/** Print identifying vertex information to output stream.
 *
 * \param[in] target Pre-initialized output stream.
 */

void Vertex::print (std::ostream & target) const
{
  char str[1024];
  sprintf(str,"Vertex <obj> %s <ind> %d",o->getName().c_str(),index);
  sprintf(str,"%s [%.15g %.15g %.15g]\n",str,p.p[0],p.p[1],p.p[2]);
  target << str;
}

/** Print to output stream the vertex position
 *  in DREaMM custom points format.
 *
 * \param[in] target Pre-initialized output stream.
 */

void Vertex::printCP (std::ostream &target) const
{
  char str[1024];
  sprintf(str,"%.15g %.15g %.15g 1 0 0 1\n",p.p[0],p.p[1],p.p[2]);
  target << str;
}

/** Calculate and store vertex normal vector computed as a
 *  weighted sum of adjacent face normals.
 */

void Vertex::setNormal (void)
{
  double theta,thetaT=0.0,L;
  n *= 0.0;
  // for each adjacent face
  for (fp_cit i=f.begin();i!=f.end();++i)
  {
    // get coordinates of polygon normal
    vector3 const * t = (*i)->getNormal();
    L=sqrt(t->dot(*t));
    theta=(*i)->getAngleProxy(this);
    thetaT+=theta;
    // and add to sum
    double a = theta/L;
    n.p[0] += (*t).p[0]*a;
    n.p[1] += (*t).p[1]*a;
    n.p[2] += (*t).p[2]*a;
  }
  double a = 1.0/(f.size()*thetaT);
  n *= a;
}

/** Calculate and return the vertex force vector and energy
 *  due to the vertex extracellular width.
 *
 *  \param[out] force The force vector due to the extracellular width
 *  is added to the input force.
 *  \param[in] compute_force If true, then compute force and energy;
 *  otherwise, compute only the energy. 
 */

double Vertex::getEcwForceEnergy (vector3 & force,bool compute_force) const
{
  Controls & cs(Controls::instance());  
  // if vertex has a closest face
  if (cl!=NULL)
  {
    // get closest point and squared extracellular width
    // between this vertex and closest point.
    vector3 closest_point;
    double sqd;
    Container::instance().findClosestPtInFaceToLocation(*getPos(),cl,closest_point,sqd);
    // compute extracellular width
    double sd = sqrt(sqd);
    // compute ecw error (signed value)
    // se > 0  means extracellular width is too large
    // se == 0 means extracellular width is perfect
    // se < 0  means extracellular width is too small
    double se;
    if (Nice::instance().vertexIsNice(this)==false)
    {
      se=sd+cs.get_target_ecw();
    }
    else
    {
      se=sd-cs.get_target_ecw();
    }
    // sanity check
    if (cs.get_target_ecw()==0)
    {
      cout << "Vertex::getEcwForceEnergy: we got problems."
            << "Divide by TARGET_ECW, but TARGET_ECW==0.\n";
      assert(cs.get_target_ecw()!=0);
      exit(1);
    }
    // optimization -> perform division only once
    se = se/cs.get_target_ecw();
    // compute ecw vector
    vector3 s(closest_point-p);
    // if curretn vertex is on surface of neighbor object
    // then use current vertex outward normal as separation vector
    // recompute extracellular width since it is used as normalization factor
    if (sd==0.0)
    {
      s = getNormal();
      sd = sqrt(s.dot(s));
    }
    if (compute_force==true)
    {
      // spring_force = spring constant * stretch
      // let scaled_spring force = spring_force/ecw
      // Note ECW_WEIGHT has already been normalized to range 0 to 1.
      double scaled_spring_force = cs.get_ecw_gain()*cs.get_ecw_weight()*se/sd;
      //    force cartesian component = spring_force * unit vector
      // where unit vector = cartesian ecw component
      //                   / ecw
      // so force cartesian component = scaled_spring_force
      //                           * cartesian ecw component
      for (int i=0;i<3;++i) force.p[i]+=scaled_spring_force*s.p[i];
    }
    // return energy
    return 0.5*cs.get_ecw_gain()*cs.get_ecw_weight()*se*se;
  }
  else
  {
    // cl==NULL, i.e. no closest point
    // therefore no contribution to force
    // peg energy at maximum possible value, i.e. max ecw error
    double max_se = sqrt(cs.get_search_radius_sq());
    return 0.5*cs.get_ecw_gain()*cs.get_ecw_weight()*max_se*max_se;
  }
}

/** Calculate and return force and/or energy at this vertex
 * due to stretch of adjacent edges.
 * \param[out] force Input force plus calculated edge stretch contribution.
 * \param[in] compute_force If true then compute force and energy;
 * otherwise compute energy only.
 * \return Current energy stored in adjacent edge stretch.
 */

double Vertex::getEdgeStretchForceEnergy (vector3 & force,
                                          bool compute_force) const
{
  double en_ergy=0;
  vec_ep e;
  e.reserve(Controls::instance().get_vector_reserve());
  getAdjacentEdges(e);
  double sum = 0.0;
  // compute mean length of adjacent edges
  // for each adjacent edge of current vertex
  for (ep_it j=e.begin();j!=e.end();++j)
  {
    sum += sqrt((*j)->getSqLength());
  }
  double mean = sum/e.size();
  // for each adjacent edge of current vertex
  for (ep_it j=e.begin();j!=e.end();++j)
  {
    // compute force and energy of edge
    en_ergy+=(*j)->getStretchForceEnergy(this,force,mean,compute_force);
  }
  return en_ergy;
}

/** Calculate and return force and/or energy at this vertex
 * due to distribution of edge angles in adjacent faces of this vertex.
 * \param[out] force Input force plus calculated edge stretch contribution.
 * \param[in] compute_force If true then compute force and energy;
 * otherwise compute energy only.
 * \return Current energy stored in adjacent face edge angles.
 */

double Vertex::getEdgeAngleForceEnergy (vector3 & force,
                                        bool compute_force) const
{
  // NOTE : I previously decided NOT to use nonadjacent edges
  // so oe is not being used at this time
  vec_ep ae,oe;
  ae.reserve(Controls::instance().get_vector_reserve());
  oe.reserve(Controls::instance().get_vector_reserve());

  // for each adjacent face of current vertex
  for (fp_cit i=f.begin();i!=f.end();++i)
  {
    // for each face edge
    for (int j=0;j<3;++j)
    {
      // if either edge vertex is current vertex
      if ((*i)->getEdge(j)->getV1()==this ||
          (*i)->getEdge(j)->getV2()==this)
      {
        // add edge to adjacent edge vector, ae
        ae.push_back((*i)->getEdge(j));
      }
      else
      {
        // add edge to non-adjacent edge vector, oe
        oe.push_back((*i)->getEdge(j));
      }
    }
  }
  // keep unique edges
  sort(ae.begin(),ae.end());
  ae.assign(ae.begin(),unique(ae.begin(),ae.end()));
  sort(oe.begin(),oe.end());
  oe.assign(oe.begin(),unique(oe.begin(),oe.end()));
  // for each adjacent edge
  double en_ergy=0;
  for (ep_it i=ae.begin();i!=ae.end();++i)
  {
    // force and energy contribution
    // arbitrarily passed '0' since it does not matter for energy calculation
    en_ergy+=(*i)->getAngleForceEnergy(0,force,compute_force);
  }
  return en_ergy;
}

/** Calculate and return force and/or energy at this vertex
 * due to adjacent face intersections.
 * \param[out] force Input force plus
 * calculated adjacent face intersection contribution.
 */

void Vertex::getAdjFaceIntForce (vector3 & force) const
{
  // for each adjacent face of current vertex
  for (fp_cit i=f.begin();i!=f.end();++i)
  {
    // if adjacent face is intersected
    if (Intersecting_Faces::instance().faceIsIntersectedLHS(*i)==true)
    {
      // add intersection_force
      Intersecting_Faces::instance().getFaceIntersectionForce(*i,force);
    }
  }
}

/** Calculate and return total force at this vertex, but energy is not computed.
 * \param[out] force Input force plus calculated force from all contributions.
 */

void Vertex::getTotalForceEnergy (vector3 & force) const
{
  ///// get force and energy contributions from extracellular width /////
  // true -> compute force
  getEcwForceEnergy(force,true);
  // DEBUG
  if (false) {cout << Container::instance().o.size() << endl;}
  bool print_switch = false;
  double fmag = 0.0;
  vector3 f_old;
  if (print_switch==true)
  {
    double diff = sqrt(force.dot(force))-fmag;
    fmag = fmag+diff;
    f_old = force;
    cout << "Vertex::getTotalForceEnergy: sep force       ["
          << force.p[0] << " "
          << force.p[1] << " "
          << force.p[2] << "]"
          << ", diff = " << diff << endl;
  }

  // DEBUG
  ///// get force and energy contributions from adjacent edge stretching /////
  getEdgeStretchForceEnergy(force,true);
  // DEBUG
  if (print_switch==true)
  {
    double diff = sqrt(force.dot(force))-fmag;
    fmag = fmag+diff;
    cout << "Vertex::getTotalForceEnergy: stretch force   ["
          << force.p[0]-f_old.p[0] << " "
          << force.p[1]-f_old.p[1] << " "
          << force.p[2]-f_old.p[2] << "]"
          << ", diff = " << diff << endl;
    f_old = force;
  }
  // DEBUG
  ////////// get force contribution from edge angles //////////
  getEdgeAngleForceEnergy(force,true);
  // DEBUG
  if (print_switch==true)
  {
    double diff = sqrt(force.dot(force))-fmag;
    fmag = fmag+diff;
    cout << "Vertex::getTotalForceEnergy: angle force     ["
          << force.p[0]-f_old.p[0] << " "
          << force.p[1]-f_old.p[1] << " "
          << force.p[2]-f_old.p[2] << "]"
          << ", diff = " << diff << endl;
    cout << "Vertex::getTotalForceEnergy: resultant force ["
          << force.p[0] << " "
          << force.p[1] << " "
          << force.p[2] << "]\n\n";
  }
  // DEBUG
  ////////// get force contribution from face intersections //////////
  getAdjFaceIntForce(force);
  // DEBUG
  if (print_switch==true)
  {
    double diff = sqrt(force.dot(force))-fmag;
    fmag = fmag+diff;
    cout << "Vertex::getTotalForceEnergy: face intersection force     ["
          << force.p[0]-f_old.p[0] << " "
          << force.p[1]-f_old.p[1] << " "
          << force.p[2]-f_old.p[2] << "]"
          << ", diff = " << diff << endl;
    cout << "Vertex::getTotalForceEnergy: resultant force ["
          << force.p[0] << " "
          << force.p[1] << " "
          << force.p[2] << "]\n\n";
  }
  // DEBUG
}

/** Calculate new vertex position.
 * \param[in] gain Proportionality between force and displacement.
 * \return New vertex location.
 */

vector3 Vertex::getNewPos (double gain)
{
  // get force
  vector3 force;
  getTotalForceEnergy(force);
  // compute new vertex coordinates
  // DEBUG
//  cout << "Vertex::getNewPos: "
//        << "force = ["
//        << force.p[0] << " "
//        << force.p[1] << " "
//        << force.p[2] << "]"
//        << ", gain = " << gain << endl;
  // DEBUG
  force *= gain;
  // DEBUG
//  cout << "Vertex::getNewPos: "
//        << "force = ["
//        << force.p[0] << " "
//        << force.p[1] << " "
//        << force.p[2] << "]\n";
  // DEBUG
  return force+p;
}

/** Collect and return vertices adjacent to this vertex.
 * \param[out] adjacent_vertices Adjacent vertices of this vertex.
 */

void Vertex::getAdjVertices (vec_vp & adjacent_vertices) const
{
  adjacent_vertices.clear();
  vec_ep e;
  getAdjacentEdges(e);
  // for each adjacent edge
  for (ep_it i=e.begin();i!=e.end();++i)
  {
    // find vertex different from self and add different vertex to vector
    if      ((*i)->getV1()!=this) {adjacent_vertices.push_back((*i)->getV1());}
    else if ((*i)->getV2()!=this) {adjacent_vertices.push_back((*i)->getV2());}
    else
    {
      printf("Error. both vertices of edge are equal to current vertex.\n");
      assert((*i)->getV1()!=this || (*i)->getV2()!=this);
      exit(1);
    }
  }
}

/** Calculare and return squared extracellular width of this vertex.
 * \return Squared extracellular width of this vertex.
 */

double Vertex::getSqSepDist (void) const
{
  if (cl==NULL) return Controls::instance().get_search_radius_sq();
  // get closest point
  vector3 closest_point;
  double sqd;
  Container::instance().findClosestPtInFaceToLocation(*getPos(),cl,closest_point,sqd);
  return sqd;
}

/** Set new location for this vertex.
 * \param[in] new_position New candidate position for this vertex.
 */

void Vertex::setNewPos (vector3 const * const new_position)
{
  double d=0;
  ///// if no adjacent vertex has same position as p /////
  vec_vp v;
  getAdjVertices(v);
  // for each adjacent vertex
  for (vp_it i=v.begin();i!=v.end();++i)
  {
    // if adjacent vertex is indistinguishable from p, then displace p
    if ( !distinguishable((*i)->p.p[0],new_position->p[0])&&
         !distinguishable((*i)->p.p[1],new_position->p[1])&&
         !distinguishable((*i)->p.p[2],new_position->p[2]))
    {
      d=Controls::instance().get_my_double_epsilon();
    }
  }
  p = *new_position+d;
}

/** Calculate and return vertex bounding box as bounding box
 * of collection of adjacent faces.
 * \param[out] lower Lower corner of bounding box of collection of vertex adjacent faces.
 * \param[out] upper Upper corner of bounding box of collection of vertex adjacent faces.
 */

void Vertex::getBoundingBox (vector3 & lower,vector3 & upper) const
{
  // assume bb=bounding box = [xmin xmax ymin ymax zmin zmax]
  lower.p[0]=lower.p[1]=lower.p[2]=1E30;
  upper.p[0]=upper.p[1]=upper.p[2]=-1E30;
  // for each adjacent face
  for (fp_cit i=f.begin();i!=f.end();++i)
  {
    // for each vertex of face
    for (int j=0;j<3;++j)
    {
      Vertex *vv=(*i)->getVertex(j);
      if      (vv->p.p[0]<lower.p[0]){lower.p[0]=vv->p.p[0];}
      else if (vv->p.p[0]>upper.p[0]){upper.p[0]=vv->p.p[0];}
      if      (vv->p.p[1]<lower.p[1]){lower.p[1]=vv->p.p[1];}
      else if (vv->p.p[1]>upper.p[1]){upper.p[1]=vv->p.p[1];}
      if      (vv->p.p[2]<lower.p[2]){lower.p[2]=vv->p.p[2];}
      else if (vv->p.p[2]>upper.p[2]){upper.p[2]=vv->p.p[2];}
    }
  }
}

/** Collect and return edges adjacent to this vertex.
 * \param[out] adjacent_edges Adjacent edges of this vertex.
 */

void Vertex::getAdjacentEdges (vec_ep & adjacent_edges) const
{
  adjacent_edges.clear();
  // for each adjacent face
  for (fp_cit i=f.begin();i!=f.end();++i)
  {
    // for each face edge
    for (int j=0;j<3;++j)
    {
      // if either edge vertex is current vertex
      if ((*i)->getEdge(j)->getV1()==this ||
          (*i)->getEdge(j)->getV2()==this)
      {
        // keep edge
        adjacent_edges.push_back((*i)->getEdge(j));
      }
    }
  }
}

/** Calculate and return the squared virtual displacement of this vertex.
 * \param[in] gain Proportionality constant between force and displacement.
 * \return Squared virtual displacement of this vertex.
 */

double Vertex::getSqVirtualDisp (double gain)
{
  // compute new vertex coords
  // q = new holding position coordinates (x,y,z)
  vector3 q = getNewPos(gain);
  vector3 diff = q-p;
  // DEBUG
//  cout << "Vertex::getSqVirtualDisp: "
//        << "sqd = " << diff.dot(diff)
//        << ", gain = " << gain << endl;
  // DEBUG
  return diff.dot(diff);
}

/** Find and return an adjacent face to this vertex
 * that does not contain the input vertex.
 * \param[in] avoided_vertex Vertex that adjacent face must not reference.
 * \return Pointer to adjacent face of this vertex.
 */

Face * Vertex::getFaceNotAdjToVertex (Vertex const * const avoided_vertex) const
{
  // for each adjacent face of tv
  for (fp_cit i=begin();i!=end();++i)
  {
    // if face is not adjacent to vertex cv
    if (avoided_vertex->faceIsAdjacent(*i)==false)
    {
      return *i;
    }
  }
  cout << "Vertex::getFaceNotAdjToVertex: Error:"
        << "no acceptable face was found!\n";
  assert(false);
  exit(1);
}

/** Recalculate normal vectors of faces adjacent to this vertex.
 */

void Vertex::updateAdjacentFaceNormals (void)
{
  // for each adjacent face
  for (fp_cit i=f.begin();i!=f.end();++i)
  {
    (*i)->updateNormal();
  }
}

/** Record bounding box of each adjacent face of this vertex.
 */

void Vertex::recordAdjFaceBoundingBoxes (hxa7241_graphics::Vector3r * const adjacent_face_lower,
                                         hxa7241_graphics::Vector3r * const adjacent_face_upper)
{
  int j=0;
  // for each adjacent face
  for (fp_cit i=f.begin();i!=f.end();++i)
  {
    (*i)->getBoundingBox(adjacent_face_lower[j],adjacent_face_upper[j]);
    j++;
  }
}

/** Calculate region around this vertex that is influenced by moving this vertex.
 * \param[out] lower Lower corner of local region.
 * \param[out] upper Upper corner of local region.
 */

void Vertex::defineLocalRegion (vector3 & lower,vector3 & upper)
{
  Controls & cs(Controls::instance());
  getBoundingBox(lower,upper);
  lower -= cs.get_update_region_size();
  upper += cs.get_update_region_size();
}

