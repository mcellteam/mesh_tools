// Author: Justin Kinney
// Date: Sep 2008

#include "edge.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "controls.h"

using std::cout;
using std::endl;

Edge::Edge (Face * const f,Vertex * const va,Vertex * const vb,Vertex * const vc)
:f1(f),f2(NULL),l(0),v1(va),v2(vb),o1(vc),o2(NULL)
{
  // compute original edge length
  vector3 a(*va->getPos()-*vb->getPos());
  l = sqrt(a.dot(a));
}

/** Write edge description to specified stream.
 * \param[in] target Output stream;
 */

void Edge::print (std::ostream &target) const
{
  //Vertex_Group vg = getVertices();
  cout.precision(12);
  if ( v1==NULL || v2==NULL ||
       o1==NULL || o2==NULL ||
       f1==NULL || f2==NULL )
  {
    cout << "\nEdge::print: NULL POINTERS IN EDGE!\n";
    if (v1==NULL){cout << "v1 is NULL\n";}
    if (v2==NULL){cout << "v2 is NULL\n";}
    if (o1==NULL){cout << "o1 is NULL\n";}
    if (o2==NULL){cout << "o2 is NULL\n";}
    if (f1==NULL){cout << "f1 is NULL\n";}
    if (f2==NULL){cout << "f2 is NULL\n";}
    assert(v1!=NULL && v2!=NULL &&
           o1!=NULL && o2!=NULL &&
           f1!=NULL && f2!=NULL);
    exit(1);
  }
  target << "\nprintEdge: v1:\n";
  v1->print(target);
  target << "printEdge: v2:\n";
  v2->print(target);
  target << "printEdge: o1:\n";
  o1->print(target);
  target << "printEdge: o2:\n";
  o2->print(target);
  target << "printEdge: f1:\n";
  f1->print(target);
  target << "printEdge: f2:\n";
  f2->print(target);
  target << "printEdge: orignal length " << l << endl;
}

/** Compute and return the force vector and energy
 * due to the stretch of the edge.
 * \param[in] v Vertex of interest.
 * \param[out] force Cumulative force.
 * \param[in] compute_force If true then compute force vector;
 * otherwise compute energy only.
 */

double Edge::getStretchForceEnergy (Vertex const * const v,
                                    vector3 & force,
                                    double const & mean,
                                    bool compute_force) const
{
  Controls & cs(Controls::instance());
  // get pointer to edge vertex
  // that is not input vertex
  Vertex * av = (v1==v) ? v2 : v1;
  // compute edge vector
  vector3 s(*av->getPos()-*v->getPos());
  // compute edge length
  double sd = sqrt(s.dot(s));
  // compute normalized edge stretch (signed value)
  double se;
  if (cs.get_use_edge_reference_length()==true)
  {
    // use original edge length as reference
    se = (sd-l)/l;
  }
  else
  {
    // compute reference length as mean of edge lengths
    // of first adjacent face to edge
    //double ref = (f1->getEdge(0)->l+f1->getEdge(1)->l+f1->getEdge(2)->l)*0.33333333;
    //double ref = (f1->getEdge(0)->l+f1->getEdge(1)->l+f1->getEdge(2)->l+f2->getEdge(0)->l+f2->getEdge(1)->l+f2->getEdge(2)->l)*0.166666667;
    //se = (sd-ref)/ref;
    se = (sd-mean)/mean;
  }
  double force_magn = cs.get_edge_length_weight()*cs.get_edge_length_gain()*se;
  if (compute_force==true)
  {
    // force contribution
    // spring_force = spring constant * stretch
    // scaled_spring force = spring_force/edge length
    double scaled_spring_force = force_magn/sd;
    //    force cartesian component = spring_force * unit vector
    // where unit vector = (adjacent vertex position
    //                   - target vertex position) / edge length
    // so force cartesian component = scaled_spring_force 
    //                              * cartesian component difference
    for (int i=0;i<3;++i) force.p[i]+=scaled_spring_force*s.p[i];
  }
  // energy contribution
  return 0.5*force_magn*se;
  //return cs.EDGE_STRETCH_WEIGHT/cs.STRETCH_EXPONENT*force_magn*se;
}

/** Calculate and return edge angle.
 * \return Edge angle in radians.
 */

double Edge::getAngle (void) const
{
  Controls & cs(Controls::instance());
  // get outward normals of edge faces
  vector3 const * n1 = f1->getNormal();
  vector3 const * n2 = f2->getNormal();
  // compute the cosine of angle between normals
  double normal_angle_cosine=n1->dot(*n2)/sqrt(n1->dot(*n1))/sqrt(n2->dot(*n2));
  // compute angle between normals 
  if (normal_angle_cosine >= 1)
  {
    return cs.get_pi();
  }
  else if (normal_angle_cosine <= -1)
  {
    return 0;
  }
  else
  {
    // use the edge itself as a reference vector
    vector3 refvec(*(v2->getPos())-*(o2->getPos()));
    // dot product of refvec and n1
    double d = refvec.dot(*n1);
    if (!d) { return cs.get_pi();}
    else   { return cs.get_pi()+d/fabs(d)*acos(normal_angle_cosine); }
    // I EXPECT 0 <= angle < 2*PI
  }
}

/** Get perpendicular length from outer vertex to edge.
 * \param[in] o Outer vertex of edge.
 * \param[in] vv1 One edge vertex.
 * \param[in] vv2 Other edge vertex.
 * \param[in] edge_length Edge Length.
 * \return Perpendicular distance from outer vertex to edge.
 */

double Edge::getCurvatureLength (Vertex const * const o,
                                 Vertex const * const vv1,
                                 Vertex const * const vv2,
                                 double const & edge_length) const
{
  vector3 o_v1(*o->getPos()-*vv1->getPos());
  vector3 v2_v1(*vv2->getPos()-*vv1->getPos());
  double unum = o_v1.dot(v2_v1);
  v2_v1 *= unum/(edge_length*edge_length);
  vector3 I(*o->getPos()-*vv1->getPos()+v2_v1);
  return sqrt(I.dot(I));
}

/** Calculate and return angle force vector and energy.
 * \param[in] o2IsRequesting If nonzero, calculate force and energy at o2 vertex;
 * otherwise calculate for o1 vertex.
 * \param[out] force Cumulative force.
 * \param[in] compute_force If true then compute force vector;
 * otherwise compute energy only.
 * \return Angle energy.
 */

double Edge::getAngleForceEnergy (int const & o2IsRequesting,
                                  vector3 & force,
                                  bool compute_force) const
{
  Controls & cs(Controls::instance());
  // get polygon outward normal
  // if vertex requesting force is o2
  // else vertex requesting force is o1
  vector3 const * n = (o2IsRequesting!=0) ? f2->getNormal() : f1->getNormal();
  // compute normal length
  double L = sqrt(n->dot(*n));
  // compute angle between edge adjacent face normals
  double angle = getAngle();
  // determine force direction
  // when angle is less than PI, then angle_error is negative
  // so force will be opposite direction of face normal
  // when angle is greater than PI, then angle_error is positive
  // so force will be in same direction as face normal
  double angle_error = angle-cs.get_pi();
  int sign = (angle_error<0) ? 1 : -1;
  // get perpendicular lengths L1 and L2 from o1 and o2 vertices to edge
  double edge_length = sqrt(getSqLength());
  double L1 = getCurvatureLength(o1,v1,v2,edge_length);
  double L2 = getCurvatureLength(o2,v1,v2,edge_length);
  // compute normalized force gain, i.e. force metric
  // DEBUG
  double ae = fabs(angle_error/cs.get_pi());
  // DEBUG
  //double ae = fabs(angle_error/cs.get_pi())/(0.5*(L1+L2));
  double force_metric = cs.get_edge_angle_gain()*cs.get_edge_angle_weight()*ae;
  // force 
  if (compute_force==true)
  {
    // further normalize force by appropriate moment arm
    double force_magn;
    if (o2IsRequesting==true)
    {
      force_magn = force_metric/L2;
      //force_magn = force_metric;
    }
    else
    {
      force_magn = force_metric/L1;
      //force_magn = force_metric;
    }
    // divide by normal vector length to account
    // for decomposition into xyz components
    double force_norm = sign*force_magn/L;
    for (int i=0;i<3;++i) force.p[i]+=force_norm*(*n).p[i];
  }
  // energy
  return 0.5*force_metric*ae;
}

/** Calculate and return angle reaction force vector and energy.
 * \param[out] force Input force plus reaction force.
 * \param[in] compute_force If true then compute reaction force vector;
 * otherwise compute energy only.
 * \return Reaction energy.
 */

void Edge::getAngleReForceEnergy (vector3 & force,
                                    bool compute_force) const
{
  Controls & cs(Controls::instance());
  // get polygon outward normals
  vector3 const * n1 = f1->getNormal();
  vector3 const * n2 = f2->getNormal();
  // compute normal length
  double L1 = sqrt(n1->dot(*n1));
  double L2 = sqrt(n2->dot(*n2));
  // compute cosine of angle between normals
  double angle=getAngle();
  // determine force direction
  // when angle is less than PI, then angle_error is negative
  // so force will be opposite direction of face normal
  // when angle is greater than PI, then angle_error is positive
  // so force will be in same direction as face normal
  double angle_error = angle-cs.get_pi();
  double ae = angle_error/cs.get_pi();
  double force_metric = ae;
  if (compute_force==true)
  {
    // force 
    double b = cs.get_edge_angle_gain()*cs.get_edge_angle_weight()*force_metric;
    double force_magn1 = b/L1;
    double force_magn2 = b/L2;
    double force1[3]={force_magn1*(*n1).p[0],force_magn1*(*n1).p[1],force_magn1*(*n1).p[2]};
    double force2[3]={force_magn2*(*n2).p[0],force_magn2*(*n2).p[1],force_magn2*(*n2).p[2]};
    double reaction[3]={-(force1[0]+force2[0]),
      -(force1[1]+force2[1]),
      -(force1[2]+force2[2])};
    for (int i=0;i<3;++i) force.p[i] += 0.5*reaction[i];
  }
}

/** Add input Face pointer to this Edge.
 * \param[in] f Face pointer.
 * \param[in] vc Pointer to vertex of second adjacent face to this edge
 * that is not part of this edge.
 */

void Edge::update (Face * const f,Vertex * const vc)
{
  //add face to edge
  if (f1==NULL) {f1=f;}
  else if (f2==NULL) {f2=f;}
  else
  {
    cout << "Error. Tried to add third face to edge.\n";
    f->print(std::cout);
    assert(f1==NULL || f2==NULL);
    exit(1); 
  }
  o2 = vc;
}

