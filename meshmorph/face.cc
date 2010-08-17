// Author: Justin Kinney
// Date: Sep 2008

#include "face.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "controls.h"
#include "opttritri.h"

using std::cout;
using std::endl;

Face::Face (const Face& rhs)
:index(rhs.index),flag(rhs.flag),e(),v(),n(rhs.n)
{
  for (int i=0;i<3;++i)
  {
    v[i]=rhs.v[i];
    e[i]=rhs.e[i];
  }
  //for (int i=0;i<6;++i) bb[i]=rhs.bb[i];
}

Face& Face::operator = (const Face& rhs)
{
  cout << "Copy assignment operator prohibited on instances of Face class.\n";
  cout << "Face " << rhs.index << endl;
  exit(1);
}

/** Creat instance of Face class.
 * \param[in] line Single line from input mesh file.
 * \param[in] vp Implicit mapping from vertex index to vertex pointer.
 */

Face::Face (int const & in,int const & v1,int const & v2,int const & v3,vec_vp &vp)
:index(in),flag(false),
  e(),
  // '-1' since pointers stored by in 0-based array
  // but mesh format is 1-based vertex indexing
  //v(vp[v1-1],vp[v2-1],vp[v3-1]),
  n()

{
  // '-1' since pointers stored by in 0-based array
  // but mesh format is 1-based vertex indexing
  v[0] = vp[v1-1];
  v[1] = vp[v2-1];
  v[2] = vp[v3-1];
  // update face bounding box
  updateNormal();
}

///** Creat instance of Face class.
// * \param[in] line Single line from input mesh file.
// * \param[in] vp Implicit mapping from vertex index to vertex pointer.
// */
//
//Face::Face (char *line,vec_vp &vp)
//:index(0),flag(false),
//  e(),
//  v(),
//  n()
//
//{
//  for (int i=0;i<3;++i) e[i]=NULL;
//  char val[80];
//  char *eptr;
//  int i;
//
//  // get past 'Face'
//  while (strchr("Face",*line)!=NULL) {line++;}
//
//  // grab Face index
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
//    v[0]=v[1]=v[2]=0;
//    printf("Error in reading face index\n");
//    assert(val!=eptr);
//    exit(1);
//  }
//
//  // grab first vertex index
//  while (strchr(" \t,",*line)!=NULL) {line++;}
//  i=0;
//  while (strchr("0123456789+-eE.",*line)!=NULL)
//  {
//    val[i++] = *line++;
//  }
//  val[i]=0;
//  // '-1' since pointers stored by in 0-based array
//  // but mesh format is 1-based vertex indexing
//  v[0] = vp[(int)strtod(val,&eptr)-1];
//  if (val==eptr)
//  {
//    v[0]=v[1]=v[2]=0;
//    printf("Error in reading vertex index\n");
//    assert(val!=eptr);
//    exit(1);
//  }
//
//  // grab second vertex index
//  while (strchr(" \t,",*line)!=NULL) {line++;}
//  i=0;
//  while (strchr("0123456789+-eE.",*line))
//  {
//    val[i++] = *line++;
//  }
//  val[i]=0;
//  // '-1' since pointers stored by in 0-based array
//  // but mesh format is 1-based vertex indexing
//  v[1] = vp[(int)strtod(val,&eptr)-1];
//  if (val==eptr)
//  {
//    v[0]=v[1]=v[2]=0;
//    printf("Error in reading vertex index\n");
//    assert(val!=eptr);
//    exit(1);
//  }
//
//  // grab third vertex index
//  while (strchr(" \t,",*line)!=NULL) {line++;}
//  i=0;
//  while (strchr("0123456789+-eE.",*line))
//  {
//    val[i++] = *line++;
//  }
//  val[i]=0;
//  // '-1' since pointers stored by in 0-based array
//  // but mesh format is 1-based vertex indexing
//  v[2] = vp[(int)strtod(val,&eptr)-1];
//  if (val==eptr)
//  {
//    v[0]=v[1]=v[2]=0;
//    printf("Error in reading vertex index\n");
//    assert(val!=eptr);
//    exit(1);
//  }
//  // update face bounding box
//  updateNormal();
//}

/** Print identifying face information to output stream.
 *
 * \param[in] target Pre-initialized output stream.
 */

void Face::print (std::ostream & target) const
{
//  char str[1024];
//  std::string s = v[0]->getObject()->getName();
//  sprintf(str,"Face <obj> %s <ind> %d\n",s.c_str(),index);
//  target << str;
  target << "Face <obj> "
         << v[0]->getObject()->getName() << " <ind> "
         << index << "\n";
  for (int i=0;i<3;++i) v[i]->print(target);
}

/** Print to output stream the position of each vertex
 * of this face in DREaMM custom points format.
 *
 * \param[in] target Pre-initialized output stream.
 */

void Face::printCP (std::ostream & target) const
{
  v[0]->printCP(target);
  v[1]->printCP(target);
  v[2]->printCP(target);
}

/** Get pointers to face vertex positions.
 * \param[out] v0 Pointer to position of first vertex.
 * \param[out] v1 Pointer to position of second vertex.
 * \param[out] v2 Pointer to position of third vertex.
 */

void Face::getVertexCoord (vector3 const * & v0,
                           vector3 const * & v1,
                           vector3 const * & v2) const
{
  v0 = v[0]->getPos();
  v1 = v[1]->getPos();
  v2 = v[2]->getPos();
}

/** Calculate and retun a ratio of face edge lengths that behaves
 * similarly to interior face angle.
 * \param[in] vv Vertex of interest.
 * \return Ratio of length of opposite edge
 * to sum of adjacent edges to vertex of interest.
 */

double Face::getAngleProxy (Vertex const * const vv) const
{
  // identify face vertices
  Vertex *vA=const_cast<Vertex*>(vv),*vB=NULL,*vC=NULL;
  if      (v[0]!=vv){vB=v[0];}
  else if (v[1]!=vv){vB=v[1];}
  else if (v[2]!=vv){vB=v[2];}
  if      (v[0]!=vv && v[0]!=vB){vC=v[0];}
  else if (v[1]!=vv && v[1]!=vB){vC=v[1];}
  else if (v[2]!=vv && v[2]!=vB){vC=v[2];}
  vector3 AB(*vB->getPos()-*vA->getPos());
  vector3 AC(*vC->getPos()-*vA->getPos());
  vector3 BC(*vC->getPos()-*vB->getPos());
  // lengths
  double acL= AC.dot(AC);
  double abL= AB.dot(AB);
  double cbL= BC.dot(BC);
  return cbL/(acL+abL);
}

/** Add input edge pointer to this face.
 * \param[in] ptr Edge of interest.
 */

void Face::addEdge (Edge * ptr)
{
  if     (e[0]==NULL){e[0]=ptr;}
  else if (e[1]==NULL){e[1]=ptr;}
  else if (e[2]==NULL){e[2]=ptr;}
  else
  {
    cout << "Error. Tried to add fourth edge to face.\n"
          << "Face " << index 
          << " " << (int)v[0]->getIndex()
          << " " << (int)v[1]->getIndex()
          << " " << (int)v[2]->getIndex()
          << endl;
    assert(e[0]!=NULL || e[1]!=NULL || e[2]!=NULL);
    exit(1); 
  }
}

/** Calculate and record face normal (not guaranteed to be unit vector).
*/

void Face::updateNormal (void)
{
  // compute vectors 01 and 12
  vector3 uu(*(v[1]->getPos())-*(v[0]->getPos()));
  vector3 vv(*(v[2]->getPos())-*(v[0]->getPos()));
  // compute cross product (u x v)
  n = uu.cross(vv);
}

/** Compare input face index and object name
 *  to this face's index and object name.
 *
 * \param[in] i Input face index.
 * \param[in] name Input face parent object name.
 * \return True if indices and names match; false otherwise.
 */

bool Face::isMatch (int i, std::string const & name) const
{
  return i==index && name==v[0]->getObject()->getName();
}


/** Calculate bounding box of this face.
 * \param[out] lower Lower corner of bounding box of this face.
 * \param[out] upper Upper corner of bounding box of this face.
 * \return Limiting extent of face in ith direction.
 */

void Face::getBoundingBox (hxa7241_graphics::Vector3r & lower,
                           hxa7241_graphics::Vector3r & upper) const
{
  Minmax x(v[0]->getCoord(0),v[1]->getCoord(0),v[2]->getCoord(0));
  Minmax y(v[0]->getCoord(1),v[1]->getCoord(1),v[2]->getCoord(1));
  Minmax z(v[0]->getCoord(2),v[1]->getCoord(2),v[2]->getCoord(2));
  lower.set(*x.min,*y.min,*z.min); 
  upper.set(*x.max,*y.max,*z.max); 
}

/** Calculate and return aspect ratio of this face.
 * \return Face aspect ratio as longest edge
 * divided by perpendicular triangle height.
 */

double Face::getAspectRatio (vector3 & alt, Vertex const * & worst) const
{
  worst = NULL;
  // Make triangle edge vectors
  vector3 va = *v[1]->getPos()-*v[0]->getPos();
  vector3 vb = *v[2]->getPos()-*v[1]->getPos();
  vector3 vc = *v[0]->getPos()-*v[2]->getPos();
  vector3 vbase,vopp; 

  // Find length of longest edge
  double lmax=-1e30;
  double la=sqrt(va.dot(va));
  double lb=sqrt(vb.dot(vb));
  double lc=sqrt(vc.dot(vc));
  if (la>lmax)
  {
    lmax=la;
    vbase=va;
    vc *= -1.0;
    vopp = vc;
    worst = v[2];
  }
  if (lb>lmax)
  {
    lmax=lb;
    vbase=vb;
    va *= -1.0;
    vopp = va;
    worst = v[0];
  }
  if (lc>lmax)
  {
    lmax=lc;
    vbase=vc;
    vb *= -1.0;
    vopp = vb;
    worst = v[1];
  }

  assert(worst!=NULL);

  // Find shortest altitude
  vbase *= 1.0/sqrt(vbase.dot(vbase));
  vbase *= vbase.dot(vopp);
  alt = vopp-vbase;
  double amin=1.0/sqrt(alt.dot(alt));
  // effectively return unit vector pointing in direction
  alt *= amin;
  return lmax*amin;
}

/** Compute and return the force vector and energy
 * due to the aspect ratio of this face.
 * \param[in] v Vertex of interest.
 * \param[out] force Cumulative force.
 * \param[in] compute_force If true then compute force vector;
 * otherwise compute energy only.
 */

double Face::getAspectRatioForceEnergy (Vertex const * const vv,
                                    vector3 & force,
                                    bool compute_force) const
{
  Controls & cs(Controls::instance());
  vector3 force_direction;
  Vertex const * worst_vertex = NULL;
  // compute aspect ratio of this face
  double aspect_ratio = getAspectRatio(force_direction,worst_vertex);
  if (worst_vertex == vv && aspect_ratio > cs.get_aspect_ratio_threshold())
  {
    // check force direction
    if (worst_vertex != vv)
    {
      force_direction *= -1.0;
    }
    // compute aspect ratio error, where best (i.e. lowest) aspect ratio
    // is that of an equilateral triangle, i.e. 1/sin(60 degrees)~=1.1547
    double aspect_ratio_error = aspect_ratio-1.1547;
    // for now use edge length gain and weight
    // should change name to avoid confusion
    double force_magn = cs.get_aspect_ratio_gain()*aspect_ratio_error;
    if (compute_force==true)
    {
      force_direction *= force_magn;
      force += force_direction;
      //for (int i=0;i<3;++i) force.p[i]+=force_magn*force_direction.p[i];
    }
  }
  // energy contribution
  //return 0.5*force_magn*aspect_ratio_error;
  //not using energy so return 0
  return 0.0;
  //return cs.EDGE_STRETCH_WEIGHT/cs.STRETCH_EXPONENT*force_magn*se;
}

/** Return pointer to parent object of this face.
 * \return Pointer to parent object.
 */

Object const * Face::getObject (void) const
{
  return v[0]->getObject();
}

