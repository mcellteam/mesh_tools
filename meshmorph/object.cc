// Author: Justin Kinney
// Date: Sep 2008

#include "object.h"

#include <cassert>
#include <iostream>

#include "edge.h"
#include "log.h"

using std::cout;
using std::endl;

Object::Object (const Object& rhs)
  :name(rhs.name),
  fv(rhs.fv),ff(rhs.ff),fe(rhs.fe),
  v(rhs.v),f(rhs.f),e(rhs.e)
{
}

Object& Object::operator = (const Object& rhs)
{
  cout << "Copy assignment operator prohibited on instances of Object class.\n";
  cout << "Object " << rhs.name << endl;
  exit(1);
}

Object::Object (std::string s)
  :name(s),
  fv(NULL),ff(NULL),fe(NULL),
  v(),f(),e()
{
}

/** Create a unique string from two input integers.
 * \param[in] a First integer.
 * \param[in] b Second integer.
 * \param[in] num_digits Number of digits used in maximum integer value. 
 * \return Unique string including two input integers.
 */

std::string Object::keyPair (int const & a,int const & b,int const & num_digits) const
{
//  char str[128],format[32];
//  sprintf(format,"%%0%dd%%0%dd",num_digits,num_digits);
//  if (a<b){ sprintf(str,format,a,b);}
//  else { sprintf(str,format,b,a); }
//  return str;
  if (a<b)
  {
    return Log::instance().format("%0*d%0*d",num_digits,a,num_digits,b);
  }
  else
  {
    return Log::instance().format("%0*d%0*d",num_digits,b,num_digits,a);
  }
}

/** Look for matching recorded edge instance.
 * \param[in] va One vertex of edge from parent face.
 * \param[in] vb Other vertex of edge from parent face.
 * \param[in] hm Map for finding existing edges with vertex indices as key.
 * \param[in] num_digits Number of digits used in maximum integer value. 
 * \return Edge pointer if match found; NULL otherwise;
 */

Edge* Object::findEdge (Vertex const * const va,
                        Vertex const * const vb,
                        map_s_ep const &hm,
                        int const & num_digits) const
{
  std::string s = keyPair(va->getIndex(),vb->getIndex(),num_digits);
  // if element exists given key, then get Edge pointer
  if (hm.count(s)>0)
  {
    msep_cit i = hm.find(s);	
    return (*i).second;
  }
  else { return NULL;}
}

/** Create and record new instance of Edge class.
 * \param[in] face Parent face.
 * \param[in] va One vertex of edge from parent face.
 * \param[in] vb Other vertex of edge from parent face.
 * \param[in] vc Vertex from parent face not on edge.
 * \param[in] hm Map for finding existing edges with vertex indices as key.
 * \param[in] num_digits Number of digits used in maximum integer value. 
 */

void Object::createEdge (Face * const face,
                         Vertex* const va,
                         Vertex* const vb,
                         Vertex* const vc,
                         map_s_ep &hm,
                         int const & num_digits)
{
  // check for space in vector
  if (e.capacity()>e.size())
  {
    // add edge pointer to object
    e.push_back(Edge(face,va,vb,vc));
    // store edge pointer in hash table
    hm[keyPair(va->getIndex(),vb->getIndex(),num_digits)]=&e.back();;
    // add edge pointer to face
    face->addEdge(&e.back());
    // store first edge pointer
    if (e.size()==1)
    {
      fe = &e[0];
    }
  }
  else
  {
    cout << "\n\nObject::createEdge: Error. "
          << "Must increase capacity of edge vector "
          << "which is currently equal to "
          << e.capacity() << endl << endl;
    assert(e.capacity()>e.size());
    exit(1);
  }
}

/** Record input face edge information as Edge class instance. 
 * \param[in] face Parent face.
 * \param[in] va One vertex of edge from parent face.
 * \param[in] vb Other vertex of edge from parent face.
 * \param[in] vc Vertex from parent face not on edge.
 * \param[in] hm Map for finding existing edges with vertex indices as key.
 * \param[in] num_digits Number of digits used in maximum integer value. 
 */

void Object::processEdge (Face * const face,
                          Vertex * const va,
                          Vertex * const vb,
                          Vertex * const vc,
                          map_s_ep & hm,
                          int const & num_digits)
{
  assert(va!=NULL);
  assert(vb!=NULL);
  Edge *ee=findEdge(va,vb,hm,num_digits);
  if (ee!=NULL)
  {
    // add face pointer to edge 
    ee->update(face,vc);
    // add edge pointer to face
    face->addEdge(ee);
  }
  else {createEdge(face,va,vb,vc,hm,num_digits);}
}

/** Determine the length in digits of largest vertex index.
 * \return Number of digits required for larget vertex index.
 */

int Object::setNumDigits (void) const
{
  int max=0;
  // for each vertex in object
  for (v_cit i=v.begin();i!=v.end();++i)
  {
    if (i->getIndex()>max){max=i->getIndex();}
  }
//  char phrase[64];
//  sprintf(phrase,"%d",max);
//  std::string s = phrase;
//
//  std::string s = Log::instance().format("%d",max);
//  return s.length();
  return Log::instance().format("%d",max).length();
}

/** Record in vertex class adjacent faces to each vertex in this object.
 */

void Object::findVertAdj (void)
{
  // for each face in object
  for (f_it j=f.begin();j!=f.end();++j)
  {
    // for each vertex of face
    for (int i=0;i<3;++i)
    {
      Vertex *vv=j->getVertex(i);
      // add face* to vertex
      vv->addAdjacentFace(&(*j));
    }
  }
}

/** Compute and return bounding box of this object along principal directions.
 * \param[out] bounding_box Bounding box of object
 * as (xmin,xmax,ymin,ymax,zmin,zmax).
 */

void Object::boundObject (double * const & bounding_box) const
{
  double xmin,xmax,ymin,ymax,zmin,zmax,x,y,z;
  //initialize mins and maxes
  xmin = *v[0].getCoord(0);
  xmax = *v[0].getCoord(0);
  ymin = *v[0].getCoord(1);
  ymax = *v[0].getCoord(1);
  zmin = *v[0].getCoord(2);
  zmax = *v[0].getCoord(2);
  // for each vertex in object
  for (v_cit i=v.begin();i!=v.end();++i)
  {
    x = *i->getCoord(0);
    y = *i->getCoord(1);
    z = *i->getCoord(2);
    if      (x>xmax) {xmax = x;}
    else if (x<xmin) {xmin = x;}
    if      (y>ymax) {ymax = y;}
    else if (y<ymin) {ymin = y;}
    if      (z>zmax) {zmax = z;}
    else if (z<zmin) {zmin = z;}
  }
  bounding_box[0]=xmin;bounding_box[1]=xmax;
  bounding_box[2]=ymin;bounding_box[3]=ymax;
  bounding_box[4]=zmin;bounding_box[5]=zmax;
}

/** Get recorded name of this object.
 * \return Object name.
 */

std::string const & Object::getName (void) const
{
  return name;
}

/** Record pointer to first face in face vector.
 * \param[in] fptr Pointer to first face.
 */

void Object::setFF (Face * fptr)
{
  ff = fptr;
}

/** Record pointer to first vertex in vertex vector.
 * \param[in] vv Pointer to first vertex.
 */

void Object::setFV (Vertex * vv)
{
  fv = vv;
}

/** Return pointer to recorded first edge in edge vector.
 * \return Pointer to first edge.
 */

Edge * Object::getFE (void) const
{
  return fe;
}

/** Return pointer to recorded first face in face vector.
 * \return Pointer to first face.
 */

Face * Object::getFF (void) const
{
  return ff;
}

/** Return pointer to recorded first vertex in vertex vector.
 * \return Pointer to first vertex.
 */

Vertex * Object::getFV (void) const
{
  return fv;
}

