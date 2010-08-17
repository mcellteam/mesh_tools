#include "object.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <string.h>

using std::cerr;
using std::cout;
using std::endl;

#include "controls.h"
#include "edge.h"
#include "face.h"
#include "vertex.h"

Object::Object (std::string filename)
  :v(),f(),e(),num_digits(0),hm()
{
  scanFile(filename);
}

Object::~Object (void)
{
  vl_iterator i;
  f_iterator j;
  e_iterator k;
  for (i=v.begin();i!=v.end();i++) { delete *i; }
  for (j=f.begin();j!=f.end();j++) { delete *j; }
  for (k=e.begin();k!=e.end();k++) { delete *k; }
  v.clear();
  f.clear();
  e.clear();
}

std::string Object::keyPair (int a,int b) const
{
  char str[128],format[32];
  sprintf(format,"%%0%dd%%0%dd",num_digits,num_digits);
  if (a<b) { sprintf(str,format,a,b); }
  else     { sprintf(str,format,b,a); }
  return str;
}

Edge* Object::findEdge (Vertex* va,Vertex* vb)
{
  Edge *ee=NULL;
  assert(va!=NULL);
  assert(vb!=NULL);
  std::string s = keyPair(va->getIndex(),vb->getIndex());
  // if element exists given key, then get Edge pointer
  if (hm.count(s)>0){ ee=hm[s]; }
  return ee;
}

void Object::createEdge (Face *ff,Vertex* va,Vertex* vb)
{
  // new edge
  Edge *en = new Edge (ff,va,vb);
  // store edge pointer in hash table
  hm[keyPair(va->getIndex(),vb->getIndex())]=en;
  // add edge pointer to face
  ff->addEdge(en);
  // add edge pointer to object
  e.push_back(en);
}

void Object::buildEdge (Face *ff,Vertex *va,Vertex *vb) 
{
  Edge * ee = findEdge(va,vb);
  // NOTE: 
  if (ee!=NULL)
  {
    ee->addFace(ff);
  }
  else 
  {
    createEdge(ff,va,vb);
  }
}

void Object::setNumDigits (void)
{
  int max=0;
  // for each vertex in object
  for (vl_iterator i=v.begin();i!=v.end();i++)
  {
    if ((*i)->getIndex()>max) { max=(*i)->getIndex(); }
  }
  char str[64];
  sprintf(str,"%d",max);
  std::string s = str;
  num_digits = s.length();
}

void Object::createEdges (void) 
{
  hm.clear();
  // determine number of digits in largest vertex index
  setNumDigits();
  // for each face
  for (f_iterator i=f.begin();i!=f.end();i++) 
  {
    buildEdge(*i,(*i)->ptr_vertex(0),(*i)->ptr_vertex(1));
    buildEdge(*i,(*i)->ptr_vertex(1),(*i)->ptr_vertex(2));
    buildEdge(*i,(*i)->ptr_vertex(2),(*i)->ptr_vertex(0));
  }
}

void Object::scanFile (const std::string & filename) 
{
  char line[2048];
  // open file
  FILE *F = fopen(filename.c_str(),"r");
  if (!F) 
  {
    cout <<"Couldn't open input file "
          << filename << endl;
    exit(1);
  }
  else 
  {
    cerr << "\n\n" << "/* ********************** "
          << "OBJECT ********************** */\n";
    cerr << "name: " << filename << endl;
    cerr.flush();
  }
  // create map: vertex_index -> vertex_pointer
  mmap_iv vp;
  // for every line in file
  for (char *str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F)) 
  {
    // skip leading whitespace
    while (strchr(" \t,",*str)!=NULL) { str++;}
    // if first character is V for Vertex, add new linked list class instance
    if (strchr("V",*str)!=NULL)
    {
      Vertex *vv=new Vertex(str);
      addVertex(vv);
      vp.insert(std::make_pair(vv->getIndex(),vv));
    }
    // if first character is F for Face, add new linked list class instance
    else if (strchr("F",*str)!=NULL)
    {
      Face *ff=new Face(str,vp);
      addFace(ff);
    }
  }
  fclose(F);
}

void Object::gatherFreeVertices (vec_v & free_verts) const
{
  free_verts.clear();
  // for each edge
  for (ce_iterator i=e.begin();i!=e.end();i++) 
  {
    // if edge is free
    if ((*i)->isBorder())
    {
      free_verts.push_back((*i)->ptr_vv1());
      free_verts.push_back((*i)->ptr_vv2());
    }
  }
  // sort and keep unique Vertex*
  sort(free_verts.begin(),free_verts.end(),my_ltv());
  free_verts.assign(free_verts.begin(),unique(free_verts.begin(),free_verts.end()));
}

c_dv_iterator Object::getIterator (const map_dv & projections,
                                   Vertex const * current_vertex,
                                   const double & current_projection) const
{
  std::pair<c_dv_iterator,c_dv_iterator> i;
  i = projections.equal_range(current_projection);
  assert(i.first!=i.second);
  for(c_dv_iterator j=i.first;j!=i.second;++j)                                                  
  {
    if ((*j).second==current_vertex) return j;
  }
  cout << "Object::getIterator: Error: "
        << "Matching vertex not found.\n";
  assert(0);
}

void Object::computeDistances (Closest & distances,
                               const vec_v & free_vertices)
{
  // populate hash table
  map_dv projections;
  for (cv_iterator i=free_vertices.begin();i!=free_vertices.end();++i)
  {
    projections.insert(std::make_pair((*i)->getProjection(),*i));
  }
  distances.clear();
  // for each free vertex
  for (cv_iterator i=free_vertices.begin();i!=free_vertices.end();++i)
  {
    getClosestFreeVerts(*i,projections,distances);
  }
  if (free_vertices.size()>0) assert(distances.good_vertex!=NULL);
}

void Object::addFinalFace (const vec_v & free_vertices)
{
  cv_iterator a = free_vertices.begin();
  cv_iterator b = a; b++;
  cv_iterator c = b; c++;
  int max = 0;
  for (f_iterator i=f.begin();i!=f.end();i++) 
  {
    if ((*i)->getIndex()>max) max = (*i)->getIndex();
  }
  max++;
  // find edge associated with vertices a and b
  Edge * ee = findEdge(*a,*b);
  assert(ee!=NULL);
  // from single face associated with edge
  // determine proper orientation of final face
  Face * ff = ee->ptr_f1();
  assert(ff!=NULL);
  assert(ee->ptr_f2()==NULL);
  if (ff->vertSeqFound(*a,*b))
  {
    ff=new Face(max,*c,*b,*a);
  }
  else
  {
    ff=new Face(max,*a,*b,*c);
  }
  addFace(ff);
}

void Object::removeVertex (Vertex const * bad_vertex)
{
  for (vl_iterator i=v.begin();i!=v.end();i++) 
  {
    if (*i==bad_vertex)
    {
      cerr << "Removing vertex = ";
      (*i)->print(cerr);
      cerr << endl;
      v.erase(i);
      return;
    }
  }
  cerr << "\n\nObject::removeVertex: ERROR."
        << "Bad vertex not found in vertex list.\n";
  cerr << "  bad vertex = ";
  bad_vertex->print(cerr);
  cerr << endl;
  for (vl_iterator i=v.begin();i!=v.end();i++) 
  {
    cerr << "  list vertex = ";
    (*i)->print(cerr);
    cerr << endl;
  }
  assert(0);
}

void Object::fixFaces (Vertex * bad_vertex,
                       Vertex * good_vertex)
{
  for (f_iterator i=f.begin();i!=f.end();i++) 
  {
    (*i)->replaceVertices(bad_vertex,good_vertex);
  }
}

void Object::rebuildEdges (void)
{
  e.clear();
  for (f_iterator i=f.begin();i!=f.end();i++) 
  {
    (*i)->clearEdges();
  }
  createEdges();
}

void Object::projectVerts (const vec_v & free_vertices)
{
  for (cv_iterator i=free_vertices.begin();i!=free_vertices.end();i++) 
  {
    (*i)->project();
  }
}

void Object::getClosestFreeVerts (Vertex * current_vertex,
                                  const map_dv & projections,
                                  Closest & distances)
{
  // get projections map iterator pointing to current_vertex
  const double current_projection = current_vertex->getProjection();
  c_dv_iterator low_iterator = getIterator(projections,current_vertex,current_projection); 
  c_dv_iterator high_iterator = low_iterator;
  c_dv_iterator myprojectionsend = projections.end();
  myprojectionsend--;
  while (1)
  {
    // shift iterators by one, if possible
    if (low_iterator!=projections.begin()) low_iterator--;
    if (high_iterator!=myprojectionsend) high_iterator++;
    Vertex * low_vertex = (*low_iterator).second;
    Vertex * high_vertex = (*high_iterator).second;
    // get distances
    const double low_projection  = (*low_iterator).first;
    const double high_projection = (*high_iterator).first;
    const double low_projected_distance  = current_projection - low_projection;
    const double high_projected_distance = high_projection - current_projection;
    const double low_distance = current_vertex->get3Ddistance(low_vertex);
    const double high_distance = current_vertex->get3Ddistance(high_vertex);
    assert ( low_projection <= current_projection);
    assert (high_projection >= current_projection);
    // add to distances
    if ((low_vertex != current_vertex) && 
        (findEdge(low_vertex,current_vertex)==NULL))
    {
      distances.add(low_distance,current_vertex,low_vertex);
    }
    if ((high_vertex != current_vertex) &&
        (findEdge(high_vertex,current_vertex)==NULL))
    {
      distances.add(high_distance,current_vertex,high_vertex);
    }
    // if distance between free vertices is less than projected distance
    if ((distances.distance < low_projected_distance) &&
        (distances.distance < high_projected_distance))
    {
      return;
    }
    if (high_iterator==myprojectionsend &&
        low_iterator==projections.begin()) break;
  }
}

