#include "edge.h"

#include <cmath>
#include <iostream>
#include <stdlib.h>

#include "object.h"

using std::cout;
using std::endl;

//bool Edge::valid (void)
//{
//  // if edge has two adjacent faces
//  if (f1!=NULL && f2!=NULL)
//  {
//    // and the faces are of interest
//    if ((f1->index==91824 && f2->index==91825 ) || 
//        (f1->index==91825 && f2->index==91824 ) )
//    {
//      cout << "Edge::valid: <" << this << "> vv1=" << vv1->index 
//            << ", vv2=" << vv2->index << endl;
//      printEdge(vv1->o->name);
//      // if the edge vertex indices match the mesh
//      if ((vv1->index==9430 && vv2->index==9632)==false &&
//          (vv1->index==9632 && vv2->index==9430)==false 
//         ){return false;}
//    }
//  }
//  return true;
//}


void Edge::printCP (std::ostream & target)
{
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  getVertices(v1,v2,o1,o2);
  target.precision(12);
  target << "Edge::print: " << this << endl;
  target << "Edge::print: v1:\n";
  if (v1==NULL){target << "v1 is NULL\n";}
  else { v1->printCP(target); }
  target << "Edge::print: v2:\n";
  if (v2==NULL){target << "v2 is NULL\n";}
  else { v2->printCP(target); }
  if (v1!=vv1)
  {
    target << "Edge::print: "
          << "Error vertices (v1,vv1) do not match:\n";
    target << "Edge::print: v1:\n";
    if (v1==NULL){target << "v1 is NULL\n";}
    else { v1->print(target); }
    target << "Edge::print: vv1:\n";
    if (vv1==NULL){target << "vv1 is NULL\n";}
    else { vv1->print(target); }
  }
  if (v2!=vv2)
  {
    target << "Edge::print: "
          << "Error vertices (v2,vv2) do not match:\n";
    target << "Edge::print: v2:\n";
    if (v2==NULL){target << "v2 is NULL\n";}
    else { v2->print(target); }
    target << "Edge::print: vv2:\n";
    if (vv2==NULL){target << "vv2 is NULL\n";}
    else { vv2->print(target); }
  }
  if ((v1!=vv1)|| (v2!=vv2))
  {
    assert(v1==vv1 && v2==vv2);
    exit(1);
  }
}

Edge::Edge (Face *f,Vertex *va,Vertex *vb)
  :vv1(va),vv2(vb),f1(f),f2(NULL),fvec(),l(0)
{
  // compute original edge length
  l=sqrt((va->getpN(0)-vb->getpN(0))*(va->getpN(0)-vb->getpN(0))+
         (va->getpN(1)-vb->getpN(1))*(va->getpN(1)-vb->getpN(1))+
         (va->getpN(2)-vb->getpN(2))*(va->getpN(2)-vb->getpN(2)));
}

double Edge::getSqLength (void)
{
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  getVertices(v1,v2,o1,o2);
  return (v1->getpN(0)-v2->getpN(0))*(v1->getpN(0)-v2->getpN(0))
        +(v1->getpN(1)-v2->getpN(1))*(v1->getpN(1)-v2->getpN(1))
        +(v1->getpN(2)-v2->getpN(2))*(v1->getpN(2)-v2->getpN(2));
}

void Edge::print (std::ostream & target)
{
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  getVertices(v1,v2,o1,o2);
  target.precision(12);
  target << "Edge::print: " << this << endl;
  target << "Edge::print: v1: ";
  if (v1==NULL){target << "v1 is NULL\n";}
  else { v1->print(target); target << "\n";}
  target << "Edge::print: v2: ";
  if (v2==NULL){target << "v2 is NULL\n";}
  else { v2->print(target); target << "\n";}
  target << "Edge::print: f1:\n";
  if (f1==NULL){target << "f1 is NULL\n";}
  else { f1->print(target); }
  target << "Edge::print: f2:\n";
  if (f2==NULL){target << "f2 is NULL\n";}
  else { f2->print(target); }
  if (v1!=vv1)
  {
    target << "Edge::print: "
          << "Error vertices (v1,vv1) do not match:\n";
    target << "Edge::print: v1:\n";
    if (v1==NULL){target << "v1 is NULL\n";}
    else { v1->print(target); }
    target << "Edge::print: vv1:\n";
    if (vv1==NULL){target << "vv1 is NULL\n";}
    else { vv1->print(target); }
  }
  if (v2!=vv2)
  {
    target << "Edge::print: "
          << "Error vertices (v2,vv2) do not match:\n";
    target << "Edge::print: v2:\n";
    if (v2==NULL){target << "v2 is NULL\n";}
    else { v2->print(target); }
    target << "Edge::print: vv2:\n";
    if (vv2==NULL){target << "vv2 is NULL\n";}
    else { vv2->print(target); }
  }
  if ((v1!=vv1)|| (v2!=vv2))
  {
    assert(v1==vv1 && v2==vv2);
    exit(1);
  }
}

double Edge::getAngle (void) 
{
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  getVertices(v1,v2,o1,o2);
  // get outward normals of edge faces
  double n1[3],n2[3];
  f1->getNormal(n1);
  f2->getNormal(n2);
  // compute the cosine of angle between normals
  double normal_angle_cosine=dot(n1,n2)/sqrt(dot(n1,n1))/sqrt(dot(n2,n2));
  // compute angle between normals 
  if 		(normal_angle_cosine >= 1)	
  {
    return Controls::instance().get_pi();
  }
  else if (normal_angle_cosine <= -1)
  {
    // normal_angle = PI;
    // gamma == 0 or 2PI
    return 0;
  }
  else 
  {
    // normal_angle = acos(normal_angle_cosine);
    // use the edge itself as a reference vector
    double refvec[3] = {v2->getpN(0)-o2->getpN(0),v2->getpN(1)-o2->getpN(1),v2->getpN(2)-o2->getpN(2)};
    // dot product of refvec and n1
    double d = dot(refvec,n1);
    if (!d)
    {
      return Controls::instance().get_pi();
    }
    else 
    {
      return Controls::instance().get_pi()+d/fabs(d)*acos(normal_angle_cosine);
    }
    // I EXPECT 0 <= angle < 2*PI
  }
}

void Edge::getVertices (Vertex *&v1,Vertex *&v2,Vertex *&o1,Vertex *&o2)
{
  if (f1!=NULL && f2!=NULL)
  {
    v1=vv1;
    v2=vv2;
    // find o1 on f1
    if      (f1->ptr_vertex(0)!=vv1 && f1->ptr_vertex(0)!=vv2) { o1=f1->ptr_vertex(0); }
    else if (f1->ptr_vertex(1)!=vv1 && f1->ptr_vertex(1)!=vv2) { o1=f1->ptr_vertex(0); }
    else if (f1->ptr_vertex(2)!=vv1 && f1->ptr_vertex(2)!=vv2) { o1=f1->ptr_vertex(0); }
    else 
    {
      cout << "\n\nEdge::getVertices: o1 not identified!\n";
      cout << "	v1:\n";
      v1->print(cout);
      cout << endl << "	v2:\n";
      v2->print(cout);
      cout << endl << "	o1:\n";
      o1->print(cout);
      cout << endl << "	o2:\n";
      o2->print(cout);
      cout << endl << "	vv1:\n";
      vv1->print(cout);
      cout << endl << "	vv2:\n";
      vv2->print(cout);
      cout << endl;
      cout << endl;
      exit(1);
    }
    // find o2 on f2
    if      (f2->ptr_vertex(0)!=vv1 && f2->ptr_vertex(0)!=vv2) { o2=f2->ptr_vertex(0); }
    else if (f2->ptr_vertex(1)!=vv1 && f2->ptr_vertex(1)!=vv2) { o2=f2->ptr_vertex(0); }
    else if (f2->ptr_vertex(2)!=vv1 && f2->ptr_vertex(2)!=vv2) { o2=f2->ptr_vertex(0); }
    else 
    {
      cout << "\n\nEdge::getVertices: o2 not identified!\n";
      cout << "	v1:\n";
      v1->print(cout);
      cout << endl << "	v2:\n";
      v2->print(cout);
      cout << endl << "	o1:\n";
      o1->print(cout);
      cout << endl << "	o2:\n";
      o2->print(cout);
      cout << endl << "	vv1:\n";
      vv1->print(cout);
      cout << endl << "	vv2:\n";
      vv2->print(cout);
      cout << endl;
      cout << endl;
      exit(1);
    }
  }
  else 
  {
    // use f1 only
    v1=vv1;
    v2=vv2;
    // find o1 on f1
    if      (f1->ptr_vertex(0)!=vv1 && f1->ptr_vertex(0)!=vv2) { o1=f1->ptr_vertex(0); }
    else if (f1->ptr_vertex(1)!=vv1 && f1->ptr_vertex(1)!=vv2) { o1=f1->ptr_vertex(0); }
    else if (f1->ptr_vertex(2)!=vv1 && f1->ptr_vertex(2)!=vv2) { o1=f1->ptr_vertex(0); }
    else 
    {
      cout << "\n\nEdge::getVertices: o1 not identified!\n";
      cout << "	v1:\n";
      v1->print(cout);
      cout << endl << "	v2:\n";
      v2->print(cout);
      cout << endl << "	o1:\n";
      o1->print(cout);
      cout << endl << "	o2:\n";
      o2->print(cout);
      cout << endl << "	vv1:\n";
      vv1->print(cout);
      cout << endl << "	vv2:\n";
      vv2->print(cout);
      cout << endl;
      cout << endl;
      exit(1);
    }
  }
}

void Edge::update (Face *f)
{
  //add face to edge
  if      (f1==NULL) { f1=f; }
  else if (f2==NULL) { f2=f; }
  else               { fvec.push_back(f); }
  // add edge pointer to face
  f->addEdge(this);
}

bool Edge::isConsistent (void)
{
  if (f2!=NULL)
  { 
    // not a border edge
    bool forward=false;
    Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
    getVertices(v1,v2,o1,o2);
    // if match v1->v2
    if ((f1->ptr_vertex(0)==v1 && f1->ptr_vertex(1)==v2 )||
        (f1->ptr_vertex(2)==v1 && f1->ptr_vertex(0)==v2 )||
        (f1->ptr_vertex(1)==v1 && f1->ptr_vertex(2)==v2 )
       ){forward=true;}
    // if match v1->v2
    if ((f2->ptr_vertex(0)==v1 && f2->ptr_vertex(1)==v2 )||
        (f2->ptr_vertex(2)==v1 && f2->ptr_vertex(0)==v2 )||
        (f2->ptr_vertex(1)==v1 && f2->ptr_vertex(2)==v2 )
       )
    {
      if (forward){ return false; }
    }
    else 
    {
      if (!forward){ return false; }
    }
  } 
  return true;
}

bool Edge::getStartingFace (Face* &sf)
{
  // if edge is manifold
  if (isManifold()==true)
  {
    // grab starting face
    if (f1!=NULL)
    {
      sf = f1;
    }
    else if (f2!=NULL) 
    {
      sf = f2;
    }
    else 
    {
      cout << "Error. Both edge faces are NULL!\n"; exit(1);
    }
    return true;
  }
  else 
  {
    return false;
  }
}

Face* Edge::getNewFace(Face *old)
{
  // return the adjacent face that is
  // different from old adjacent face
  if (old!=f1 && old!=f2)
  {
    cout << "Error. Neither edge adjacent face matches old face.\n";
    exit(1);
  }
  if (old==f1) { return f2; } 
  else         { return f1; }
}

bool Edge::isManifold(void)
{
  // if three or more faces share an edge then NOT manifold
  // fvec stores adjacent faces beyond first and second faces
  // so if fvec is empty, then edge is manifold
  return fvec.empty()==true;
}

