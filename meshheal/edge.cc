#include "edge.h"

#include <iostream>
#include <stdlib.h>

#include "object.h"

using std::cerr;
using std::cout;
using std::endl;

Edge::Edge (Face *f,Vertex *va,Vertex *vb)
  :vv1(va),vv2(vb),f1(f),f2(NULL)
{
}

void Edge::print (std::ostream & target) const
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

void Edge::getVertices (Vertex *&v1,Vertex *&v2,Vertex *&o1,Vertex *&o2) const
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

void Edge::addFace (Face *f)
{
  //add face to edge
  if      (f1==NULL) { f1=f; }
  else if (f2==NULL) { f2=f; }
  else
  {
    cerr << "\nEdge::addFace: ERROR."
          << "Attempted to add third face to edge.\n";
    print(cerr);
    cerr << endl << "Third face:\n";
    f->print(cerr);
    cerr << endl;
    assert(0);
  }
  // add edge pointer to face
  f->addEdge(this);
}
