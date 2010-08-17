#include "edge.h"

#include <cmath>
#include <iostream>

#include "object.h"

using std::cout;
using std::endl;

bool Edge::valid (void)
{
  // if edge has two adjacent faces
  if (f1!=NULL && f2!=NULL)
  {
    // and the faces are of interest
    if ((f1->index==91824 && f2->index==91825 ) || 
        (f1->index==91825 && f2->index==91824 ) )
    {
      cout << "Edge::valid: <" << this << "> vv1=" << vv1->index 
            << ", vv2=" << vv2->index << endl;
      printEdge(vv1->o->name);
      // if the edge vertex indices match the mesh
      if ((vv1->index==9430 && vv2->index==9632)==false &&
          (vv1->index==9632 && vv2->index==9430)==false 
         ){return false;}
    }
  }
  return true;
}

void Edge::printEdgeCP (void)
{
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  getVertices(v1,v2,o1,o2);
  if ((v1!=vv1) || (v2!=vv2))
  {
    cout << "Edge::printEdge: "
          << "vertices don't match:\n"
          << "	v1:\n";
    v1->printVertex(v1->o->name);
    cout << endl << "	v2:\n";
    v2->printVertex(v2->o->name);
    cout << endl << "	o1:\n";
    o1->printVertex(o1->o->name);
    cout << endl << "	o2:\n";
    o2->printVertex(o2->o->name);
    cout << endl << "	vv1:\n";
    vv1->printVertex(vv1->o->name);
    cout << endl << "	vv2:\n";
    vv2->printVertex(vv2->o->name);
    cout << endl;
    exit(1);
  }
  cout.precision(12);
  cout << "printEdge: " << this << endl;
  cout << "printEdge: <obj>" << v1->o->name << endl;
  cout << v1->pN[0] << " "
        << v1->pN[1] << " "
        << v1->pN[2] << " 1 0 0 1\n";
  cout << v2->pN[0] << " "
        << v2->pN[1] << " "
        << v2->pN[2] << " 1 0 0 1\n";
}

Edge::Edge (Face *f,Vertex *va,Vertex *vb)
  :vv1(va),vv2(vb),f1(f),f2(NULL),fvec(),l(0)
{
  // compute original edge length
  l=sqrt((va->pN[0]-vb->pN[0])*(va->pN[0]-vb->pN[0])+
         (va->pN[1]-vb->pN[1])*(va->pN[1]-vb->pN[1])+
         (va->pN[2]-vb->pN[2])*(va->pN[2]-vb->pN[2]));
}

double Edge::getSqLength (void)
{
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  getVertices(v1,v2,o1,o2);
  return (v1->pN[0]-v2->pN[0])*(v1->pN[0]-v2->pN[0])
        +(v1->pN[1]-v2->pN[1])*(v1->pN[1]-v2->pN[1])
        +(v1->pN[2]-v2->pN[2])*(v1->pN[2]-v2->pN[2]);
}

void Edge::printEdge (std::string s)
{
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  getVertices(v1,v2,o1,o2);
  if ((v1!=vv1) || (v2!=vv2))
  {
    cout << "Edge::printEdge: "
          << "vertices don't match:\n"
          << "	v1:\n";
    v1->printVertex(v1->o->name);
    cout << endl << "	v2:\n";
    v2->printVertex(v2->o->name);
    cout << endl << "	o1:\n";
    o1->printVertex(o1->o->name);
    cout << endl << "	o2:\n";
    o2->printVertex(o2->o->name);
    cout << endl << "	vv1:\n";
    vv1->printVertex(vv1->o->name);
    cout << endl << "	vv2:\n";
    vv2->printVertex(vv2->o->name);
    cout << endl;
    exit(1);
  }
  cout.precision(12);
  cout << "printEdge: " << this << endl;
  cout << "printEdge: <obj>" << s << endl;
  cout << "printEdge:"
        << " v1 "<< v1->index << " ["
        << v1->pN[0] << " "
        << v1->pN[1] << " "
        << v1->pN[2] << "]\n";
  cout << "printEdge:"
        << " v2 "<< v2->index << " ["
        << v2->pN[0] << " "
        << v2->pN[1] << " "
        << v2->pN[2] << "]\n"
        << "printEdge:"
        << " f1 "<< f1->index << " ["
        << f1->v[0]->index << " "
        << f1->v[1]->index << " "
        << f1->v[2]->index << "],";
  if (f2!=NULL)
  {
    cout << " f2 "<< f2->index << " ["
          << f2->v[0]->index << " "
          << f2->v[1]->index << " "
          << f2->v[2]->index << "]\n";
  }
  else
  {
    cout << " f2 is NULL\n";
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
    double refvec[3] = {v2->pN[0]-o2->pN[0],v2->pN[1]-o2->pN[1],v2->pN[2]-o2->pN[2]};
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
    if      (f1->v[0]!=vv1 && f1->v[0]!=vv2) { o1=f1->v[0]; }
    else if (f1->v[1]!=vv1 && f1->v[1]!=vv2) { o1=f1->v[0]; }
    else if (f1->v[2]!=vv1 && f1->v[2]!=vv2) { o1=f1->v[0]; }
    else 
    {
      cout << "\n\nEdge::getVertices: o1 not identified!\n";
      cout << "	v1:\n";
      v1->printVertex(v1->o->name);
      cout << endl << "	v2:\n";
      v2->printVertex(v2->o->name);
      cout << endl << "	o1:\n";
      o1->printVertex(o1->o->name);
      cout << endl << "	o2:\n";
      o2->printVertex(o2->o->name);
      cout << endl << "	vv1:\n";
      vv1->printVertex(vv1->o->name);
      cout << endl << "	vv2:\n";
      vv2->printVertex(vv2->o->name);
      cout << endl;
      cout << endl;
      exit(1);
    }
    // find o2 on f2
    if      (f2->v[0]!=vv1 && f2->v[0]!=vv2) { o2=f2->v[0]; }
    else if (f2->v[1]!=vv1 && f2->v[1]!=vv2) { o2=f2->v[0]; }
    else if (f2->v[2]!=vv1 && f2->v[2]!=vv2) { o2=f2->v[0]; }
    else 
    {
      cout << "\n\nEdge::getVertices: o2 not identified!\n";
      cout << "	v1:\n";
      v1->printVertex(v1->o->name);
      cout << endl << "	v2:\n";
      v2->printVertex(v2->o->name);
      cout << endl << "	o1:\n";
      o1->printVertex(o1->o->name);
      cout << endl << "	o2:\n";
      o2->printVertex(o2->o->name);
      cout << endl << "	vv1:\n";
      vv1->printVertex(vv1->o->name);
      cout << endl << "	vv2:\n";
      vv2->printVertex(vv2->o->name);
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
    if      (f1->v[0]!=vv1 && f1->v[0]!=vv2) { o1=f1->v[0]; }
    else if (f1->v[1]!=vv1 && f1->v[1]!=vv2) { o1=f1->v[0]; }
    else if (f1->v[2]!=vv1 && f1->v[2]!=vv2) { o1=f1->v[0]; }
    else 
    {
      cout << "\n\nEdge::getVertices: o1 not identified!\n";
      cout << "	v1:\n";
      v1->printVertex(v1->o->name);
      cout << endl << "	v2:\n";
      v2->printVertex(v2->o->name);
      cout << endl << "	o1:\n";
      o1->printVertex(o1->o->name);
      cout << endl << "	o2:\n";
      o2->printVertex(o2->o->name);
      cout << endl << "	vv1:\n";
      vv1->printVertex(vv1->o->name);
      cout << endl << "	vv2:\n";
      vv2->printVertex(vv2->o->name);
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
    if ((f1->v[0]==v1 && f1->v[1]==v2 )||
        (f1->v[2]==v1 && f1->v[0]==v2 )||
        (f1->v[1]==v1 && f1->v[2]==v2 )
      ){forward=true;}
    // if match v1->v2
    if ((f2->v[0]==v1 && f2->v[1]==v2 )||
        (f2->v[2]==v1 && f2->v[0]==v2 )||
        (f2->v[1]==v1 && f2->v[2]==v2 )
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

