#include "vertex.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <string.h>

#include "edge.h"
#include "object.h"

using std::cout;
using std::endl;

void Vertex::setVertexNiceness (int val)
{
  o->setVertexNiceness(this,val);
}

int Vertex::getVertexNiceness (void)
{
  return o->getVertexNiceness(this);
}

bool Vertex::vertexIsNice (void)
{
  return o->vertexIsNice(this);	
}

Vertex::Vertex (char* triplet,Object *q)
  :index(0),cl(NULL),o(q),e(),f(),nf()
{

  char val[80];
  char *eptr;
  int i;

  char *cp=triplet;

  // get past 'Vertex'
  while (strchr("Vertx",*triplet)!=NULL) {triplet++;}

  // grab vertex index
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  index = atoi(val);
  if (val==eptr)
  {
    index=0;
    printf("Error in reading vertex index\n");
    return;
  }

  // grab x coord
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  pN[0]=strtod(val,&eptr);
  if (val==eptr)
  {
    printf("Error in reading vertex\n");
    printf("Error in reading vertex: string %s\n",cp);
    return;
  }

  // grab y coord
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  pN[1]=strtod(val,&eptr);
  if (val==eptr)
  {
    printf("Error in reading vertex\n");
    printf("Error in reading vertex: string %s\n",cp);
    return;
  }

  // grab z coord
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  pN[2]=strtod(val,&eptr);
  if (val==eptr)
  {
    printf("Error in reading vertex\n");
    printf("Error in reading vertex: string %s\n",cp);
    return;
  }
  pC[0]=pN[0];
  pC[1]=pN[1];
  pC[2]=pN[2];
}

void Vertex::getNormal (double *n) 
{
  f_iterator i;
  double t[3],theta,thetaT=0,L;
  n[0]=n[1]=n[2]=0;
  // for each adjacent face
  for (i=f.begin();i!=f.end();i++) 
  {
    // get coordinates of polygon normal
    (*i)->getNormal(t);
    L=sqrt(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]);
    theta=(*i)->getAngle(this);
    thetaT+=theta;
    // and add to sum
    n[0] += t[0]/L*theta;
    n[1] += t[1]/L*theta;
    n[2] += t[2]/L*theta;
  }
  n[0] = n[0]/f.size()/thetaT;
  n[1] = n[1]/f.size()/thetaT;
  n[2] = n[2]/f.size()/thetaT;
}

void Vertex::getAdjacentVertices (std::vector<Vertex*> &a)
{
  a.clear();
  e_iterator i;
  // for each adjacent edge
  for (i=e.begin();i!=e.end();i++) 
  {
    Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
    (*i)->getVertices(v1,v2,o1,o2);
    // find vertex different from self and add different vertex to vector
    if      (v1->index!=index){a.push_back(v1);}
    else if (v2->index!=index) {a.push_back(v2);}
    else { printf("Error. both vertices of edge are equal to current vertex.\n"); exit(1); }
  }
}

void Vertex:: getAdjacentFaces (hset_f &fset)
{
  e_iterator i;
  // for each adjacent edge
  for (i=e.begin();i!=e.end();i++) 
  {
    // add edge faces to set
    fset.insert((*i)->ptr_f1());
    fset.insert((*i)->ptr_f2());
  }
}

bool Vertex::scanAdjFaces (Edge *se,Face *sf,bool &nonman)
{
  // collect touched edges
  vec_e te;
  te.push_back(se);
  // collect touched faces
  vec_f tf;
  tf.push_back(sf);
  // get new edge
  se=sf->getNewEdge(se,this);
  if (se->isManifold()==false)
  {
    // vertex manifoldness cannot be determined
    // because it's complicated, so just return flag
    nonman=true;
    return false;
  }
  // while new edge has not already been touched
  while (find(te.begin(),te.end(),se)==te.end())
  {
    // keep new edge
    te.push_back(se);
    // get new face
    sf=se->getNewFace(sf);
    if (sf==NULL) {break;}
    // keep new face
    tf.push_back(sf);
    // get new edge
    se=sf->getNewEdge(se,this);
    if (se->isManifold()==false)
    {
      // vertex manifoldness cannot be determined
      // because it's complicated, so just return flag
      nonman=true;
      return false;
    }
  }
  // if number of touched faces != number of vertex adjacent faces
  return tf.size()==f.size();
}

bool Vertex::isManifold (bool flag)
{
  //	cout << "Vertex::isManifold: @start: flag=" << flag << endl;
  // try to "walk" around vertex adjacent faces
  // using adjacent face edges
  // if an edge is nonmanifold then abort mission
  // stop when the walk returns to starting adjacent face
  if (e.empty()==true)
  {
    cout << "\n\nVertex::isManifold: Error."
          << " Vertex was not an 'orphan', but has no edges.\n";
    cout << "\n\nVertex::isManifold: Confused vertex:\n";
    print(cout);
    cout << endl;
    exit(1);
  }
  // grab starting edge
  Edge *se=e.front();
  // grab starting face
  Face *sf=NULL,*cw=NULL,*ccw=NULL;
  // if edge is manifold
  if (se->isManifold()==true)
  {
    // grab starting face
    if (se->ptr_f1()!=NULL){ cw = se->ptr_f1(); }
    if (se->ptr_f2()!=NULL) { 
      if (cw==NULL) { cw = se->ptr_f2(); }
      else          { ccw = se->ptr_f2(); }
    }
    if (cw==NULL)
    {
      cout << "Error. Both edge faces are NULL!\n"; exit(1);
    }
  }
  else 
  {
    // vertex manifoldness cannot be determined
    // because the starting edge is not manifold
    // thus determining manifoldness of vertex
    // is complicated, so just return flag
    return flag;
    // TODO: report number of vertices for which manifoldness
    // was not determined and alert user
  }

  // try clockwise
  sf=cw;
  bool nonman=false;
  // if fail i.e. not all adjacent faces touched
  if (scanAdjFaces(se,sf,nonman)==false)
  {
    // if nonmanifold edge found, then bail
    if (nonman==true){return flag;}
    // if ccw face is not NULL
    if (ccw!=NULL)
    {
      // try counter-clockwise
      sf=ccw;
      // if fail i.e. all adjacent faces still not touched
      if (scanAdjFaces(se,sf,nonman)==false)
      {
        // record offending edge
        //o->nonman_v.push_back(this);
        o->addNonmanVertices(this);
        return false;
      }
      else 
      {
        return flag;
      }
    }
    else 
    {
      // record offending edge
      //o->nonman_v.push_back(this);
      o->addNonmanVertices(this);
      return false;
    }
  }
  else 
  {
    // all adjacent faces touched
    return flag;
  }
  /*	if (se->getStartingFace(sf)==false)
        {
  // vertex manifoldness cannot be determined
  // because the starting edge is not manifold
  // thus determining manifoldness of vertex
  // is complicated, so just return flag
  return flag;
  // TODO: report number of vertices for which manifoldness
  // was not determined and alert user
  }
  */
}

