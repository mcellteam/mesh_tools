#include "face.h"

#include <cfloat>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>

using std::cout;
using std::endl;

#include "edge.h"
#include "object.h"
#include "vertex.h"

bool Face::match (int i)
{
  return i==index;
}

Face::Face (Vertex *v1,Vertex *v2,Vertex *v3)
  :index(0),b()
{
  v[0]=v1;
  v[1]=v2;
  v[2]=v3;
}

Face::Face (char *triplet,Object *obj)
  :index(0),b()
{

  std::pair<iv_iterator,iv_iterator> pp;
  std::pair<ib_iterator,ib_iterator> qq;

  e[0]=e[1]=e[2]=NULL;

  char val[80];
  char *eptr;
  int i;

  // get past 'Face'
  while (strchr("Face",*triplet)!=NULL) {triplet++;}

  // grab Face index
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  index = (int) strtod(val,&eptr);
  if (val==eptr)
  {
    index=0;
    v[0]=v[1]=v[2]=0;
    printf("Error in reading face index\n");
    return;
  }

  // grab first vertex index
  while (strchr(" \t,",*triplet)!=NULL) { triplet++; }
  i=0;
  while (strchr("0123456789+-eE.",*triplet)!=NULL)
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  int zz = (int)strtod(val,&eptr);
  // search for index in multimap
  //qq=obj->found.equal_range(zz);
  qq = obj->findIndexInMap2(zz);
  if (qq.first!=qq.second)
  {
    // set all matching elements to true
    for(ib_iterator k=qq.first;k!=qq.second;k++)
    {
      // set flag to true;
      (*k).second=true;
    }
  }
  //pp=obj->vp.equal_range(zz);
  pp = obj->findIndexInMap(zz);
  if (pp.first!=pp.second)
  {
    v[0] = (*(pp.first)).second;
  } else
  {
    v[0] = NULL;
    vi[0]=zz;
    // add index to cs.missing_v
    obj->addMissingVertex(zz);
    // add face to cs.missing_f
    obj->addMissingFace(this);
  }
  if (val==eptr)
  {
    v[0]=v[1]=v[2]=NULL;
    printf("Error in reading vertex index\n");
    return;
  }

  // grab second vertex index
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  zz = (int)strtod(val,&eptr);
  // search for index in multimap
  //qq=obj->found.equal_range(zz);
  qq = obj->findIndexInMap2(zz);
  if (qq.first!=qq.second)
  {
    // set all matching elements to true
    for(ib_iterator k=qq.first;k!=qq.second;k++)
    {
      // set flag to true;
      (*k).second=true;
    }
  }
  //pp=obj->vp.equal_range(zz);
  pp = obj->findIndexInMap(zz);
  if (pp.first!=pp.second)
  {
    v[1] = (*(pp.first)).second;
  } else
  {
    v[1] = NULL;
    vi[1]=zz;
    // add index to cs.missing_v
    obj->addMissingVertex(zz);
    // add face to cs.missing_f
    obj->addMissingFace(this);
  }
  if (val==eptr)
  {
    v[0]=v[1]=v[2]=NULL;
    printf("Error in reading vertex index\n");
    return;
  }

  // grab third vertex index
  while (strchr(" \t,",*triplet)!=NULL) triplet++;
  i=0;
  while (strchr("0123456789+-eE.",*triplet))
  {
    val[i++] = *triplet++;
  }
  val[i]=0;
  zz = (int)strtod(val,&eptr);
  // search for index in multimap
  //qq=obj->found.equal_range(zz);
  qq = obj->findIndexInMap2(zz);
  if (qq.first!=qq.second)
  {
    // set all matching elements to true
    for(ib_iterator k=qq.first;k!=qq.second;k++)
    {
      // set flag to true;
      (*k).second=true;
    }
  }
  //pp=obj->vp.equal_range(zz);
  pp = obj->findIndexInMap(zz);
  if (pp.first!=pp.second)
  {
    v[2] = (*(pp.first)).second;
  } else
  {
    v[2] = NULL;
    vi[2]=zz;
    // add index to cs.missing_v
    obj->addMissingVertex(zz);
    // add face to cs.missing_f
    obj->addMissingFace(this);
  }
  if (val==eptr)
  {
    v[0]=v[1]=v[2]=NULL;
    printf("Error in reading vertex index\n");
    return;
  }
}

void Face::addFaceToTable_intf (void)
{
  // create new face vector in table
  vec_f *nfv = new std::vector<Face*>();
  v[0]->getObject()->addIntFaceLHS(this,nfv);
}

bool Face::faceInTable_intf (void)
{
  return v[0]->getObject()->faceInTable_intf(this);
}

ff_iterator Face::findFaceInTable_intf (void)
{
  return v[0]->getObject()->findFaceInTable_intf(this);
}

void Face::addFaceToVector (Face* f)
{
  //(*v[0]->getObject()->intf[this]).push_back(f);
  v[0]->getObject()->addIntFaceRHS(this,f);
}

void Face::clearFaceFromTable_intf (void)
{
  // if this face is in intf table 
  if (faceInTable_intf())
  {
    // delete vector<face*>*
    //delete v[0]->getObject()->intf[this];
    v[0]->getObject()->removeIntFaceLHS(this);
    // remove element from table
    //v[0]->getObject()->intf.erase(this);
  }
}

void Face::print (std::ostream & target)
{
  target.precision(12);
  target << "Face index " << index << endl;
  if (v[0]!=NULL)
  {
    v[0]->print(target);
    target << "\n";
    //cout << "[v0 "
    //      << v[0]->index << " "
    //      << v[0]->pN[0] << " "
    //      << v[0]->pN[1] << " "
    //      << v[0]->pN[2] << "]\n";
  }
  else
  {
    target << "[v0 "
          << vi[0] << "]\n";
  }

  if (v[1]!=NULL)
  {
    v[1]->print(target);
    target << "\n";
    //cout << "[v1 "
    //      << v[1]->index << " "
    //      << v[1]->pN[0] << " "
    //      << v[1]->pN[1] << " "
    //      << v[1]->pN[2] << "]\n";
  }
  else
  {
    target << "[v1 "
          << vi[1] << "]\n";
  }

  if (v[2]!=NULL)
  {
    v[2]->print(target);
    target << "\n";
    //cout << "[v2 "
    //      << v[2]->index << " "
    //      << v[2]->pN[0] << " "
    //      << v[2]->pN[1] << " "
    //      << v[2]->pN[2] << "]\n";
  }
  else
  {
    target << "[v2 "
          << vi[2] << "]\n";
  }

  if (e[0]!=NULL) { target << "[e0 " << e[0] << "]\n";}
  else            { target << "[e0 NULL]\n"; }
  if (e[1]!=NULL) { target << "[e1 " << e[1] << "]\n";}
  else            { target << "[e1 NULL]\n"; }
  if (e[2]!=NULL) { target << "[e2 " << e[2] << "]\n";}
  else            { target << "[e2 NULL]\n"; }
}

void Face::printCP (std::ostream & target)
{
  target.precision(12);
  if (v[0]!=NULL)
  {
    v[0]->printCP(target);
    //cout << v[0]->pN[0] << " "
    //      << v[0]->pN[1] << " "
    //      << v[0]->pN[2] << " 1 0 0 1\n";
  }
  else
  {
    target << "[v0 "
          << v[0] << "]\n";
  }

  if (v[1]!=NULL)
  {
    v[1]->printCP(target);
    //cout << v[1]->pN[0] << " "
    //      << v[1]->pN[1] << " "
    //      << v[1]->pN[2] << " 1 0 0 1\n";
  }
  else
  {
    target << "[v1 "
          << v[1] << "]\n";
  }

  if (v[2]!=NULL)
  {
    v[2]->printCP(target);
    //cout << v[2]->pN[0] << " "
    //      << v[2]->pN[1] << " "
    //      << v[2]->pN[2] << " 1 0 0 1\n";
  }
  else
  {
    target << "[v2 "
          << v[2] << "]\n";
  }
}

void Face::getVertexCoordinates (double *cpvc[3])
{
  cpvc[0]=v[0]->getpN_ptr();
  cpvc[1]=v[1]->getpN_ptr();
  cpvc[2]=v[2]->getpN_ptr();
}

double Face::getAngle (Vertex *vv)
{
  Vertex *vA=vv,*vB=NULL,*vC=NULL;
  double AB[3],AC[3],abL,acL,costheta;
  // identify face vertices
  if      (v[0]!=vv){vB=v[0];}
  else if (v[1]!=vv){vB=v[1];}
  else if (v[2]!=vv){vB=v[2];}
  if      (v[0]!=vv && v[0]!=vB){vC=v[0];}
  else if (v[1]!=vv && v[1]!=vB){vC=v[1];}
  else if (v[2]!=vv && v[2]!=vB){vC=v[2];}
  // AB,AC
  AB[0]=vB->getpN(0)-vA->getpN(0);
  AB[1]=vB->getpN(1)-vA->getpN(1);
  AB[2]=vB->getpN(2)-vA->getpN(2);
  AC[0]=vC->getpN(0)-vA->getpN(0);
  AC[1]=vC->getpN(1)-vA->getpN(1);
  AC[2]=vC->getpN(2)-vA->getpN(2);
  // lengths
  acL=sqrt( dot(AC,AC) );
  abL=sqrt( dot(AB,AB) );
  costheta=( dot(AB,AC) )/abL/acL;
  cout.precision(12);
  if (costheta > 1) costheta=1;
  if (costheta < -1) costheta=-1;
  return acos(costheta);
}

void Face::addEdge (Edge* ptr)
{
  if      (e[0]==NULL){e[0]=ptr;}
  else if (e[1]==NULL){e[1]=ptr;}
  else if (e[2]==NULL){e[2]=ptr;}
  else { cout << "Error. Tried to add fourth edge to face.\n"
    << "Face " << index 
          << " " << static_cast<int>(v[0]->getIndex())
          << " " << static_cast<int>(v[1]->getIndex())
          << " " << static_cast<int>(v[2]->getIndex())
          << endl;
    exit(1); 
  }
}

void Face::getNormal (double n[3]) 
{
  double uX, uY, uZ, vX, vY, vZ;
  // compute vectors 01 and 12
  if (v[1]==NULL) { cout << "v[1]==NULL\n";cout.flush(); }
  if (v[2]==NULL) { cout << "v[2]==NULL\n";cout.flush(); }
  uX = v[1]->getpN(0)-v[0]->getpN(0);
  uY = v[1]->getpN(1)-v[0]->getpN(1);
  uZ = v[1]->getpN(2)-v[0]->getpN(2);
  vX = v[2]->getpN(0)-v[0]->getpN(0);
  vY = v[2]->getpN(1)-v[0]->getpN(1);
  vZ = v[2]->getpN(2)-v[0]->getpN(2);
  // compute cross product (u x v)
  n[0] = uY*vZ-uZ*vY;
  n[1] = uZ*vX-uX*vZ;
  n[2] = uX*vY-uY*vX;
}

void Face::recordBoxes (std::vector<Box*> &ptr)
{
  b.assign(ptr.begin(),ptr.end());
}

Edge* Face::getNewEdge (Edge *old,Vertex *vv)
{
  // return face edge that contains vertex vv
  // and is different from old edge
  if ((e[0]!=old && e[1]!=old && e[2]!=old) ||
     (v[0]!=vv && v[1]!=vv && v[2]!=vv))
  {
    cout << "\n\nOld edge does not match any on face.\n";
    cout << "Current face:\n";
    print(cout);
    cout << endl;
    cout << "Old edge:\n";
    old->print(cout);
    cout << endl;
    cout << "Current vertex:\n";
    vv->print(cout);
    cout << endl << endl;

    exit(1);
  }
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  e[0]->getVertices(v1,v2,o1,o2);
  if (e[0]!=old && (v1==vv || v2==vv)){return e[0];}
  e[1]->getVertices(v1,v2,o1,o2);
  if (e[1]!=old && (v1==vv || v2==vv)){return e[1];}
  e[2]->getVertices(v1,v2,o1,o2);
  if (e[2]!=old && (v1==vv || v2==vv)){return e[2];}
  cout << "\n\nFace::getNewEdge: Error. No edge on face contains current vertex \n"
        << "and is different from old vertex.\n\n";
  cout << "current vertex:\n";
  vv->print(cout);
  cout << endl << endl;
  cout << "old edge:\n";
  old->print(cout);
  cout << endl << endl;
  cout << "current face:\n";
  print(cout);
  cout << endl << endl;
  exit(1);
}

double Face::getAspectRatio (void)
{
  /* Make triangle edge vectors */
  double va[3]={v[1]->getpN(0)-v[0]->getpN(0),v[1]->getpN(1)-v[0]->getpN(1),v[1]->getpN(2)-v[0]->getpN(2)};
  double vb[3]={v[2]->getpN(0)-v[1]->getpN(0),v[2]->getpN(1)-v[1]->getpN(1),v[2]->getpN(2)-v[1]->getpN(2)};
  double vc[3]={v[0]->getpN(0)-v[2]->getpN(0),v[0]->getpN(1)-v[2]->getpN(1),v[0]->getpN(2)-v[2]->getpN(2)};
  double vbase[3]={0,0,0};
  double vopp[3]={0,0,0};

  /* Find length of longest edge */
  double lmax=-DBL_MAX;
  double la=sqrt(dot(va,va));
  double lb=sqrt(dot(vb,vb));
  double lc=sqrt(dot(vc,vc));
  if (la>lmax)
  {
    lmax=la;
    vbase[0]=va[0];
    vbase[1]=va[1];
    vbase[2]=va[2];
    vc[0]=v[2]->getpN(0)-v[0]->getpN(0);
    vc[1]=v[2]->getpN(1)-v[0]->getpN(1);
    vc[2]=v[2]->getpN(2)-v[0]->getpN(2);
    vopp[0]=vc[0];
    vopp[1]=vc[1];
    vopp[2]=vc[2];
  }
  if (lb>lmax)
  {
    lmax=lb;
    vbase[0]=vb[0];
    vbase[1]=vb[1];
    vbase[2]=vb[2];
    va[0]=v[0]->getpN(0)-v[1]->getpN(0);
    va[1]=v[0]->getpN(1)-v[1]->getpN(1);
    va[2]=v[0]->getpN(2)-v[1]->getpN(2);
    vopp[0]=va[0];
    vopp[1]=va[1];
    vopp[2]=va[2];
  }
  if (lc>lmax)
  {
    lmax=lc;
    vbase[0]=vc[0];
    vbase[1]=vc[1];
    vbase[2]=vc[2];
    vb[0]=v[1]->getpN(0)-v[2]->getpN(0);
    vb[1]=v[1]->getpN(1)-v[2]->getpN(1);
    vb[2]=v[1]->getpN(2)-v[2]->getpN(2);
    vopp[0]=vb[0];
    vopp[1]=vb[1];
    vopp[2]=vb[2];
  }

  /* Find shortest altitude */
  double ll = sqrt(dot(vbase,vbase));
  vbase[0]=vbase[0]/ll;
  vbase[1]=vbase[1]/ll;
  vbase[2]=vbase[2]/ll;
  double dot_prod = dot(vbase,vopp);
  double alt[3]={vopp[0]-(dot_prod*vbase[0]),
    vopp[1]-(dot_prod*vbase[1]),
    vopp[2]-(dot_prod*vbase[2])};
  double amin=sqrt(dot(alt,alt));

  return lmax/amin;
}

