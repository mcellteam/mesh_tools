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

#include "boundary.h"
#include "box.h"
#include "container.h"
#include "controls.h"
#include "edge.h"
#include "face.h"
#include "face_pair.h"
#include "meshalyzer.h"
#include "space.h"
#include "vertex.h"

Object::Object (std::string s)
  :name(s),v(),f(),e(),
  //closed(false),consistent(false),
  outward(false),//manifold(false),
  iv(),intf(),nice(),flip(),num_sep(0),num_bou(0),vol(0),genus(0),vp(),found(),
  contig_f(false),contig_v(false),area(),aspect_ratio(),edge_length(),edge_angle(),
  adjacent_face(),border(),nonman_e(),flipped(),indistin_v(),indistin(),dupl_v_index(),
  dupl_f_index(),nonman_v(),orphan(),missing_f(),missing_v(),degen(),bad_aspect(),
  bad_angle(),bad_length(),zero_area()

{
  // cumulative statistics
  bb[0]=bb[1]=bb[2]=bb[3]=bb[4]=bb[5]=0;
  vec_cv[0]=vec_cv[1]=vec_cv[2]=vec_cv[3]=vec_cv[4]=0;
  vec_cf[0]=vec_cf[1]=vec_cf[2]=vec_cf[3]=vec_cf[4]=0;
}

Object::~Object (void)
{
  v_iterator i;
  f_iterator j;
  e_iterator k;
  for (i=v.begin();i!=v.end();i++) { delete *i; }
  for (j=f.begin();j!=f.end();j++) { delete *j; }
  for (k=e.begin();k!=e.end();k++) { delete *k; }
  v.clear();
  f.clear();
  e.clear();
}

void Object::clearFlipTable (void)
{
  flip.clear();
}

void Object::addEdgeFlip (Edge* ee,int i)
{
  // if no element in hashtable has this key
  if (flip.find(ee)==flip.end()){flip[ee]=i;}
  else {cout << "Edge element already in hashtable.\n";exit(1);}
}

int Object::getEdgeFlip (Edge* ee)
{
  // if element with this edge key exists in hashtable
  if (flip.find(ee)!=flip.end()){return flip[ee];}
  else {cout << "Edge element does not exist in hashtable.\n";exit(1);}
}

bool Object::vertexIsNice (Vertex *vv)
{
  return nice.find(vv)==nice.end();	
}

int Object::getVertexNiceness (Vertex *vv)
{
  if (vertexIsNice(vv))
  {
    return 0;
  }
  else
  {
    return nice[vv];
  }
}

void Object::setVertexNiceness (Vertex *vv,int val)
{
  if (val==0)
  {
    nice.erase(vv);
  }
  else
  {
    nice[vv]=val;
  }
}

bool Object::faceInTable_intf (Face *ff)
{
  return intf.find(ff)!=intf.end();
}

ff_iterator Object::findFaceInTable_intf (Face *ff)
{
  return intf.find(ff);
}

Edge* Object::findEdge (Vertex* va,Vertex* vb,map_se &hm,int num_digits)
{
  Edge *ee=NULL;
  std::string s = keyPair(va->getIndex(),vb->getIndex(),num_digits);
  // if element exists given key, then get Edge pointer
  if (hm.count(s)>0){ ee=hm[s]; }
  return ee;
}

void Object::createEdge (Face *ff,Vertex* va,Vertex* vb,map_se &hm,int num_digits)
{
  // new edge
  Edge *en = new Edge (ff,va,vb);
  // store edge pointer in hash table
  hm[keyPair(va->getIndex(),vb->getIndex(),num_digits)]=en;
  // add edge pointer to face
  ff->addEdge(en);
  // add edge pointer to object
  e.push_back(en);
}

void Object::buildEdge (Face *ff,Vertex *va,Vertex *vb,map_se &hm,int num_digits) 
{
  Edge *ee=NULL;
  ee=findEdge(va,vb,hm,num_digits);
  // NOTE: 
  if (ee!=NULL)
  {
    ee->update(ff);
  }
  else 
  {
    createEdge(ff,va,vb,hm,num_digits);
  }
}

int Object::setNumDigits (void)
{
  int max=0;
  // for each vertex in object
  for (v_iterator i=v.begin();i!=v.end();i++)
  {
    if ((*i)->getIndex()>max) { max=(*i)->getIndex(); }
  }
  char str[64];
  sprintf(str,"%d",max);
  std::string s = str;
  return s.length();
}

void Object::createEdges (void) 
{
  cerr << "Create edges for object [" << name << "]...............................";
  cerr.flush();
  // determine number of digits in largest vertex index
  int num_digits = setNumDigits();
  // create map for finding edges
  map_se hm;
  f_iterator i;
  // for each face
  for (i=f.begin();i!=f.end();i++) 
  {
    buildEdge(*i,(*i)->ptr_vertex(0),(*i)->ptr_vertex(1),hm,num_digits);
    buildEdge(*i,(*i)->ptr_vertex(1),(*i)->ptr_vertex(2),hm,num_digits);
    buildEdge(*i,(*i)->ptr_vertex(2),(*i)->ptr_vertex(0),hm,num_digits);
  }
  cerr << "complete.\n";cerr.flush();
}

//void Object::verifyEdges (void)
//{
//  for (e_iterator i=e.begin();i!=e.end();i++)
//  {
//    Edge *ee=*i;
//    if (ee->valid()==false)
//    {
//      cout << "Object::verifyEdges: Error! edge is invalid.\n";
//      ee->printEdge(ee->ptr_vv1()->o->name);
//      cout << endl << endl;
//      exit(1);
//    }
//  }
//}

double Object::getMeanEdgeLength (void)
{
  e_iterator i;
  double L = 0;
  // for each edge
  for (i=e.begin();i!=e.end();i++) 
  {
    L += sqrt((*i)->getSqLength());
  }
  // compute mean edge length
  return (L/(double)e.size());
}

bool Object::processEdge (Edge *ee,hmap_vd &hood,vec_e &bucket,Vertex *vv)
{
  bool empty = true;
  bool f1 =false,f2=false;
  Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
  ee->getVertices(v1,v2,o1,o2);
  Vertex *vp1=v1,*vp2=v2;
  // if v1 found in hashtable, then f1 = true;
  if (hood.find(vp1)!=hood.end()){f1=true;}
  // if v2 found in hashtable, then f1 = true;
  if (hood.find(vp2)!=hood.end()){f2=true;}
  // if neither v1 nor v2 are in struct list
  if (!f1&&!f2)
  {
    // then add edge* to bucket
    bucket.push_back(ee);
    empty = false;
  }
  else if (f1 && f2)
  {
  }
  else 
  {
    if (f1)
    {
      // if v1 found and v2 not found in hashtable
      // then add v2 to hashtable
      hood[vp2]=sqrt( (vp2->getpN(0)-vv->getpN(0))*(vp2->getpN(0)-vv->getpN(0))+
                      (vp2->getpN(1)-vv->getpN(1))*(vp2->getpN(1)-vv->getpN(1))+
                      (vp2->getpN(2)-vv->getpN(2))*(vp2->getpN(2)-vv->getpN(2)));
    }
    if (f2)
    {
      // if v2 found and v1 not found in hashtable
      // then add v1 to hashtable
      hood[vp1]=sqrt( (vp1->getpN(0)-vv->getpN(0))*(vp1->getpN(0)-vv->getpN(0))+
                      (vp1->getpN(1)-vv->getpN(1))*(vp1->getpN(1)-vv->getpN(1))+
                      (vp1->getpN(2)-vv->getpN(2))*(vp1->getpN(2)-vv->getpN(2)));
    }
  }
  return empty;
}

void Object::collectFaces (hmap_vd &hood,set_v &disabled,vec_f &new_faces)
{
  new_faces.clear();
  // for each vertex in hood
  for (vd_iterator i=hood.begin();i!=hood.end();i++)
  {
    // if vertex is thawed and not disabled
    if (!ifFrozen(hood,(*i).first) && disabled.find((*i).first)==disabled.end())
    {
      // for each adjacent face of thawed vertex
      for (f_iterator j=(*i).first->f.begin();j!=(*i).first->f.end();j++)
      {
        // if any face vertex is thawed, then add face to collection
        if (!ifFrozen(hood,(*j)->ptr_vertex(0)) ||
            !ifFrozen(hood,(*j)->ptr_vertex(1)) ||
            !ifFrozen(hood,(*j)->ptr_vertex(2))){new_faces.push_back(*j);}
      }
      // add vertex to disabled list
      disabled.insert((*i).first);
    }
  }
}

bool Object::ifFrozen (hmap_vd &neighborhood,Vertex *vv)
{
  // if vertex is in hashtable
  if (neighborhood.find(vv)!=neighborhood.end())
  {
    if (neighborhood[vv]<Controls::instance().get_neighborhood_radius()) { return false; }
    else {return true;}
  }
  return true;
}

bool Object::thawedAndAble (hmap_vd &hood,set_v &disabled)
{
  // for each vertex in hood
  for (vd_iterator i=hood.begin();i!=hood.end();i++)
  {
    // if vertex is thawed and not disabled
    if (!ifFrozen(hood,(*i).first) && disabled.find((*i).first)==disabled.end())
    {
      // then affirm that at least one vertex is thawed and able
      return true;
    }
  }
  return false;
}

void Object::newFindNeighborhoods (void)
{
  // for each vertex in object
  for (v_iterator i=v.begin();i!=v.end();i++) 
  {
    // initialize vertex*->double hash table
    // represents a neighbor vertex and
    // the shortest cumulative edge length
    // to reach it from current vertex
    hmap_vd hood;
    hood.clear();
    // add current vertex to neighbor list
    // naturally, assign it zero length
    hood[*i]=0.0;
    // initialize set of vertices to constitute disabled list
    set_v disabled;
    disabled.clear();
    // init collection of neighborhood faces
    vec_f c;
    ///// initial round /////
    // init collection of new faces
    vec_f new_faces;
    collectFaces(hood,disabled,new_faces);
    // for each face in new collection
    for (f_iterator k=new_faces.begin();k!=new_faces.end();k++)
    {
      // initialize container for edges
      // with neither vertex in hood
      vec_e bucket;
      bucket.clear();
      bool bucket_empty = true;
      // for each edge in face
      for (int j=0;j<3;j++)
      {
        if (!processEdge((*k)->ptr_edge(j),hood,bucket,*i))
        {
          bucket_empty = false;
        }
      }
      if (!bucket_empty)
      {
        // initialize another container for edges
        // with neither vertex in hood
        vec_e pail;
        pail.clear();
        bool pail_empty = true;
        // for each edge in bucket
        for (e_iterator j=bucket.begin();j!=bucket.end();j++)
        {
          if (!processEdge(*j,hood,pail,*i))
          {
            pail_empty = false;
          }
        }
        if (!pail_empty)
        {
          cout << "Error. Multiple rounds of bucket use required.\n";
          exit(1);
        }
      }
    }

    ///// all subsequent rounds /////

    // while there are thawed vertices in neighbor list, hood,
    // that are also not disabled
    while (thawedAndAble(hood,disabled))
    {
      // init collection of new faces
      vec_f new_faces_too;
      collectFaces(hood,disabled,new_faces_too);
      // for each face in new collection
      for (f_iterator k=new_faces_too.begin();k!=new_faces_too.end();k++)
      {
        // initialize container for edges
        // with neither vertex in hood
        vec_e bucket;
        bucket.clear();
        bool bucket_empty = true;
        // for each edge in face
        for (int j=0;j<3;j++)
        {
          if (!processEdge((*k)->ptr_edge(j),hood,bucket,*i))
          {
            bucket_empty = false;
          }
        }
        if (!bucket_empty)
        {
          // initialize another container for edges
          // with neither vertex in hood
          vec_e pail;
          pail.clear();
          bool pail_empty = true;
          // for each edge in bucket
          for (e_iterator j=bucket.begin();j!=bucket.end();j++)
          {
            if (!processEdge(*j,hood,pail,*i))
            {
              pail_empty = false;
            }
          }
          if (!pail_empty)
          {
            cout << "Error. Multiple rounds of bucket use required.\n";
            exit(1);
          }
        }
        // add face to neighborhood, c
        c.push_back(*k);
      }
    }
    // copy local neighborhood to Object class neighborhood faces vector
    (*i)->nf.assign(c.begin(),c.end());
    // sort vectors
    sort((*i)->nf.begin(),(*i)->nf.end());
  }
}

void Object::findVertexAdjacencies (void)
{
  cerr << "Finding vertex adjacencies for object [" << name << "].................";
  cerr.flush();
  e_iterator i;
  f_iterator j;
  // for each edge, add edge* to both edge vertices
  for (i=e.begin();i!=e.end();i++)
  {
    Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
    (*i)->getVertices(v1,v2,o1,o2);
    v1->e.push_back(*i);
    v2->e.push_back(*i);
  }
  // for each face, add face* to each face vertex
  for (j=f.begin();j!=f.end();j++)
  {
    ((*j)->ptr_vertex(0))->f.push_back(*j);
    ((*j)->ptr_vertex(1))->f.push_back(*j);
    ((*j)->ptr_vertex(2))->f.push_back(*j);
  }
  cerr << "complete.\n";cerr.flush();
}

void Object::boundObject (double r[6]) 
{
  v_iterator i;
  double xmin,xmax,ymin,ymax,zmin,zmax,x,y,z;
  //initialize mins and maxes
  xmin = v[0]->getpN(0);
  xmax = v[0]->getpN(0);
  ymin = v[0]->getpN(1);
  ymax = v[0]->getpN(1);
  zmin = v[0]->getpN(2);
  zmax = v[0]->getpN(2);
  // for each vertex in object
  for (i=v.begin();i!=v.end();i++) 
  {
    ///////// extract coordinates //////////
    x = (*i)->getpN(0);
    y = (*i)->getpN(1);
    z = (*i)->getpN(2);
    if (x>xmax) {xmax = x;}
    else if (x<xmin) {xmin = x;}
    if (y>ymax) {ymax = y;}
    else if (y<ymin) {ymin = y;}
    if (z>zmax) {zmax = z;}
    else if (z<zmin) {zmin = z;}
  }
  // bounding box: [xmin,ymin,zmin,xmax,ymax,zmax]
  r[0]=xmin;r[3]=xmax;
  r[1]=ymin;r[4]=ymax;
  r[2]=zmin;r[5]=zmax;
}

void Object::evalAttributes (Space &s)
{
  cerr << "Checking if object [" << name << "] is closed..........................";
  cerr.flush();
  //closed=isClosed();
  cerr << "complete.\n";
  cerr.flush();
  cerr << "Checking if object [" << name << "] is manifold........................";
  cerr.flush();
  //manifold=isManifold();
  cerr << "complete.\n";cerr.flush();
  if (isManifold())
  {
    cerr << "Checking if object [" << name << "] faces are consistently oriented....";
    cerr.flush();
    //consistent=isConsistent();
    cerr << "complete.\n";cerr.flush();
    if (isConsistent() && isClosed())
    {
      cerr << "Checking if object [" << name << "] faces are oriented outward.........";
      cerr.flush();
      outward=isOutward(s);
      cerr << "complete.\n";cerr.flush();
    }
  }
}

bool Object::rayIntersectsSide (const char *str,double lp[2][3],double mybb[6],double n[3])
{
  double x1,y1,z1,a;
  if (!strcmp(str,"-x"))
  {
    x1=mybb[0];
    a=(x1-lp[0][0])/n[0];
    y1=a*n[1]+lp[0][1];
    z1=a*n[2]+lp[0][2];
    if (a>0 && y1>mybb[2] && y1<mybb[3] && z1>mybb[4] && z1<mybb[5] ){return true;}
  }
  else if (!strcmp(str,"+x"))
  {
    x1=mybb[1];
    a=(x1-lp[0][0])/n[0];
    y1=a*n[1]+lp[0][1];
    z1=a*n[2]+lp[0][2];
    if (a>0 && y1>mybb[2] && y1<mybb[3] && z1>mybb[4] && z1<mybb[5] ){return true;}
  }
  else if (!strcmp(str,"-y"))
  {
    y1=mybb[2];
    a=(y1-lp[0][1])/n[1];
    x1=a*n[0]+lp[0][0];
    z1=a*n[2]+lp[0][2];
    if (a>0 && x1>mybb[0] && x1<mybb[1] && z1>mybb[4] && z1<mybb[5] ){return true;}
  }
  else if (!strcmp(str,"+y"))
  {
    y1=mybb[3];
    a=(y1-lp[0][1])/n[1];
    x1=a*n[0]+lp[0][0];
    z1=a*n[2]+lp[0][2];
    if (a>0 && x1>mybb[0] && x1<mybb[1] && z1>mybb[4] && z1<mybb[5] ){return true;}
  }
  else if (!strcmp(str,"-z"))
  {
    z1=mybb[4];
    a=(z1-lp[0][2])/n[2];
    x1=a*n[0]+lp[0][0];
    y1=a*n[1]+lp[0][1];
    if (a>0 && x1>mybb[0] && x1<mybb[1] && y1>mybb[2] && y1<mybb[3] ){return true;}
  }
  else if (!strcmp(str,"+z"))
  {
    z1=mybb[5];
    a=(z1-lp[0][2])/n[2];
    x1=a*n[0]+lp[0][0];
    y1=a*n[1]+lp[0][1];
    if (a>0 && x1>mybb[0] && x1<mybb[1] && y1>mybb[2] && y1<mybb[3] ){return true;}
  }
  return false;
}

bool Object::rayIntersectsBB (double lp[2][3],Face *ff,double n[3])
{
  ///// compute face bounding box /////
  vec_d xv,yv,zv;
  // identify face bounding box limits
  xv.push_back(ff->ptr_vertex(0)->getpN(0));
  xv.push_back(ff->ptr_vertex(1)->getpN(0));
  xv.push_back(ff->ptr_vertex(2)->getpN(0));
  yv.push_back(ff->ptr_vertex(0)->getpN(1));
  yv.push_back(ff->ptr_vertex(1)->getpN(1));
  yv.push_back(ff->ptr_vertex(2)->getpN(1));
  zv.push_back(ff->ptr_vertex(0)->getpN(2));
  zv.push_back(ff->ptr_vertex(1)->getpN(2));
  zv.push_back(ff->ptr_vertex(2)->getpN(2));
  sort(xv.begin(),xv.end());
  sort(yv.begin(),yv.end());
  sort(zv.begin(),zv.end());
  // grab face 3D location range
  double mybb[6];
  mybb[0] = xv[0];  // -x
  mybb[1] = xv[2];  //  +x
  mybb[2] = yv[0];  // -y
  mybb[3] = yv[2];  //  +y
  mybb[4] = zv[0];  // -z
  mybb[5] = zv[2];  //  +z

  ///// for each side of bounding box /////
  if (rayIntersectsSide("-x",lp,mybb,n) || 
      rayIntersectsSide("+x",lp,mybb,n) || 
      rayIntersectsSide("-y",lp,mybb,n) || 
      rayIntersectsSide("+y",lp,mybb,n) || 
      rayIntersectsSide("-z",lp,mybb,n) || 
      rayIntersectsSide("+z",lp,mybb,n)
    ){return true;}
  // if ray intersects side, then return true
  return false;
}

bool Object::isOutward (Space &s)
{
  // assuming object mesh is closed, manifold, and consistent...
  //
  bool line_flag=false, poly_flag=false, poly_edge_flag=false;
  f_iterator fff=f.begin();
  double count;
  do
  {
    // set origin face
    Face *ff = *fff;
    count=0;
    // compute normal of origin face
    double n[3];
    ff->getNormal(n);
    double L=sqrt( dot(n,n) );
    n[0]=n[0]/L;
    n[1]=n[1]/L;
    n[2]=n[2]/L;
    // ray origin = centroid of origin face
    double lp[2][3];
    lp[0][0] = (ff->ptr_vertex(0)->getpN(0)+
                ff->ptr_vertex(1)->getpN(0)+
                ff->ptr_vertex(2)->getpN(0))/3.0;
    lp[0][1] = (ff->ptr_vertex(0)->getpN(1)+
                ff->ptr_vertex(1)->getpN(1)+
                ff->ptr_vertex(2)->getpN(1))/3.0;
    lp[0][2] = (ff->ptr_vertex(0)->getpN(2)+
                ff->ptr_vertex(1)->getpN(2)+
                ff->ptr_vertex(2)->getpN(2))/3.0;
    // ray end = point on normal advanced from origin
    //			 a distance equal to 2*world width along each axis
    double del[3]={fabs(s.getWorld(1)-s.getWorld(0)),fabs(s.getWorld(3)-s.getWorld(2)),fabs(s.getWorld(5)-s.getWorld(4))};
    int big;
    biggest(del,big);
    lp[1][0] = lp[0][0]+n[0]*2*del[big];
    lp[1][1] = lp[0][1]+n[1]*2*del[big];
    lp[1][2] = lp[0][2]+n[2]*2*del[big];
    // for each face in object
    for (f_iterator i=f.begin();i!=f.end();i++)
    {
      // if face in Object not same as origin face 
      if (*i!=ff)
      {
        // if ray intersects face bounding box
        if (rayIntersectsBB(lp,*i,n))
        {
          // if hit, determine face-line intersection
          checkLineFaceIntersection(*i,lp,line_flag,poly_flag,poly_edge_flag);
          // check
          if (poly_edge_flag==true) 
          {
            fff++;
            if (fff==f.end())
            {
              cout << "Object::isOutward: "
                    << "Error. Ray-tracing failed from all object <"
                    << name << "> faces.\n";
              exit(1);
            }
            break;
          }
          // does point intersect polygon
          if (poly_flag && line_flag) {count++;}
        }
      }
    }
  }while (poly_edge_flag==true);

  count = ceil(count);
  // if odd hit count then inward normal
  if (static_cast<int>(count)%2){return false;}
  // if even hit count then outward normal
  else {return true;}
}

bool Object::isConsistent (void)
{
  // assuming object mesh is manifold...
  bool flag=true;
  // for each edge in object
  for (e_iterator i=e.begin();i!=e.end();i++)
  {
    // if the one or two adjacent edge faces traverse the edge
    // in the same direction, either both v1->v2 or both v2->v1,
    // then the entire set of mesh faces is declared inconsistent
    // i.e. the first edge with a flipped face will trigger signal return
    // and no more edge checking is performed, but -p option
    // needs ALL flipped edges.
    if (!(*i)->isConsistent())
    {
      flag=false;
      // record offending edge
      flipped.push_back(*i);
    }
  }
  // update cumulative # flipped edges
  //	cs.num_flip=cs.flipped.size();
  return flag;
}

bool Object::isClosed (void)
{
  bool flag=true;
  // for each edge in object
  for (e_iterator i=e.begin();i!=e.end();i++)
  {
    if ((*i)->ptr_f2()==NULL)
    {
      flag=false;
      // record offending edge
      border.push_back(*i);
    }
  }
  // update cumulative # border edges
  //	cs.num_bor+=cs.border.size();
  return flag;
}

bool Object::verticesManifold (bool flag)
{
  ///// if vertices are manifold /////
  ////// confirm that all adjacent faces /////
  ///// are consecutively reachable by edge hops /////

  // sort orphan vertices
  sort(orphan.begin(),orphan.end());

  // for each vertex in object
  for (v_iterator i=v.begin();i!=v.end();i++)
  {
    // if vertex is not an orphan
    if (binary_search(orphan.begin(),orphan.end(),*i)==false)
    {
      flag = (*i)->isManifold(flag);
    }
  }
  return flag;
}

bool Object::edgesManifold (bool flag)
{
  // for each edge in object
  for (e_iterator i=e.begin();i!=e.end();i++)
  {
    // if edge is NOT manifold
    if ((*i)->isManifold()==false)
    {
      flag=false;
      // record offending edge
      nonman_e.push_back(*i);
    }
  }
  return flag;
}

bool Object::isManifold (void)
{
  bool flag=true;
  // check edge manifoldness
  flag=edgesManifold(flag);
  // check vertex manifoldness
  flag=verticesManifold(flag);

  // update cumulative # nonmanifold edges and vertices
  //	cs.num_nonman_e+=cs.nonman_e.size();
  //	cs.num_nonman_v+=cs.nonman_v.size();
  return flag;
}

void Object::checkIntegrity (void)
{
  // check vertex integrity
  cerr << "\nChecking vertex integrity for object [" << name << "]..................";
  cerr.flush();
  orphanMissingContig();
  cerr << "complete.\n";cerr.flush();
  // check face integrity
  cerr << "Checking face integrity for object [" << name << "]....................";
  cerr.flush();
  degenContig();
  cerr << "complete.\n";cerr.flush();
}

bool Object::goodIntegrity (void)
{
  // integrity is good if there are
  // (1) no missing vertices, face references nonexistent vertex index
  // (2) no orphan vertices, vertices referenced by no face
  // (3) no degenerate face, face references same vertex more than once
  // Note: missing data, such as a face with only two vertex indices given
  // 	would have been caught during input file parsing.
  return ((missing_v.empty()==true) &&
          (orphan.empty()==true) &&
          (degen.empty()==true) &&
          (dupl_v_index.empty()==true) &&
          (dupl_f_index.empty()==true));
}

// VERTICES
void Object::orphanMissingContig (void)
{
  // keep unique vertices
  sort(missing_v.begin(),missing_v.end());
  i_iterator new_end_v = unique(missing_v.begin(),missing_v.end());
  missing_v.assign(missing_v.begin(),new_end_v);
  // keep unique faces
  sort(missing_f.begin(),missing_f.end());
  f_iterator new_end_f = unique(missing_f.begin(),missing_f.end());
  missing_f.assign(missing_f.begin(),new_end_f);

  // update cumulative # missing vertices
  //	num_missing+=missing_v.size();

  ///// orphan ////
  std::pair<ib_iterator,ib_iterator> pp;
  std::pair<iv_iterator,iv_iterator> qq;
  // for each element of map #1
  ib_iterator i=found.begin();
  while (i!=found.end())
  {
    // if this vertex was not used in any face
    if ((*i).second==false)
    {
      // find all elements in map with this key
      pp=found.equal_range((*i).first);
      // if no matching elements were found
      if (pp.first==pp.second){ cout << "Error. Weird. how can this be?\n";exit(1);}
      else { // orphan vertex (or vertices)
        // find matching elements in map #2
        qq=vp.equal_range((*i).first);
        // if no matching elements were found
        if (qq.first==qq.second){ cout << "Error. Weird. how can this be too?\n";exit(1);}
        else 
        {
          // for each matching element
          for (iv_iterator n=qq.first;n!=qq.second;n++)
          {
            // add vertex* to orphan list
            orphan.push_back((*n).second);
          }
        }
        i=pp.second;
      }
    }
    else {i++;}
  }

  // update cumulative # orphan vertices
  //	num_orphan+=orphan.size();

  ///// contiguousness /////
  contig_v=true;
  // since map #1 is automatically sorted by integer value
  int pace = 1;
  // for each other element in map
  for (i=found.begin();i!=found.end();i++)
  {
    // if vertex indexing deviates 
    if ((*i).first!=pace++)
    {
      pace--;
      contig_v=false;
      // save five vertex indices to vec_cv
      if (pace==1)
      {
        vec_cv[0]=vec_cv[1]=-1;
        vec_cv[2]=(*i).first; i++;
        vec_cv[3]=(*i).first; i++;
        vec_cv[4]=(*i).first;
        i--;i--;
      }
      else if (pace==2) 
      {
        vec_cv[0]=-1;
        i--;
        vec_cv[1]=(*i).first; i++;
        vec_cv[2]=(*i).first; i++;
        vec_cv[3]=(*i).first; i++;
        vec_cv[4]=(*i).first;
        i--;i--;
      }
      else if (pace==static_cast<int>(v.size())-1) 
      {
        i--;i--;
        vec_cv[0]=(*i).first; i++;
        vec_cv[1]=(*i).first; i++;
        vec_cv[2]=(*i).first; i++;
        vec_cv[3]=(*i).first; i++;
        vec_cv[4]=-1;
        i--;i--;
      }
      else if (pace==static_cast<int>(v.size())) 
      {
        i--;i--;
        vec_cv[0]=(*i).first; i++;
        vec_cv[1]=(*i).first; i++;
        vec_cv[2]=(*i).first; i++;
        vec_cv[3]=vec_cv[4]=-1;
        i--;i--;
      }
      else 
      {
        i--;i--;
        vec_cv[0]=(*i).first; i++;
        vec_cv[1]=(*i).first; i++;
        vec_cv[2]=(*i).first; i++;
        vec_cv[3]=(*i).first; i++;
        vec_cv[4]=(*i).first;
        i--;i--;
      }
      break;
    }
  }

  ///// contiguousness /////
  // for each element in multimap(int->Vertex*)
  for (iv_iterator k=vp.begin();k!=vp.end();k++)
  {
    // find # elements in map with this key
    if (vp.count((*k).first)>1)
    {
      dupl_v_index.push_back((*k).second);
    }
  }
}

void Object::vertexAdjacentFaces (void)
{
  //	mmap_iv af;
  // for each vertex in object
  for (v_iterator i=v.begin();i!=v.end();i++)
  {
    int c=(*i)->f.size();
    //		af.insert(std::make_pair(c,*i));
    //		adjacent_face.n++;
    adjacent_face.add2sum(c);
    adjacent_face.add2sum2(c*c);
    adjacent_face.add2total(c);
    if (c<adjacent_face.getMin()) {adjacent_face.setMin(c);}
    if (c>adjacent_face.getMax()) {adjacent_face.setMax(c);}
    // add to vector
    adjacent_face.insertElement(c);
  }
  // build adjacent face histogram
  adjacent_face.createAdjacentFaceHistogram();
}

// FACES
void Object::areaAspectRatio (Controls &cs)
{
  // for each face in object
  for (f_iterator i=f.begin();i!=f.end();i++)
  {
    //		area.n++;
    //		aspect_ratio.n++;
    ///// face area /////
    // compute face normal vector
    double n[3];
    (*i)->getNormal(n);
    // compute face area = half normal vector length
    double aa=sqrt(dot(n,n))/2.0;
    //
    area.add2sum(aa);
    area.add2sum2(aa*aa);
    area.add2total(aa);
    if (aa<area.getMin()) {area.setMin(aa);}
    if (aa>area.getMax()) {area.setMax(aa);}
    // add to vector
    area.insertElement(aa);
    ///// aspect ratio /////
    double ar = (*i)->getAspectRatio();
    aspect_ratio.add2sum(ar);
    aspect_ratio.add2sum2(ar*ar);
    aspect_ratio.add2total(ar);
    if (ar<aspect_ratio.getMin()) aspect_ratio.setMin(ar);
    if (ar>aspect_ratio.getMax()) aspect_ratio.setMax(ar);
    if (cs.signal[0]==true)
    {
      // compare face aspect ratio to user-defined threshold
      if (ar>cs.thresholds[0]) bad_aspect[*i]=ar;
    }
    // add to vector
    aspect_ratio.insertElement(ar);
  }
  assert (area.getSize() > 0);
  // update cumulative surface area
  //	a+=area.sum;
  // build face area histogram
  area.createHistogram();
  // build aspect ratio histogram
  aspect_ratio.createAspectRatioHistogram();
}

void Object::degenContig (void)
{
  // buide map: face index (integer)->face*
  mmap_if tt;
  // for each face in object
  for (f_iterator i=f.begin();i!=f.end();i++)
  {
    tt.insert(std::make_pair((*i)->getIndex(),*i));
  }
  contig_f=true;
  int pace=1;
  // for each face in map
  for (if_iterator i=tt.begin();i!=tt.end();i++)
  {
    Face *ff=(*i).second;
    ///// duplicity /////
    // find # elements in map with this key
    if (tt.count((*i).first)>1)
    {
      dupl_f_index.push_back(ff);
    }
    ///// degeneracy /////
    // if any two of the vertex indices are the same
    if ( ff->ptr_vertex(0)==ff->ptr_vertex(1) || ff->ptr_vertex(0)==ff->ptr_vertex(2) || ff->ptr_vertex(1)==ff->ptr_vertex(2))
    {
      // then face is degenerate
      degen.push_back(ff);
      //			num_degen++;
    }
    // update cumulative # degenerate faces
    ///// contiguousness /////
    if (contig_f)
    {
      if ((*i).first!=pace++)
      {
        // face indexing NOT contiguous
        contig_f=false;
        pace--;
        // save five face indices to vec_cf
        if (pace==1)
        {
          vec_cf[0]=vec_cf[1]=-1;
          vec_cf[2]=(*i).first; i++;
          vec_cf[3]=(*i).first; i++;
          vec_cf[4]=(*i).first;
          i--;i--;
        }
        else if (pace==2) 
        {
          vec_cf[0]=-1;
          i--;
          vec_cf[1]=(*i).first; i++;
          vec_cf[2]=(*i).first; i++;
          vec_cf[3]=(*i).first; i++;
          vec_cf[4]=(*i).first;
          i--;i--;
        }
        else if (pace==static_cast<int>(f.size())-1) 
        {
          i--;i--;
          vec_cf[0]=(*i).first; i++;
          vec_cf[1]=(*i).first; i++;
          vec_cf[2]=(*i).first; i++;
          vec_cf[3]=(*i).first; i++;
          vec_cf[4]=-1;
          i--;i--;
        }
        else if (pace==static_cast<int>(f.size())) 
        {
          i--;i--;
          vec_cf[0]=(*i).first; i++;
          vec_cf[1]=(*i).first; i++;
          vec_cf[2]=(*i).first; i++;
          vec_cf[3]=vec_cf[4]=-1;
          i--;i--;
        }
        else 
        {
          i--;i--;
          vec_cf[0]=(*i).first; i++;
          vec_cf[1]=(*i).first; i++;
          vec_cf[2]=(*i).first; i++;
          vec_cf[3]=(*i).first; i++;
          vec_cf[4]=(*i).first;
          i--;i--;
        }
      }
    }
  }
}

void Object::setAll (Vertex *vv,hmap_fi &group,int &g)
{
  // for each vertex adjacent face
  for (f_iterator i=vv->f.begin();i!=vv->f.end();i++)
  {
    group[*i]=g;
  }
  g++;
}

void Object::getGroups (Vertex *vv,hmap_fi &group,set_i &s)
{
  s.clear();
  // for each vertex adjacent face
  for (f_iterator i=vv->f.begin();i!=vv->f.end();i++)
  {
    s.insert(group[*i]);
  }
}

int Object::getLowest (set_i &s)
{
  int i = 100000;
  // for each element in set
  for (std::set<int,lti>::iterator j=s.begin();j!=s.end();j++)
  {
    // if element is lower than i and not 0
    if (*j<i && *j){i=*j;}
  }
  return i;
}

void Object::replaceGroups (Vertex *vv,hmap_fi &group,int z)
{
  // for each vertex adjacent face
  for (f_iterator i=vv->f.begin();i!=vv->f.end();i++)
  {
    group[*i]=z;
  }
}

void Object::setZero (Edge *ee,hmap_fi &group,int z)
{
  fi_iterator i=group.find(ee->ptr_f1());
  if ((*i).second==0){group[ee->ptr_f1()]=z;}
  // if second adjacent face
  if (ee->ptr_f2()!=NULL)
  {
    i=group.find(ee->ptr_f2());
    if ((*i).second==0){group[ee->ptr_f2()]=z;}
  }
  // if more adjacent faces
  if (!ee->noExtraFaces())
  {
    // for each adjacent face
    for (c_f_iterator j=ee->first_extra_face();j!=ee->one_past_last_extra_face();j++)
    {
      // if adjacent face has no group
      i=group.find(*j);
      if ((*i).second==0){group[*j]=z;}
    }
  }
}

bool Object::removeSelectedFaces (vec_f &sf,vec_f &fv)
{
  bool flag=false;
  // for each selected face
  f_iterator i=sf.begin();
  while (i!=sf.end())
  {
    f_iterator j= find(fv.begin(),fv.end(),*i);
    // if selected face found in face vector
    if (j!=fv.end())
    {
      // remove face from face vector
      fv.erase(j);
      flag=true;
      i++;
    }
    else 
    {
      // remove face from selected face vector
      //  so that the face is not used next round
      sf.erase(i);
    }
  }
  return flag;
}

void Object::getSelectedFaces (vec_f &sf,vec_f fv)
{
  // if selected faces vector is empty
  if (sf.empty()==true)
  {
    // then grab first face from vector
    sf.push_back(fv.front());
  }
  else 
  {
    // adjacent faces
    vec_f af;
    // for each selected face
    for (f_iterator i=sf.begin();i!=sf.end();i++)
    {
      // grab adjacent faces of selected face
      if ((*i)->ptr_edge(0)->ptr_f1()!=*i) { af.push_back((*i)->ptr_edge(0)->ptr_f1()); }
      else                    { af.push_back((*i)->ptr_edge(0)->ptr_f2()); }
      if ((*i)->ptr_edge(1)->ptr_f1()!=*i) { af.push_back((*i)->ptr_edge(1)->ptr_f1()); }
      else                    { af.push_back((*i)->ptr_edge(1)->ptr_f2()); }
      if ((*i)->ptr_edge(2)->ptr_f1()!=*i) { af.push_back((*i)->ptr_edge(2)->ptr_f1()); }
      else                    { af.push_back((*i)->ptr_edge(2)->ptr_f2()); }
    }
    // keep unique faces
    sort(af.begin(),af.end());
    f_iterator i=unique(af.begin(),af.end());
    sf.assign(af.begin(),i);
  }
}

int Object::countComponents (void)
{
  // initialize group number
  int g=0;
  ///// create hashed map face*->integer (group #) /////
  hmap_fi group;
  // for each face
  for (f_iterator i=f.begin();i!=f.end();i++)
  {
    group[*i]=g;
  }
  // increment group number
  g++;
  // set of integers
  set_i s;
  bool changes = true;
  while (changes==true)
  {
    changes=false;
    // for each vertex in object
    for (v_iterator i=v.begin();i!=v.end();i++)
    {
      // get vertex adjacent face groups
      getGroups(*i,group,s);
      std::set<int,lti>::iterator j=s.begin();
      // if no adjacent face group has been set
      if (s.size()==1 && *j==0)
      {
        // set all adjacent faces to next available group #
        setAll(*i,group,g);
        changes=true;
      }
      else if (s.size()>1)
      { 
        // more than one group
        // identify lowest group # larger than 0 in set
        int z=getLowest(s);
        // replace the larger group # with lowest group # in all faces
        replaceGroups(*i,group,z);
        changes=true;
      } 
    }
  }
  // analyze map
  s.clear();
  // for each element in group
  for (fi_iterator i=group.begin();i!=group.end();i++)
  {
    s.insert((*i).second);
  }
  return s.size();
}

void Object::findIntersectingFaces (Container *c,Space &s)
{
  vec_f dummy;
  // for each box
  for (c_b_iterator i=s.first_box();i!=s.one_past_last_box();i++)
  {
    // if box is not empty
    if ((*i)->noFaces()==false)
    {
      // check for intersections
      (*i)->getFaceIntersection(c);
      // intersecting faces end up in o->intf (hashtable_f_face)
    }
  }
}

// EDGES
void Object::countBoundaries (void)
{
  num_bou=0;
  Boundary boundary;
  // keep unique border edges
  sort(border.begin(),border.end());
  e_iterator new_end = unique(border.begin(),border.end());
  border.assign(border.begin(),new_end);
  // copy the border vector to ws
  vec_e ws;
  ws.assign(border.begin(),border.end());
  // while there are border edges in ws
  while (ws.empty()==false)
  {
    // if the boundary is closed
    if (boundary.isOpen()==false)
    {
      // init boundary
      boundary.init(ws.back());
      // pop last edge from ws
      ws.pop_back();
    }
    bool myfound = false;
    // for each edge in ws
    e_iterator i=ws.begin();
    while (i!=ws.end())
    {
      // if edge extends boundary
      if (boundary.edgeExtendsBoundary(*i))
      {
        myfound = true;
        // remove edge from ws
        ws.erase(i);
        // if boundary is closed
        if (boundary.closed())
        {
          // increment boundary count
          num_bou++;
          // break
          break;
        }
      }
      else { i++;}
    }
    // if no edges added to boundary
    // then assume boundary is due to dangling face
    if (myfound==false)
    {
      boundary.setOpen(false);
      num_bou++;
    }
  }
  // update cumulative number boundaries
  //	cs.num_bou+=num_bou;
}

void Object::processEdgeLengths (Controls &cs)
{
  // for each edge in object
  for (e_iterator i=e.begin();i!=e.end();i++)
  {
    double l = (*i)->getOrigLength();
    //		edge_length.n++;
    edge_length.add2sum(l);
    edge_length.add2sum2(l*l);
    edge_length.add2total(l);
    if (l<edge_length.getMin()) edge_length.setMin(l);
    if (l>edge_length.getMax()) edge_length.setMax(l);
    if (cs.signal[3]==true)
    {
      if (l<cs.thresholds[3]) bad_length[*i]=l;
    }
    if (cs.signal[4]==true)
    {
      if (l>cs.thresholds[4]) bad_length[*i]=l;
    }
    // check distinguishability
    Vertex *v1=NULL,*v2=NULL,*o1=NULL,*o2=NULL;
    (*i)->getVertices(v1,v2,o1,o2);
    if (!distinguishable(v1->getpN(0),v2->getpN(0)) &&
        !distinguishable(v1->getpN(1),v2->getpN(1)) &&
        !distinguishable(v1->getpN(2),v2->getpN(2)) 
      ){ indistin.push_back(*i); }
    // add to vector
    edge_length.insertElement(l);
  }
  edge_length.createHistogram();
}

void Object::computeGenus (void)
{
  int num = v.size()-e.size()+f.size();
  if (num%2) 
  {
    cout << "Round off error in genus computation. "
          << "(#vertices + #faces - #edges) is not evenly divisible by 2.\n";
    cout << "v-e+f = " << num << endl;
    cout << v.size() << "-" << e.size() << "+" << f.size() << endl;
    exit(1);
  }
  else 
  {
    genus = num_sep-num/2;
  }
}

void Object::computeEdgeAngles (void)
{
  Controls & cs(Controls::instance());
  // for each edge in object
  for (e_iterator i=e.begin();i!=e.end();i++)
  {
    // if edge has exactly two adjacent faces
    if ((*i)->ptr_f2()!=NULL && (*i)->noExtraFaces())
    {
      double angle = (*i)->getAngle()*180/cs.get_pi(); // degrees
      //			edge_angle.n++;
      edge_angle.add2sum(angle);
      edge_angle.add2sum2(angle*angle);
      edge_angle.add2total(angle);
      if (angle<edge_angle.getMin()) edge_angle.setMin(angle);
      if (angle>edge_angle.getMax()) edge_angle.setMax(angle);
      if (cs.signal[1]==true)
      {
        if (angle<cs.thresholds[1]) bad_angle[*i]=angle;
      }
      if (cs.signal[2]==true)
      {
        if (angle>cs.thresholds[2]) bad_angle[*i]=angle;
      }
      // add to vector
      edge_angle.insertElement(angle);
    }
  }
  edge_angle.createHistogram();
}

void Object::computeVolume (void)
{
  vol = 0.0;
  // for each face in object
  for (f_iterator i=f.begin();i!=f.end();i++)
  {
    double x1=(*i)->ptr_vertex(0)->getpN(0);
    double y1=(*i)->ptr_vertex(0)->getpN(1);
    double z1=(*i)->ptr_vertex(0)->getpN(2);
    double x2=(*i)->ptr_vertex(1)->getpN(0);
    double y2=(*i)->ptr_vertex(1)->getpN(1);
    double z2=(*i)->ptr_vertex(1)->getpN(2);
    double x3=(*i)->ptr_vertex(2)->getpN(0);
    double y3=(*i)->ptr_vertex(2)->getpN(1);
    double z3=(*i)->ptr_vertex(2)->getpN(2);
    /* compute determinant of oriented triangle */
    double det=x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1);
    vol+=det;
  }
  vol=vol/6.0;
  // update cumulative volume
  //	vol+=o->vol;
}

void Object::vertexDistin (void)
{
  ///// check vertex distinguishability /////
  // multimap: double -> Vertex*
  mmap_dv mm;
  // random vector
  double rand_vec[3]={0.236416584579274058342,
    0.927225593011826276779,
    0.389099507126957178116};
  // for each vertex
  for (v_iterator i=v.begin();i!=v.end();i++)
  {
    // add dot product of vertex and random vector to multimap
    mm.insert(std::make_pair(dot((*i)->getpN_ptr(),rand_vec),*i));
  }
  // for each multimap element
  dv_iterator j;
  for (dv_iterator i=mm.begin();i!=mm.end();i++)
  {
    j=i;j++;
    if (j!=mm.end())
    {
      // if sequential pair of multimap elements is not distinguishable
      if (!distinguishable((*i).second->getpN(0),(*j).second->getpN(0)) &&
          !distinguishable((*i).second->getpN(1),(*j).second->getpN(1)) &&
          !distinguishable((*i).second->getpN(2),(*j).second->getpN(2)) )
      {
        indistin_v.push_back((*i).second);
        indistin_v.push_back((*j).second);
      }
    }
  }
}

void Object::evalCharacteristics (Container* c,Controls &cs,Space &s)
{
  // vertices
  cerr << "Bound object [" << name << "]..........................................";
  cerr.flush();
  boundObject(bb);
  cerr << "complete.\n";cerr.flush();
  cerr << "Check if vertices are distinguishable for object [" << name << "]......";
  cerr.flush();
  vertexDistin();
  cerr << "complete.\n";cerr.flush();
  cerr << "Analyze vertex adjacent faces for object [" << name << "]..............";
  cerr.flush();
  vertexAdjacentFaces();
  cerr << "complete.\n";cerr.flush();
  // faces
  cerr << "Identify separate components of object [" << name << "]................";
  cerr.flush();
  num_sep=countComponents();
  //	cs.num_sep+=num_sep;
  cerr << "complete.\n";cerr.flush();
  cerr << "Compute face area and aspect ratio for object [" << name << "].........";
  cerr.flush();
  areaAspectRatio(cs);
  cerr << "complete.\n";cerr.flush();
  // if not batch mode
  if (cs.get_detect_interobject_intersections()==0)
  {
    cerr << "Find intersecting faces for object [" << name << "]....................";
    cerr.flush();
    findIntersectingFaces(c,s);
    cerr << "complete.\n";cerr.flush();
  }
  // edges
  cerr << "Identify boundaries for object [" << name << "]........................";
  cerr.flush();
  countBoundaries();
  cerr << "complete.\n";cerr.flush();
  cerr << "Analyze edge lengths for object [" << name << "].......................";
  cerr.flush();
  processEdgeLengths(cs);
  cerr << "complete.\n";cerr.flush();

  if (isManifold() && isConsistent())
  {
    // manifold and consistently oriented face normals
    cerr << "Analyze edge angles for object [" << name << "]........................";
    cerr.flush();
    computeEdgeAngles();
    cerr << "complete.\n";cerr.flush();
    if (isClosed())
    {
      // closed, manifold, and consistently oriented face normals
      cerr << "Compute volume of object [" << name << "]..............................";
      cerr.flush();
      computeVolume();
      cerr << "complete.\n";cerr.flush();
      // if number of components == 1
      // FUTURE IMPROVEMENT: only require orientable, not oriented
      // FUTURE IMPROVEMENT: separate components and compute genus of each component
      if (num_sep==1 && orphan.empty())
      {
        cerr << "Compute genus of object [" << name << "]...............................";
        cerr.flush();
        computeGenus();
        cerr << "complete.\n";cerr.flush();
      }
    }
  }
}

void Object::analyze (Container *c,Controls &cs,Space &s)
{
  // single object
  // evaluate mesh attributes
  evalAttributes(s);
  // eval mesh characteristics
  if (cs.get_compute_attributes_only()==0)
  {
    evalCharacteristics(c,cs,s);
    assert(area.getSize() > 0);
  }
}

void Object::printChars (Controls &cs)
{
  cout << "\nMESH CHARACTERISTICS\n\n";
  cout << "    # vertices: " << v.size() << endl
        << "    # faces: " << f.size() << endl
        << "    # edges: " << e.size() << endl
        << "    # components: " << num_sep << endl;
  if (isManifold()==false)
  {
    cout << "    # boundaries: Since object is nonmanifold,\n"
          << "    # boundaries: the number of boundaries may be underestimated.\n";
  }
  //////////////////// borders
  if (border.empty())
  {
    cout << "    # boundaries: none\n";
  }
  else 
  {
    cout << "    # boundaries: " << border.size() << endl;
    //	if -p option, print offending
    if (cs.get_print_detailed_info()==1)
    {
      int j=1;
      // for each border edge
      for (e_iterator i=border.begin();i!=border.end();i++)
      {
        cout << "    # boundaries: boundary edge " << j++ << endl;
        if (cs.get_use_dreamm_output_format())
        {
          (*i)->printCP(cout);
        }
        else 
        {
          (*i)->print(cout);
          cout << endl;
        }
      }
    }
  }
  ///////////////////// indistinguishable vertices
  if (indistin_v.empty())
  {
    cout << "    # indistinguishable vertices: none\n";
  }
  else 
  {
    cout << "    # indistinguishable vertices: " << indistin_v.size() << endl;
    //	if -p option, print offending
    if (cs.get_print_detailed_info()==1)
    {
      int j=1;
      // for each indistinguishable vertec
      for (v_iterator i=indistin_v.begin();i!=indistin_v.end();i++)
      {
        cout << "    #  indistinguishable vertices: vertex " << j++ << endl;
        if (cs.get_use_dreamm_output_format())
        {
          (*i)->printCP(cout);
        }
        else 
        {
          (*i)->print(cout);
          cout << endl;
        }
      }
    }
  }
  //////////////// total area
  cout << "    object area: [(data units)^2]" << endl;
  cout << "    object area: " << area.getSum() << endl;
  //////////////// if volume computed
  if (isClosed()==true && isManifold()==true && isConsistent()==true)
  {
    cout << "    object volume: [(data units)^3]" << endl;
    cout << "    object volume: " << vol << endl;
  }
  else 
  {
    cout << "    object volume: not computed, since ";
    if (isClosed()==false){cout << "not closed,";}
    if (isConsistent()==false){cout << "not consistent,";}
    if (isManifold()==false){cout << "not manifold";}
    cout << endl;
  }
  /////////////// if genus computed
  if (isClosed()==true && isManifold()==true && isConsistent()==true && num_sep==1 && orphan.empty())
  {
    cout << "    object genus: " << genus << endl;
  }
  else 
  {
    cout << "    object genus: not computed, since ";
    if (isClosed()==false){cout << "not closed,";}
    if (isConsistent()==false){cout << "not consistent,";}
    if (isManifold()==false){cout << "not manifold,";}
    if (num_sep>1){cout << "#components=" << num_sep << ",";}
    if (!orphan.empty()){cout << "orphan vertices were found";}
    cout << endl;
  }
  //////////////// bounding box
  cout << "    bounding box: [data units]\n";
  cout << "    bounding box: [xmin,ymin,zmin][xmax,ymax,zmax]\n";
  cout << "    bounding box: ["
        << bb[0] << ","
        << bb[1] << ","
        << bb[2] << "]["
        << bb[3] << ","
        << bb[4] << ","
        << bb[5] << "]" << endl;
  ////////////////// edges with indistinguishable vertices
  if (indistin.empty()==true)
  {
    cout << "    # edges with indistinguishable vertices: none\n";
  }
  else 
  {
    cout << "    # edges with indistinguishable vertices: "
          << indistin.size() << endl;
  }
  //	if -p option, print offending
  if (cs.get_print_detailed_info()==1)
  {
    if (indistin.empty()==false)
    {
      // for each afflicted edge
      for (e_iterator i=indistin.begin();i!=indistin.end();i++)
      {
        if (cs.get_use_dreamm_output_format())
        {
          (*i)->printCP(cout);
        }
        else 
        {
          (*i)->print(cout);
          cout << endl;
        }
      }
    }
  }
  /////////////// intersecting faces
  // if not batch mode
  if (cs.get_detect_interobject_intersections()==0)
  {
    // intersecting faces
    if (intf.empty())
    {
      cout << "    # intersecting faces: none\n\n";
    }
    else 
    {
      cout << "    # intersecting faces: " << intf.size() << endl;
      //	if -p option, print offending
      if (cs.get_print_detailed_info()==1)
      {
        int j=1;
        // for each intersected face
        for (ff_iterator i=intf.begin();i!=intf.end();i++)
        {
          cout << "    # intersecting faces: intersected face " << j++ << endl;
          // print intersected face
          if (cs.get_use_dreamm_output_format())
          {
            (*i).first->printCP(cout);
          }
          else 
          {
            (*i).first->print(cout);
          }
          // keep unique list of intersecting faces
          sort((*(*i).second).begin(),(*(*i).second).end());
          f_iterator new_end = unique((*(*i).second).begin(),(*(*i).second).end());
          (*(*i).second).assign((*(*i).second).begin(),new_end);
          // print intersecting faces
          for (f_iterator k=(*(*i).second).begin();k!=(*(*i).second).end();k++)
          {
            if (cs.get_use_dreamm_output_format())
            {
              (*k)->printCP(cout);
            }
            else 
            {
              (*k)->print(cout);
              cout << endl;
            }
          }
        }
        cout << endl << endl;
      }
    }
  }
  //////////////// vertex adjacent faces
  cout << "    Vertex adjacent face statistics [faces]:" << endl;
  adjacent_face.printStats();
  cout << "    Vertex adjacent face histogram [faces]:" << endl;
  adjacent_face.printAdjacentFaceHistogram();
  cout << endl;
  ///////////////// face area
  cout << "    Face area statistics [(data units)^2]:" << endl;
  //cout << "       total    " << area.getSum() << endl;
  area.printStats();
  cout << "    Face area histogram [(data units)^2]:" << endl;
  area.printHistogram();
  cout << endl;
  /////////////////// face aspect ratio
  cout << "    Face aspect ratio statistics [unitless]:" << endl;
  aspect_ratio.printStats();
  cout << "    Face aspect ratio histogram [unitless]:" << endl;
  aspect_ratio.printAspectRatioHistogram();
  cout << "      (Aspect ratio is longest edge "
        << "divided by shortest altitude)\n";
  // if face aspect ratio threshold specified
  if (cs.signal[0]==true)
  {
    if (bad_aspect.empty()==true)
    {
      cout << "    # faces with bad aspect ratios: none" << endl;
    }
    else 
    {
      cout << "    # faces with bad aspect ratios: "
            << bad_aspect.size() << endl;
    }
  }
  //	if -p option, print offending
  if (cs.get_print_detailed_info()==1)
  {
    if (bad_aspect.empty()==false)
    {
      // print faces with aspect ratios that violate threshold
      Face_Pair fp(this);
      for (fd_iterator i=bad_aspect.begin();i!=bad_aspect.end();i++)
      {
        cout << "    face aspect ratio: " << (*i).second << endl;
        if (cs.get_use_dreamm_output_format())
        {
          (*i).first->printCP(cout);
        }
        else 
        {
          (*i).first->print(cout);
          cout << endl;
        }
        fp.processBadFace((*i).first);
      }
    }
  }
  cout << "\n";
  ////////////// edge length
  cout << "    Edge length statistics [data units]:" << endl;
  edge_length.printStats();
  cout << "    Edge length histogram [data units]:" << endl;
  edge_length.printHistogram();
  cout << endl;
  // if edge length threshold specified
  if (cs.signal[3]==true || cs.signal[4]==true)
  {
    if (bad_length.empty()==true)
    {
      cout << "    # edges with bad lengths: none" << endl;
    }
    else 
    {
      cout << "    # edges with bad lengths: "
            << bad_length.size() << endl;
    }
  }
  //	if -p option, print offending
  if (cs.get_print_detailed_info()==1)
  {
    if (bad_length.empty()==false)
    {
      // print edges with lengthss that violate threshold
      for (ed_iterator i=bad_length.begin();i!=bad_length.end();i++)
      {
        cout << "    edge length: " << (*i).second << endl;
        if (cs.get_use_dreamm_output_format())
        {
          (*i).first->printCP(cout);
        }
        else 
        {
          (*i).first->print(cout);
          cout << endl;
        }
      }
    }
  }
  //////////////////// if edge angles computed
  if (isManifold()==true && isConsistent()==true)
  {
    cout << "    Edge angle statistics [degress]:" << endl;
    edge_angle.printStats();
    cout << "    Edge angle histogram [degress]:" << endl;
    edge_angle.printHistogram();
    cout << endl;
    // if edge angle threshold specified
    if (cs.signal[1]==true || cs.signal[2]==true)
    {
      if (bad_angle.empty()==true)
      {
        cout << "    # edges with bad angles: none" << endl;
      }
      else 
      {
        cout << "    # edges with bad angles: "
              << bad_angle.size() << endl;
      }
    }
    //	if -p option, print offending
    if (cs.get_print_detailed_info()==1)
    {
      if (bad_angle.empty()==false)
      {
        // print edges with angles that violate threshold
        for (ed_iterator i=bad_angle.begin();i!=bad_angle.end();i++)
        {
          cout << "    edge angle: " << (*i).second << endl;
          if (cs.get_use_dreamm_output_format())
          {
            (*i).first->printCP(cout);
          }
          else 
          {
            (*i).first->print(cout);
            cout << endl;
          }
        }
      }
    }
  }
  else 
  {
    cout << "    edge angles: not computed, since ";
    if (isConsistent()==false){cout << "not consistent,";}
    if (isManifold()==false){cout << "not manifold";}
    cout << endl;
  }

}

void Object::printIntegrity (Controls &cs)
{
  cout << "\nMESH FILE INTEGRITY\n\n";
  // orphan vertices
  if (orphan.empty())
  {
    cout << "    # orphan vertices: none\n";
  }
  else 
  {
    cout << "    # orphan vertices: " << orphan.size() << endl;
    //	if -p option, print offending
    if (cs.get_print_detailed_info()==1)
    {
      int j=1;
      // for each orphan vertex
      for (v_iterator i=orphan.begin();i!=orphan.end();i++)
      {
        cout << "    # orphan vertices: orphan vertex " << j++ << endl;
        if (cs.get_use_dreamm_output_format())
        {
          (*i)->printCP(cout);
        }
        else 
        {
          (*i)->print(cout);
          cout << endl;
        }
      }
    }
  }
  // missing vertices
  if (missing_v.empty())
  {
    cout << "    # missing vertices: none\n";
  }
  else 
  {
    cout << "    # missing vertices: " << missing_v.size() << endl;
    //	if -p option, print offending
    if (cs.get_print_detailed_info()==1)
    {
      int j=1;
      // for each missing vertex
      for (i_iterator i=missing_v.begin();i!=missing_v.end();i++)
      {
        cout << "    # missing vertices: #" << j++
              << "-> missing vertex index " << *i << endl;
      }
      j=1;
      // for each missing face
      for (f_iterator i=missing_f.begin();i!=missing_f.end();i++)
      {
        cout << "    # missing vertices: affected face " << j++ << endl;
        if ((*i)->ptr_vertex(0)!=NULL)
        {
          if (cs.get_use_dreamm_output_format())
          {
            (*i)->printCP(cout);
          }
          else 
          {
            (*i)->print(cout);
            cout << endl;
          }
        }
        else if ((*i)->ptr_vertex(1)!=NULL)
        {
          if (cs.get_use_dreamm_output_format())
          {
            (*i)->printCP(cout);
          }
          else 
          {
            (*i)->print(cout);
            cout << endl;
          }
        }
        else if ((*i)->ptr_vertex(2)!=NULL)
        {
          if (cs.get_use_dreamm_output_format())
          {
            (*i)->printCP(cout);
          }
          else 
          {
            (*i)->print(cout);
            cout << endl;
          }
        }
      }
    }
  }
  // degenerate faces
  if (degen.empty())
  {
    cout << "    # degenerate faces: none\n";
  }
  else 
  {
    cout << "    # degenerate faces: " << degen.size() << endl;
    //	if -p option, print offending
    if (cs.get_print_detailed_info()==1)
    {
      int j=1;
      // for each degenerate face
      for (f_iterator i=degen.begin();i!=degen.end();i++)
      {
        cout << "    # degenerate faces: affected face " << j++ << endl;
        if (cs.get_use_dreamm_output_format())
        {
          (*i)->printCP(cout);
        }
        else 
        {
          (*i)->print(cout);
          cout << endl;
        }
      }
    }
  }

  // duplicate indices
  if (dupl_v_index.empty())
  {
    cout << "    # duplicate vertex indices: none\n";
  }
  else 
  {
    cout << "    # duplicat vertex indices: " << dupl_v_index.size() << endl;
    //	if -p option, print offending
    if (cs.get_print_detailed_info()==1)
    {
      int j=1;
      // for each vertex with a duplicate vertex index
      for (v_iterator i=dupl_v_index.begin();i!=dupl_v_index.end();i++)
      {
        cout << "    # duplicate vertex indices: affected vertex " << j++ << endl;
        if (cs.get_use_dreamm_output_format())
        {
          (*i)->printCP(cout);
        }
        else 
        {
          (*i)->print(cout);
          cout << endl;
        }
      }
    }
  }
  if (dupl_f_index.empty())
  {
    cout << "    # duplicate face indices: none\n";
  }
  else 
  {
    cout << "    # duplicat face indices: " << dupl_f_index.size() << endl;
    //	if -p option, print offending
    if (cs.get_print_detailed_info()==1)
    {
      int j=1;
      // for each face with a duplicate face index
      for (f_iterator i=dupl_f_index.begin();i!=dupl_f_index.end();i++)
      {
        cout << "    # duplicate face indices: affected face " << j++ << endl;
        if (cs.get_use_dreamm_output_format())
        {
          (*i)->printCP(cout);
        }
        else 
        {
          (*i)->print(cout);
          cout << endl;
        }
      }
    }
  }

  // contiguous numbering
  if (contig_v)
  {
    cout << "    contiguous vertex indexing from 1: yes\n";
  }
  else 
  {
    cout << "    contiguous vertex indexing from 1: no\n";
    cout << "    contiguous vertex indexing from 1: bad index and +/- 2\n";
    cout << "    contiguous vertex indexing from 1: ";
    if (vec_cv[0]!=-1){ cout << vec_cv[0] << " ";}
    else 			{cout << "NA ";}
    if (vec_cv[1]!=-1){ cout << vec_cv[1] << " ";}
    else 			{cout << "NA ";}
    if (vec_cv[2]!=-1){ cout << vec_cv[2] << " ";}
    else 			{cout << "NA ";}
    if (vec_cv[3]!=-1){ cout << vec_cv[3] << " ";}
    else 			{cout << "NA ";}
    if (vec_cv[4]!=-1){ cout << vec_cv[4] << endl;}
    else 			{cout << "NA\n";}
  }
  if (contig_f)
  {
    cout << "    contiguous face indexing from 1: yes\n";
  }
  else 
  {
    cout << "    contiguous face indexing from 1: no\n";
    cout << "    contiguous face indexing from 1: bad index and +/- 2\n";
    cout << "    contiguous face indexing from 1: ";
    if (vec_cf[0]!=-1){ cout << vec_cf[0] << " ";}
    else 			{cout << "NA ";}
    if (vec_cf[1]!=-1){ cout << vec_cf[1] << " ";}
    else 			{cout << "NA ";}
    if (vec_cf[2]!=-1){ cout << vec_cf[2] << " ";}
    else 			{cout << "NA ";}
    if (vec_cf[3]!=-1){ cout << vec_cf[3] << " ";}
    else 			{cout << "NA ";}
    if (vec_cf[4]!=-1){ cout << vec_cf[4] << endl;}
    else 			{cout << "NA\n";}
  }
}

void Object::printAttr (Controls &cs)
{
  cout << "\nMESH ATTRIBUTES\n\n";
  // closed
  if (isClosed()==true)
  {
    cout << "    mesh is closed: yes\n";
  }
  else 
  {
    cout << "    mesh is closed: no\n";
    //	if -p option, print offending
    if (cs.get_print_detailed_info()==1)
    {
      int j=1;
      cout << "    mesh is closed: # border edges - " << border.size() << endl;
      // for each border edge
      for (e_iterator i=border.begin();i!=border.end();i++)
      {
        cout << "    mesh is closed: border edge # " << j++ << endl;
        if (cs.get_use_dreamm_output_format())
        {
          (*i)->printCP(cout);
        }
        else 
        {
          (*i)->print(cout);
          cout << endl;
        }
      }
    }
  }
  // manifold
  if (isManifold()==true)
  {
    cout << "    mesh is manifold: yes\n";
  }
  else if (isManifold()==false && isClosed()==true)
  {
    cout << "    mesh is manifold: no\n";
    //	if -p option, print offending
    if (cs.get_print_detailed_info()==1)
    {
      if (nonman_v.empty())
      {
        cout << "    mesh is manifold: # nonmanifold vertices - none\n";
      }
      else 
      {
        int j=1;
        cout << "    mesh is manifold: # nonmanifold vertices - " << nonman_v.size() << endl;
        // for each nonmanifold vertex
        for (v_iterator i=nonman_v.begin();i!=nonman_v.end();i++)
        {
          cout << "    mesh is manifold: nonmanifold vertex # " << j++ << endl;
          if (cs.get_use_dreamm_output_format())
          {
            (*i)->printCP(cout);
          }
          else 
          {
            (*i)->print(cout);
            cout << endl;
          }
        }
      }
      if (nonman_e.empty())
      {
        cout << "    mesh is manifold: # nonmanifold edges - none\n";
      }
      else 
      {
        int j=1;
        cout << "    mesh is manifold: # nonmanifold edges - " << nonman_e.size() << endl;
        // for each nonmanifold edge
        for (e_iterator i=nonman_e.begin();i!=nonman_e.end();i++)
        {
          cout << "    mesh is manifold: nonmanifold edge # " << j++ << endl;
          if (cs.get_use_dreamm_output_format())
          {
            (*i)->printCP(cout);
          }
          else 
          {
            (*i)->print(cout);
            cout << endl;
          }
        }
      }
    }
  }
  else 
  {
    cout << "    mesh is manifold: undefined since mesh is open\n";
  }
  // consistent
  if (isManifold()==false)
  {
    cout << "    mesh has consistently oriented face normals: undefined since not manifold\n";
  }
  else 
  {
    if (isConsistent()==true)
    {
      cout << "    mesh has consistently oriented face normals: yes\n";
    }
    else 
    {
      cout << "    mesh has consistently oriented face normals: no\n";
      //	if -p option, print offending
      if (cs.get_print_detailed_info()==1)
      {
        int j=1;
        cout << "    mesh has consistently oriented face normals: # flipped edges - " << flipped.size() << endl;
        // for each flipped edge
        // FUTURE IMPROVEMENT: if edge is nonmanifold then exclude from flipped list
        for (e_iterator i=flipped.begin();i!=flipped.end();i++)
        {
          cout << "    mesh has consistently oriented face normals: flipped edge # " << j++ << endl;
          if (cs.get_use_dreamm_output_format())
          {
            (*i)->printCP(cout);
          }
          else 
          {
            (*i)->print(cout);
            cout << endl;
          }
        }
      }
    }
  }
  // outward
  if (isManifold()==false || isConsistent()==false || isClosed()==false)
  {
    cout << "    mesh has outward oriented face normals: uncomputable since ";
    if (isClosed()==false){cout << "not closed,";}
    if (isConsistent()==false){cout << "not consistent,";}
    if (isManifold()==false){cout << "not manifold";}
    cout << endl;
  }
  else 
  {
    if (outward==true)
    {
      cout << "    mesh has outward oriented face normals: yes\n";
    }
    else 
    {
      cout << "    mesh has outward oriented face normals: no\n";
    }
  }
}

void Object::print (Controls &cs)
{
  /*	cout << "\n\n" << "************************ "
        << "OBJECT *************************\n";
  //	print object name 
  cout << "name: " << name << endl;
  */
  //	print Integrity
  printIntegrity(cs);
  if (goodIntegrity()==false)
  {
    cout << "\n\nWarning: Attributes and "
          << "characteristics were not evaluated,\n"
          << " since mesh file failed the integrity check.\n\n";
  }
  else 
  {
    //	print attributes
    printAttr(cs);
    //	print characteristics
    if (cs.get_compute_attributes_only()==0)
    { 
      printChars(cs);
    }
  }
}

void Object::store (Container & c)

{
  Controls & cs(Controls::instance());
  if (cs.get_compute_attributes_only()==0)
  {
    //c.area.insertElementRange(                   area.first_element(),         area.one_past_last_element());
    c.addAreaRange        (         area.first_element(),         area.one_past_last_element());
    c.addAspectRatioRange ( aspect_ratio.first_element(), aspect_ratio.one_past_last_element());
    c.addEdgeLengthRange  (  edge_length.first_element(),  edge_length.one_past_last_element());
    c.addEdgeAngleRange   (   edge_angle.first_element(),   edge_angle.one_past_last_element());
    c.addAdjacentFaceRange(adjacent_face.first_element(),adjacent_face.one_past_last_element());
  }
}

