// Author: Justin Kinney
// Date: Sep 2008

#include "edgeblock.h"
#include "void_list.h"

#include <iostream>
#include <cassert>
#include <cstring>
#include <cstdio>
#include <cstdlib>

EdgeBlock::EdgeBlock (const EdgeBlock& rhs)           
  :e(rhs.e),n(rhs.n),p(rhs.p),ht(rhs.ht)
{                                                        
}                                             
                                                                               
EdgeBlock& EdgeBlock::operator = (const EdgeBlock& rhs)
{                                                                              
  std::cout << "Assignment operator prohibited on instances of EdgeBlock class.\n";    
  std::cout << "EdgeBlock " << rhs.n << std::endl;                                      
  exit(1);                                               
}                                                        
   
EdgeBlock::EdgeBlock (int m)
  :e(NULL),n(m),p(0),ht(NULL)
{
  e = new Edge[m];
  ht = new HashTable(n);
}

EdgeBlock::~EdgeBlock(void)
{
  delete[] e;
  delete ht;
}

void EdgeBlock::getEdges (void_list * flh)
{
  // for each face
  for (void_list * q=flh;q!=NULL;q=q->next)
  {
    Face * f=(Face*)q->data;
    checkEdge(f,f->v1,f->v2);
    checkEdge(f,f->v2,f->v3);
    checkEdge(f,f->v3,f->v1);
  }
}

void EdgeBlock::checkEdge (Face *F,
                           int va,
                           int vb)
{
  Edge *ee = findEdge(va,vb);
  if(ee){ updateEdge(F,ee,va,vb);}
  else {createEdge(F,va,vb);}
}

Edge* EdgeBlock::findEdge (const int & va,
                           const int & vb)
{
  // compute hashval of va,vb
  int hashval = getEdgeHashVal(va,vb,ht->s);
  if (hashval>ht->s-1)
  {
    printf("hashval %i, size of hashtable %i\n",hashval,ht->s);
  }
  // for each edge in hashtable element pointed to by hashval
  void_list * pp = ht->t[hashval];
  bool found = false;
  Edge * ee = NULL;
  while (pp !=NULL && !found)
  {
    ee=(Edge*)pp->data;
    if  ((ee->v1==va||ee->v2==va)&&(ee->v1==vb||ee->v2==vb)){found=true;}
    else {pp=pp->next;}
  }
  if (found) return ee;
  return NULL;
}

void EdgeBlock::updateEdge (Face * F,
                            Edge * ee,
                            const int & va,
                            const int & vb)
{
  assert(F!=NULL);
  assert(ee!=NULL);
  //add face to edge
  void_list * qq = new void_list();
  qq->next = ee->f;
  qq->data = (void*)F;
  ee->f = qq;
  if  (ee->v1==va && ee->v2==vb) {ee->c12+=1;}
  else {ee->c21+=1;}
}

void EdgeBlock::createEdge (Face *F,
                            int va,
                            int vb)
{
  Edge *ee=NULL;
  bool found=false;
  // is there an open Edge slot
  if (p<(n-1))
  {
    // if yes grab Edge pointer
    ee=&e[p];
    p++;
    found = true;
  }
  if (!found)
  {
    printf("Increase number of edges in EdgeBlock!\n");
    exit(1);
  }
  // load Edge
  ee->v1 = va;
  ee->v2 = vb;
  ee->c12 = 1;
  ee->c21 = 0;
  // add face to edge
  // set link data to point to Face class instance
  ee->f = new void_list();
  (ee->f)->next = NULL;
  (ee->f)->data = (void*)F;
  // store edge pointer in hash table
  // compute hashval of va,vb
  int hashval = getEdgeHashVal(va,vb,ht->s);
  // load hash table
  void_list * pp=new void_list();
  pp->next = ht->t[hashval];
  pp->data= (void*)ee;
  ht->t[hashval] = pp;
}

int EdgeBlock::getEdgeHashVal (int va,int vb,int mask)
{
  u4 hashval1 = computeHashValue(va);
  u4 hashval2 = computeHashValue(vb);
  u4 result = hashval1+hashval2;
  return (int)(result&(mask-1));
}


/* The whole new hash function */
u4 EdgeBlock::hash (register u1* k, u4 length, u4 initval)
  //register u1 *k;        /* the key */
  //u4           length;   /* the length of the key in bytes */
  //u4           initval;  /* the previous hash, or an arbitrary value */
{
  register u4 a,b,c;  /* the internal state */
  u4          len;    /* how many key bytes still need mixing */

  /* Set up the internal state */
  len = length;
  a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
  c = initval;         /* variable initialization of internal state */

  /*---------------------------------------- handle most of the key */
  while (len >= 12)
  {
    a=a+(k[0]+((u4)k[1]<<8)+((u4)k[2]<<16) +((u4)k[3]<<24));
    b=b+(k[4]+((u4)k[5]<<8)+((u4)k[6]<<16) +((u4)k[7]<<24));
    c=c+(k[8]+((u4)k[9]<<8)+((u4)k[10]<<16)+((u4)k[11]<<24));
    mix(a,b,c);
    k = k+12; len = len-12;
  }

  /*------------------------------------- handle the last 11 bytes */
  c = c+length;
  switch(len)              /* all the case statements fall through */
  {
    case 11: c=c+((u4)k[10]<<24);
    case 10: c=c+((u4)k[9]<<16);
    case 9 : c=c+((u4)k[8]<<8);
             /* the first byte of c is reserved for the length */
    case 8 : b=b+((u4)k[7]<<24);
    case 7 : b=b+((u4)k[6]<<16);
    case 6 : b=b+((u4)k[5]<<8);
    case 5 : b=b+k[4];
    case 4 : a=a+((u4)k[3]<<24);
    case 3 : a=a+((u4)k[2]<<16);
    case 2 : a=a+((u4)k[1]<<8);
    case 1 : a=a+k[0];
             /* case 0: nothing left to add */
  }
  mix(a,b,c);
  /*-------------------------------------------- report the result */
  return c;
}

u4 EdgeBlock::computeHashValue (int key)
{
  char cval[10];
  u4 result;
  sprintf(cval,"%i",key);
  result = hash((u1*)cval,(u4)strlen(cval),0xa205b064);
  return result;
}

double EdgeBlock::computeLongestEdge (Vertex ** vert_array,
                                      double epsilon)
{
  double max=0.0;
  // for each edge in block
  for (int i=0;i<n;i++)
  {
    Edge * ee=&e[i];
    // if edge is valid
    if (ee->isValid())
    {
      // compute squared edge length
      double d = ee->getEdgeLengthSq(vert_array,epsilon);
      // if edge is longer than max, then save d
      if (d>max) max=d;
    }
  }
  return max; // return threshold squared
}

void_list * EdgeBlock::gatherFreeVertices(int max_vert,
                                          Vertex **vert_array)
{
  void_list * vfree = NULL;
  // for each edge in block
  for (int i=0;i<n;i++)
  {
    Edge * ee = &e[i];
    // if edge is free
    if ((ee->c12+ee->c21)==1)
    {
      void_list * pp = new void_list();
      pp->next = vfree;
      pp->data = (void*)vert_array[ee->v1];
      vfree = pp;
      pp = new void_list();
      pp->next = vfree;
      pp->data = (void*)vert_array[ee->v2];
      vfree = pp;
    }
  }

  // go through linked list backwards and add previous pointers
  vfree->addPrevious();

  ///// gather unique free vertices /////
  int v_array[max_vert+1];
  for (int i=0;i<max_vert+1;i++) {v_array[i]=0;}
  void_list * ptr=vfree;
  void_list * pp=vfree;
  while (pp!=NULL)
  {
    Vertex * v =(Vertex*)pp->data;
    if(!v_array[v->index]){v_array[v->index]=1;pp=pp->next;}
    else {
      // adjust list pointer
      if (pp==ptr) { ptr=pp->next;}
      // remove free face link
      void_list * qq = pp;
      pp = pp->removeLink();
      delete qq;
    }
  }
  // adjust pointer
  return ptr;
}

