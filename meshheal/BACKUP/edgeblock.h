// Author: Justin Kinney
// Date: Sep 2008

#ifndef EDGEBLOCK_H
#define EDGEBLOCK_H 1

#include "edge.h"
#include "hashtable.h"
#include "vertex.h"

typedef  unsigned long int  u4;   /* unsigned 4-byte type */
typedef  unsigned     char  u1;   /* unsigned 1-byte type */

/* The mixing step */
#define mix(a,b,c) \
{ \
  a=a-b;  a=a-c;  a=a^(c>>13); \
  b=b-c;  b=b-a;  b=b^(a<<8);  \
  c=c-a;  c=c-b;  c=c^(b>>13); \
  a=a-b;  a=a-c;  a=a^(c>>12); \
  b=b-c;  b=b-a;  b=b^(a<<16); \
  c=c-a;  c=c-b;  c=c^(b>>5);  \
  a=a-b;  a=a-c;  a=a^(c>>3);  \
  b=b-c;  b=b-a;  b=b^(a<<10); \
  c=c-a;  c=c-b;  c=c^(b>>15); \
}

class EdgeBlock
{
public:
  Edge *e;            // pointer to edge
  int n;              // number of edges in block
  int p;              // next open edge in block
  HashTable * ht;
public:
  EdgeBlock(int m);
  ~EdgeBlock(void);
  EdgeBlock              (EdgeBlock const &);
  EdgeBlock & operator = (EdgeBlock const &);
  void getEdges(void_list *flh);
  void checkEdge(Face *F,int va,int vb);
  Edge* findEdge(const int & va,const int & vb);
  void updateEdge(Face *F,Edge *e,const int & va,const int & vb);
  void createEdge(Face *F,int va,int vb);
  int getEdgeHashVal(int va,int vb,int mask);
  u4 computeHashValue(int key);
  u4 hash( register u1* k, u4 length, u4 initval);
  double computeLongestEdge (Vertex **vert_array, double epsilon);
  void_list * gatherFreeVertices(int max_vert,
                                 Vertex **vert_array);
};

#endif

