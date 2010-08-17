// Author: Justin Kinney
// Date: Sep 2008

#include "meshheal.h"
                      
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>
#include <iostream>
#include <cassert>

using std::cout;
using std::endl;

#include "distance.h"
#include "edgeblock.h"
#include "face.h"
#include "vertex.h"
#include "void_list.h"

void getData (char * infile,
              void_list * & flh,
              void_list * & vlh)
{

  // open first file
  FILE * F = fopen(infile,"r");
  if (!F) { printf("Couldn't open input file %s\n",infile);return;}

  // for every line in first file
  char line[2048];
  for (char *str=fgets(line,2048,F) ; str!=NULL ; str=fgets(line,2048,F))
  {
    // skip leading whitespace
    while (strchr(" \t,",*str)!=NULL) { str++;}
    // if first character is V for Vertex, add new linked list class instance
    if (strchr("V",*str)!=NULL)
    {
      void_list * vl = new void_list();
      vl->next = vlh;
      Vertex * v = new Vertex(str);
      vl->data = (void*)v;
      vlh = vl;
    } 
    // if first character is F for Face, add new linked list class instance
    else if (strchr("F",*str)!=NULL)
    {
      void_list * fl = new void_list();
      fl->next = flh;
      Face * f = new Face(str);
      fl->data = (void*)f;
      flh = fl;
    }
  }
  fclose(F);
}

int maxVert (void_list * L)
{
  int max=0;
  for (void_list *p=L;p!=NULL;p=p->next)
  {
    if (max<((Vertex*)p->data)->index) {max=((Vertex*)p->data)->index;}
  }
  return max;
}

int maxFace (void_list * L)
{
  int max=0;
  for (void_list *p=L;p!=NULL;p=p->next)
  {
    if (max<((Face*)p->data)->index) {max=((Face*)p->data)->index;}
  }
  return max;
}

Vertex ** initVertArray (void_list * vlh,
                         int max_vert)
{
  Vertex ** vert_array = new Vertex*[max_vert+1];
  for (void_list *p=vlh;p!=NULL;p=p->next)
  {
    vert_array[((Vertex*)p->data)->index]=(Vertex*)p->data;
  }
  return vert_array;
}

int compare (const void* a, const void* b )
{
  assert(a!=NULL);
  assert(b!=NULL);

  void_list *j,*k;

  j = *(void_list**)a;
  k = *(void_list**)b;

  assert(j!=NULL);
  assert(k!=NULL);
  assert((Distance*)j->data!=NULL);
  assert((Distance*)k->data!=NULL);

  double i = ((Distance*)j->data)->d - ((Distance*)k->data)->d;
  if (i<0) {	return -1;}
  else if (i>0) {	return 1;}
  else { return (0);}
}

void_list ** computeDistances (EdgeBlock * eb,
                               void_list * v,
                               int & num_free_vertex_pairs,
                               double epsilon,
                               int print_flag,
                               const double & threshold)
{
  fprintf(stderr,"Computing vertex distances...");fflush(stderr);
  void_list *dlh = NULL;
  num_free_vertex_pairs = 0;
  // for each free vertex
  for (void_list * q=v;q!=NULL;q=q->next)
  {
    Vertex *vo=(Vertex*)q->data;
    // for each free vertex
    for (void_list * qq=v;qq!=NULL;qq=qq->next)
    {
      Vertex *vi=(Vertex*)qq->data;
      // if combination is unique
      if (vo->index > vi->index)
      {
        // are vertices members of the same face?
        Edge * e=eb->findEdge(vo->index,vi->index);
        if (!e)
        {
          // do not factor out epsilon, leave as is
          double diffx = vo->x/epsilon-vi->x/epsilon;
          double diffy = vo->y/epsilon-vi->y/epsilon;
          double diffz = vo->z/epsilon-vi->z/epsilon;
          double squared_dist = diffx*diffx+diffy*diffy+diffz*diffz;
          // if distance less than longest edge, then save
          if (squared_dist<threshold)
          {
            // create distance
            void_list * dl = new void_list();
            assert(dl!=NULL);
            dl->next = dlh;
            Distance * d = new Distance(squared_dist,q,qq);
            assert(d!=NULL);
            dl->data = (void*)d;
            assert(dl->data!=NULL);
            dlh = dl;
            assert(dlh!=NULL);
            num_free_vertex_pairs++;
          }
        }
      }
    }
  }

  ///// sort distances /////
  // create array of distances
  void_list ** dist_array = new void_list*[num_free_vertex_pairs];
  int i=0;
  for (void_list * q=dlh;q!=NULL;q=q->next)
  {
    assert(i<num_free_vertex_pairs);
    assert(q!=NULL);
    assert((Distance*)q->data!=NULL);
    dist_array[i]=q;
    i++;
  }
  assert(i==num_free_vertex_pairs);
  qsort(dist_array,num_free_vertex_pairs,sizeof(void_list*),compare);
  if (!print_flag) {printf("complete.\n");fflush(stdout);}
  return dist_array;
}

int countFreeVertices (void_list *v)
{
  ///// count number free vertices /////
  int i=0;
  for (void_list *q=v;q!=NULL;q=q->next) i++;
  fprintf(stderr,"num free vertices = %i\n",i);fflush(stderr);
  return i;
}

int findVerticesToMerge (void_list ** array,
                         int count,
                         void_list * vfree,
                         double epsilon)
{
  int i=0;
  Distance * d = NULL;
  bool flag = false;
  while (i<count && !flag)
  {
    d=(Distance*)array[i]->data;
    // is ith vertex pair free?
    bool vflag1=false;
    bool vflag2=false;
    void_list * q=vfree;
    // for each free vertex
    while (q!=NULL && (!vflag1||!vflag2))
    {
      Vertex * v=(Vertex*)q->data;
      if (((Vertex*)d->vA)->index==v->index){vflag1=true;}
      if (((Vertex*)d->vB)->index==v->index){vflag2=true;}
      q=q->next;
    }
    if (vflag1&&vflag2) flag=true;
    else {i++;}
  }
  if (i==count)
  {
    fprintf(stderr,"ERROR! No eligible vertices to merge were found.\n");
  }
  else
  {
    fprintf(stderr,"vertices to merge: %i %i, dist = %.15g\n",
            ((Vertex*)d->vA)->index,((Vertex*)d->vB)->index,sqrt(d->d*epsilon*epsilon));
  }
  return i;
}

void_list * fixVertices (void_list * v,
                         Distance * d,
                         void_list * & vbad)
{
  ///// remove second vertex in pair from vertex list /////
  void_list *ptr=v;
  void_list *p=v;
  bool found = false;
  while (p!=NULL && !found)
  {
    if (((Vertex*)d->vB)->index == ((Vertex*)p->data)->index)
    {
      // add pointer to bad vertex list
      vbad = vbad->addLink((Vertex*)p->data);
      // adjust list pointer
      if (p==ptr) { ptr=p->next;}
      // remove link from list
      void_list * q = p;
      p = p->removeLink();
      delete q;
      found = true;
    }
    else {p=p->next;}
  }
  // adjust pointer
  return ptr;
}

void_list * fixFaces (void_list * flh,
                      Distance * d)
{
  ///// replace all instances of second vertex in face list with first vertex /////
  Vertex *va=(Vertex*)d->vA;
  Vertex *vb=(Vertex*)d->vB;
  for (void_list * p=flh;p!=NULL;p=p->next)
  {
    Face *f=(Face*)p->data;
    if(vb->index==f->v1){f->v1=va->index;}
    if(vb->index==f->v2){f->v2=va->index;}
    if(vb->index==f->v3){f->v3=va->index;}
  }
  // tag distance class for future vertex b deletion
  d->deleteme=true;
  return flh;
}

void deleteFreeVertices (void_list *v)
{
  void_list *p=v;
  while (p!=NULL)
  {
    void_list *q=p->next;
    delete p;
    p=q;
  }
}

void printVerticesFaces (void_list * vlh,
                         void_list * flh)
{
  void_list *vend=NULL,*fend=NULL;
  for (void_list *q=vlh;q!=NULL;q=q->next) {vend=q;}
  for (void_list *q=flh;q!=NULL;q=q->next) {fend=q;}
  // write out vertex linked list
  for (void_list *q=vend;q!=NULL;q=q->previous)
  {
    Vertex *v=(Vertex*)q->data;		
    printf("Vertex %i  %.15g %.15g %.15g\n",v->index,v->x,v->y,v->z);
  }
  // write out face linked list
  for (void_list *q=fend;q!=NULL;q=q->previous)
  {
    Face *f=(Face*)q->data;
    printf("Face %i  %i %i %i\n",f->index,f->v1,f->v2,f->v3);
  }
}

void deleteEdges (EdgeBlock *eb)
{
  delete eb;
}

void cleanup (EdgeBlock * eb,
              void_list * vf,
              void_list * v,
              void_list * f,
              void_list ** d,
              int count,
              Vertex ** vert_array,
              void_list * vbad)
{
  void_list *p,*q;
  deleteEdges(eb);
  deleteFreeVertices(vf);
  ////////// delete vertices and faces /////////
  p=v;
  while (p!=NULL) {
    q=p->next;
    delete (Vertex*)p->data;
    delete p;
    p=q;
  }
  p=f;
  while (p!=NULL) {
    q=p->next;
    delete (Face*)p->data;
    delete p;
    p=q;
  }
  ////////// delete dist_array /////////
  for(int i=0;i<count;i++){
    delete (Distance*)d[i]->data;
    delete d[i];
  }
  delete[] d;
  ///////////////////////
  delete[] vert_array;
  ////////// delete bad vertices /////////
  p=vbad;
  while (p!=NULL) {
    q=p->next;
    delete (Vertex*)p->data;
    delete p;
    p=q;
  }
}

int main (int argc,char *argv[])
{

  if (argc != 4)
  {
    printf("\nSyntax: meshheal input_file epsilon print_flag\n\n");
    printf("Description: Fills holes in meshes caused by vertices separated\n");
    printf("             by small a  distance, epsilon, and therefore not being treated\n");
    printf("             as the same vertex. This script takes a user-defined epsilon\n");
    printf("             and identifies and removes duplicate vertices,\n");
    printf("             updates polygon vertex info, and writes\n");
    printf("             output to stdout. print_flag = 1 to print\n");
    printf("             updated mesh to stdout. print_flag = 0 to not print.\n\n");
    return 1;
  }

  /*
      This program works as follows. For each vertex on a border, i.e. open edge, 
      compute the squared distances to neighboring vertices in 3D to identify the closest
      vertex that is also on a border. Declare these two vertices to be duplicates.
   */

  ////////// get data /////////
  char * eptr;
  const double epsilon = strtod(argv[2],&eptr);
  const int print_flag = (int)strtod(argv[3],&eptr);
  void_list * flh = NULL , * vlh = NULL;
  getData(argv[1],flh,vlh);
  flh->addPrevious();
  vlh->addPrevious();
  int num_verts = maxVert(vlh);
  int num_faces = maxFace(flh);

  // create array of pointers to vertex list
  Vertex ** vert_array = initVertArray(vlh,num_verts);

  ///// build edge list ////
  int num_edges = num_faces+num_verts-2+10000;
  EdgeBlock * eb = new EdgeBlock(num_edges);
  eb->getEdges(flh);

  // get max edge length squared
  const double threshold = eb->computeLongestEdge(vert_array,epsilon);

  ///// gather free vertices /////
  void_list * vfree = eb->gatherFreeVertices(num_verts,vert_array);

  ///// compute distances /////
  int num_free_vertex_pairs = 0;
  void_list ** dist_array = computeDistances(eb,vfree,num_free_vertex_pairs,epsilon,print_flag,threshold);

  ///// merge vertices /////
  void_list * vbad=NULL;
  while(1)
  {
    ///// identify free vertex pair separated by the smallest distance /////
    int closest_free_vertex_pair = findVerticesToMerge(dist_array,num_free_vertex_pairs,vfree,epsilon);
    if (closest_free_vertex_pair==num_free_vertex_pairs) break;

    ///// remove second vertex in pair from vertex list /////
    vlh = fixVertices(vlh,(Distance*)dist_array[closest_free_vertex_pair]->data,vbad);

    ///// replace all instances of second vertex in face list with first vertex /////
    flh = fixFaces(flh,(Distance*)dist_array[closest_free_vertex_pair]->data);

    ///// delete edge list /////
    deleteEdges(eb);

    ///// build edge list ////
    num_edges=num_faces+num_verts-2+10000;
    eb = new EdgeBlock(num_edges);
    eb->getEdges(flh);

    ///// delete free vertices /////
    deleteFreeVertices(vfree);

    ///// gather free vertices /////
    vfree = eb->gatherFreeVertices(num_verts,vert_array);

    ///// count number free vertices /////
    if (countFreeVertices(vfree)==0) break;
  }
  // print to stdout
  if(print_flag){printVerticesFaces(vlh,flh);}

  ///// cleanup /////
  cleanup(eb,vfree,vlh,flh,dist_array,num_free_vertex_pairs,vert_array,vbad);

  return 0;
}
