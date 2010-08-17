#include "face.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>

using std::cout;
using std::endl;

#include "edge.h"

Face::Face (int i, Vertex * v1, Vertex *v2, Vertex *v3)
  :index(i)
{
  v[0] = v1;
  v[1] = v2;
  v[2] = v3;
}

Face::Face (char *triplet,mmap_iv & vp)
  :index(0)
{
  e[0]=e[1]=e[2]=NULL;

  std::pair<iv_iterator,iv_iterator> pp;
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
  int zz = static_cast<int>(strtod(val,&eptr));
  pp = vp.equal_range(zz);
  if (pp.first!=pp.second)
  {
    v[0] = (*(pp.first)).second;
  }
  else
  {
    cout << "Face::Face: ERROR:"
          << "Vertex index not found.\n";
    assert(0);
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
  zz = static_cast<int>(strtod(val,&eptr));
  pp = vp.equal_range(zz);
  if (pp.first!=pp.second)
  {
    v[1] = (*(pp.first)).second;
  }
  else
  {
    cout << "Face::Face: ERROR:"
          << "Vertex index not found.\n";
    assert(0);
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
  zz = static_cast<int>(strtod(val,&eptr));
  pp = vp.equal_range(zz);
  if (pp.first!=pp.second)
  {
    v[2] = (*(pp.first)).second;
  }
  else
  {
    cout << "Face::Face: ERROR:"
          << "Vertex index not found.\n";
    assert(0);
  }
  if (val==eptr)
  {
    v[0]=v[1]=v[2]=NULL;
    printf("Error in reading vertex index\n");
    return;
  }
}

void Face::print (std::ostream & target) const
{
  target.precision(12);
  target << "Face index " << index << endl;
  if (v[0]!=NULL)
  {
    v[0]->print(target);
    target << "\n";
  }
  else
  {
    cout << "Face::print: ERROR:"
          << "Vertex pointer is NULL.\n";
    assert(0);
  }

  if (v[1]!=NULL)
  {
    v[1]->print(target);
    target << "\n";
  }
  else
  {
    cout << "Face::print: ERROR:"
          << "Vertex pointer is NULL.\n";
    assert(0);
  }

  if (v[2]!=NULL)
  {
    v[2]->print(target);
    target << "\n";
  }
  else
  {
    cout << "Face::print: ERROR:"
          << "Vertex pointer is NULL.\n";
    assert(0);
  }

  if (e[0]!=NULL) { target << "[e0 " << e[0] << "]\n";}
  else            { target << "[e0 NULL]\n"; }
  if (e[1]!=NULL) { target << "[e1 " << e[1] << "]\n";}
  else            { target << "[e1 NULL]\n"; }
  if (e[2]!=NULL) { target << "[e2 " << e[2] << "]\n";}
  else            { target << "[e2 NULL]\n"; }
}

void Face::print4mesh (std::ostream & target) const
{
  target.precision(12);
  target << "Face "
        << index << " "
        << v[0]->getIndex() << " "
        << v[1]->getIndex() << " "
        << v[2]->getIndex() << endl;
}

void Face::addEdge (Edge * ptr)
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
