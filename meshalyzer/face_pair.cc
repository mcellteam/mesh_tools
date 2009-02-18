#include "face_pair.h"

#include <iostream>

#include "edge.h"
#include "object.h"

using std::cout;
using std::endl;

bool Face_Pair::existingEdge (void)
{
  // for each edge in object
  for (e_iterator i=o->e.begin();i!=o->e.end();i++)
  {
    // if edge vertices are b and d
    if ( ((*i)->vv1==b && (*i)->vv2==d) || 
        ((*i)->vv1==d && (*i)->vv2==b)){return true;}
  }
  return false;
}

bool Face_Pair::aspectRatiosImprove (void)
{
  double arf1 = f1->getAspectRatio();
  double arf2 = f2->getAspectRatio();
  Face nf1(a,b,d);
  Face nf2(b,c,d);
  double arnf1 = nf1.getAspectRatio();
  double arnf2 = nf2.getAspectRatio();

  double old_max,new_max;
  if (arf1>arf2){old_max=arf1;}
  else {old_max=arf2;}
  if (arnf1>arnf2){new_max=arnf1;}
  else {new_max=arnf2;}

  if (new_max<old_max) {return true;}
  else {return false;}
}

void Face_Pair::print(void)
{
  cout << "Follow these instructions to improve face aspect ratio.\n"
        << "open " << o->name << endl
        << "remove Face "
        << f1->index << " "
        << f1->v[0]->index << " "
        << f1->v[1]->index << " "
        << f1->v[2]->index << endl
        << "remove Face "
        << f2->index << " "
        << f2->v[0]->index << " "
        << f2->v[1]->index << " "
        << f2->v[2]->index << endl
        << "add Face "
        << next_i++ << " "
        << a->index << " "
        << b->index << " "
        << d->index << endl
        << "add Face "
        << next_i++ << " "
        << b->index << " "
        << c->index << " "
        << d->index << endl;
}

void Face_Pair::findF2 (void)
{
  // for each face in object
  for (f_iterator i=o->f.begin();i!=o->f.end();i++)
  {
    // if face not f1
    if (*i!=f1)
    {
      // if face contains edge e
      if ((*i)->e[0]==e||(*i)->e[1]==e||(*i)->e[2]==e)
      {
        f2=*i;
        break;
      }
    }
  }
  //check
  if (f2==NULL)
  {
    cout << "Face_Pair::findF2: Unable to identify f2.\n";
    exit(1);
  }
  // identify d
  if (f2->v[0]!=a && f2->v[0]!=c){d=f2->v[0];}
  else if (f2->v[1]!=a && f2->v[1]!=c){d=f2->v[1];}
  else if (f2->v[2]!=a && f2->v[2]!=c){d=f2->v[2];}
  else {cout << "Face_Pair::findF2: Error: Unable to identify d.\n";
    exit(1);
  }
}

void Face_Pair::analyzeF1 (void)
{
  double ll = -1.0;
  // find longest edge
  for (int i=0;i<3;i++)
  {
    if (f1->e[i]->l>ll)
    {
      e=f1->e[i];
      ll=f1->e[i]->l;
    }
  }
  // identify a,b,c
  if	  (f1->v[0]==e->vv1 && f1->v[1]==e->vv2)
  {
    a=e->vv2;
    c=e->vv1;
    b=f1->v[2];
  }
  else if (f1->v[1]==e->vv1 && f1->v[2]==e->vv2)
  {
    a=e->vv2;
    c=e->vv1;
    b=f1->v[0];
  }
  else if (f1->v[2]==e->vv1 && f1->v[0]==e->vv2)
  {
    a=e->vv2;
    c=e->vv1;
    b=f1->v[1];
  }
  else if  (f1->v[0]==e->vv2 && f1->v[1]==e->vv1)
  {
    a=e->vv1;
    c=e->vv2;
    b=f1->v[2];
  }
  else if (f1->v[1]==e->vv2 && f1->v[2]==e->vv1)
  {
    a=e->vv1;
    c=e->vv2;
    b=f1->v[0];
  }
  else if (f1->v[2]==e->vv2 && f1->v[0]==e->vv1)
  {
    a=e->vv1;
    c=e->vv2;
    b=f1->v[1];
  }
  else
  {
    cout << "Face_Pair::analyzeF1: Error: Unable to identify a,b,c.\n";
    exit(1);
  }
}

Face_Pair::Face_Pair (Object *oo)
  :a(NULL),b(NULL),c(NULL),d(NULL),e(NULL),f1(NULL),
   f2(NULL),o(NULL),next_i(0)
{
  int max=-1;
  o=oo;
  // for each face in object
  for (f_iterator i=o->f.begin();i!=o->f.end();i++)
  {
    if ((*i)->index > max){max=(*i)->index;}
  }
  next_i=max+1;
  //
  a=b=c=d=NULL;
  f1=f2=NULL;
  e=NULL;
}

void Face_Pair::clear (void)
{
  a=b=c=d=NULL;
  f1=f2=NULL;
  e=NULL;
}

void Face_Pair::processBadFace (Face *f)
{
  clear();
  f1=f;
  analyzeF1();
  findF2();
  // if reconfiguring the faces would improve the face aspect ratios
  // i.e. the largest aspect ratio would get smaller
  if (aspectRatiosImprove()==true && existingEdge()==false)
  {
    print();
  }
}

