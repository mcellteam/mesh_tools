#include "face_pair.h"

#include <iostream>
#include <stdlib.h>

#include "edge.h"
#include "object.h"

using std::cout;
using std::endl;

bool Face_Pair::existingEdge (void)
{
  // for each edge in object
  for (e_iterator i=o->getFirstEdge();i!=o->getOnePastLastEdge();i++)
  {
    // if edge vertices are b and d
    if ( ((*i)->ptr_vv1()==b && (*i)->ptr_vv2()==d) || 
        ((*i)->ptr_vv1()==d && (*i)->ptr_vv2()==b)){return true;}
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
        << "open " << o->getName() << endl
        << "remove Face "
        << f1->getIndex() << " "
        << f1->ptr_vertex(0)->getIndex() << " "
        << f1->ptr_vertex(1)->getIndex() << " "
        << f1->ptr_vertex(2)->getIndex() << endl
        << "remove Face "
        << f2->getIndex() << " "
        << f2->ptr_vertex(0)->getIndex() << " "
        << f2->ptr_vertex(1)->getIndex() << " "
        << f2->ptr_vertex(2)->getIndex() << endl
        << "add Face "
        << next_i++ << " "
        << a->getIndex() << " "
        << b->getIndex() << " "
        << d->getIndex() << endl
        << "add Face "
        << next_i++ << " "
        << b->getIndex() << " "
        << c->getIndex() << " "
        << d->getIndex() << endl;
}

void Face_Pair::findF2 (void)
{
  // for each face in object
  for (f_iterator i=o->getFirstFace();i!=o->getOnePastLastFace();i++)
  {
    // if face not f1
    if (*i!=f1)
    {
      // if face contains edge e
      if ((*i)->ptr_edge(0)==e||(*i)->ptr_edge(1)==e||(*i)->ptr_edge(2)==e)
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
  if      (f2->ptr_vertex(0)!=a && f2->ptr_vertex(0)!=c){d=f2->ptr_vertex(0);}
  else if (f2->ptr_vertex(1)!=a && f2->ptr_vertex(1)!=c){d=f2->ptr_vertex(1);}
  else if (f2->ptr_vertex(2)!=a && f2->ptr_vertex(2)!=c){d=f2->ptr_vertex(2);}
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
    if (f1->ptr_edge(i)->getOrigLength()>ll)
    {
      e=f1->ptr_edge(i);
      ll=f1->ptr_edge(i)->getOrigLength();
    }
  }
  // identify a,b,c
  if	  (f1->ptr_vertex(0)==e->ptr_vv1() && f1->ptr_vertex(1)==e->ptr_vv2())
  {
    a=e->ptr_vv2();
    c=e->ptr_vv1();
    b=f1->ptr_vertex(2);
  }
  else if (f1->ptr_vertex(1)==e->ptr_vv1() && f1->ptr_vertex(2)==e->ptr_vv2())
  {
    a=e->ptr_vv2();
    c=e->ptr_vv1();
    b=f1->ptr_vertex(0);
  }
  else if (f1->ptr_vertex(2)==e->ptr_vv1() && f1->ptr_vertex(0)==e->ptr_vv2())
  {
    a=e->ptr_vv2();
    c=e->ptr_vv1();
    b=f1->ptr_vertex(1);
  }
  else if  (f1->ptr_vertex(0)==e->ptr_vv2() && f1->ptr_vertex(1)==e->ptr_vv1())
  {
    a=e->ptr_vv1();
    c=e->ptr_vv2();
    b=f1->ptr_vertex(2);
  }
  else if (f1->ptr_vertex(1)==e->ptr_vv2() && f1->ptr_vertex(2)==e->ptr_vv1())
  {
    a=e->ptr_vv1();
    c=e->ptr_vv2();
    b=f1->ptr_vertex(0);
  }
  else if (f1->ptr_vertex(2)==e->ptr_vv2() && f1->ptr_vertex(0)==e->ptr_vv1())
  {
    a=e->ptr_vv1();
    c=e->ptr_vv2();
    b=f1->ptr_vertex(1);
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
  for (f_iterator i=o->getFirstFace();i!=o->getOnePastLastFace();i++)
  {
    if ((*i)->getIndex() > max){max=(*i)->getIndex();}
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

