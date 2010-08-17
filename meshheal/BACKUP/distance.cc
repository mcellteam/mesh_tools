// Author: Justin Kinney
// Date: Sep 2008

#include "distance.h"
#include "void_list.h"

#include <cassert>
#include <cstring>
#include <cstdio>
#include <cstdlib>

Distance::Distance (double dist,
                    void_list * p,
                    void_list * q)
  :d(dist),vA(NULL),vB(NULL),deleteme(false)
{
  assert(p!=NULL);
  assert(q!=NULL);
  assert((Vertex*)p->data!=NULL);
  assert((Vertex*)q->data!=NULL);
  vA = (Vertex*)p->data;
  vB = (Vertex*)q->data;
}
