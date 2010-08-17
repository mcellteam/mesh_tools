// Author: Justin Kinney
// Date: Sep 2008

#include "void_list.h"
#include <cstring>

void void_list::addPrevious (void)
{
  // go through linked list backwards and add previous pointers
  void_list *p,*prev;
  prev = NULL;
  for (p=this;p!=NULL;p=p->next) {
    p->previous = prev;
    prev = p;
    if (p->next==NULL) break;
  }
}

void_list * void_list::removeLink (void)
{
  //void_list *q;
  // and remove face from candidate face list
  if (previous!=NULL) {
    if (next!=NULL) {
      // if both previous and next exist
      (previous)->next = next;
      (next)->previous = previous;
    } else {
      // if previous exists and next does not
      (previous)->next = NULL;
    }
  } else {
    if (next!=NULL) {
      // if previous does not exist and next does
      (next)->previous = NULL;
    } // else { // if neither previous nor next exists }
  }
  // update pointer
  //q=next;
  //delete this;
  //return q;
  return next;
}

void_list * void_list::addLink (Vertex *v)
{
  void_list *q;
  q = new void_list();
  q->next = this;
  q->data = (void*)v;
  return q;
}


