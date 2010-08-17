#include <string.h>

#include "void_list.h"

void_list* void_list::removeLink (void)
{
  // and remove face from candidate face list
  if (previous!=NULL)
  { 
    if (next!=NULL)
    { 
      // if both previous and next exist
      (previous)->next = next;
      (next)->previous = previous;
    }
    else
    {
      // if previous exists and next does not
      (previous)->next = NULL;
    }
  }
  else
  {
    if (next!=NULL)
    { 
      // if previous does not exist and next does
      (next)->previous = NULL;
    } // else { // if neither previous nor next exists }
  }
  // update pointer
  void_list* q=next;
  delete this;
  return q;
}

void_list* void_list::addPrevious (void)
{
  // go through linked list backwards and add previous pointers
  void_list *prev = NULL;
  for (void_list *p=this;p!=NULL;p=p->next)
  {
    p->previous = prev;
    prev = p;
    if (p->next==NULL) break;
  }
  return this;
}


