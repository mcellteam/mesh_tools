// Author: Justin Kinney
// Date: Sep 2008

#ifndef VOIDLIST_H
#define VOIDLIST_H 1

class Vertex;

class void_list
{
public:
  void_list *previous;
  void_list *next;
  void *data;
public:
  void addPrevious (void);
  void_list * addLink (Vertex *v);
  void_list * removeLink (void);
};

#endif

