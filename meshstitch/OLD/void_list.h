#ifndef VOID_LIST_H
#define VOID_LIST_H 1

class void_list
{
public:
  void_list *previous;
  void_list *next;
  void *data;
  void_list* removeLink(void); 
  void_list* addPrevious(void);
};

#endif
