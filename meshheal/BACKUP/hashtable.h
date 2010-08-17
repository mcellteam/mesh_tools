// Author: Justin Kinney
// Date: Sep 2008

#ifndef HASHTABLE_H
#define HASHTABLE_H 1

#include "meshheal.h"

class HashTable
{
public:
  void_list **t;  // pointer to void_list * // use for array of void_list *
  int s;          // size of array
  int n;          // edges per hash bin
  HashTable(int);
  ~HashTable(void);
  HashTable              (HashTable const &);
  HashTable & operator = (HashTable const &);
};

#endif    
