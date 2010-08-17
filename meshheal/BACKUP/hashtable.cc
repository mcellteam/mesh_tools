// Author: Justin Kinney
// Date: Sep 2008

#include "hashtable.h"
#include "void_list.h"

#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>

HashTable::HashTable (const HashTable& rhs)           
  :t(rhs.t),s(rhs.s),n(rhs.n)
{                                                        
}                                             
                                                                               
HashTable& HashTable::operator = (const HashTable& rhs)
{                                                                              
  std::cout << "Assignment operator prohibited on instances of HashTable class.\n";    
  std::cout << "HashTable " << rhs.s << std::endl;                                      
  exit(1);                                               
}                                                        
   
HashTable::HashTable (int a)
:t(NULL),s(0),n(0)
{
  n=4;
  // a = number of edges in block
  int bins = a/n+1;
  int i=0;
  //while(bins=bins/2){i++;}
  bins=bins/2;
  while(bins)
  {
    i++;
    bins=bins/2;
  }
  s=1;
  while(i){s=s*2;i--;}
  t=new void_list*[s];
  for(i=0;i<s;i++){
    t[i]=NULL;
  }
}

HashTable::~HashTable(void){
  void_list *p,*q;
  int i;
  for(i=0;i<s;i++){
    p=t[i];
    while (p!=NULL) {
      q=p->next;
      delete p;
      p=q;
    }
  }
  delete[] t;
}



