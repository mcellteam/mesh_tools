#ifndef SITE_H
#define SITE_H 1

#include "vertex.h"
#include "extra_vertex.h"

class Site
{
private:
  Site & operator = (Site const &);
public:
  Site (Site const &);
  int n; // number of extra vertices in site
  int orient; // orientation of l1,l2 and th
  Vertex *l1,*l2,*th;
  ExtraVertex **ev;
  Site(void);
  ~Site(void);
  void incrementN (void)
  {
    n++;
  }
  int getVertCount (void)
  {
    return n;
  }
  void initExtraVertex (void)
  {
    ev = new ExtraVertex*[n];
    // initialize ExtraVertex*
    for (int k=0;k<n;k++)
    {
      ev[k]=NULL;
    }
  }
};

#endif
