#ifndef CONTOUR_H
#define CONTOUR_H 1

#include <list>
#include <vector>

#include "face.h"
#include "site.h"

class Contour
{
private:
  Contour (Contour const &);
  Contour & operator = (Contour const &);
  int num;	                  // number of sites
public:
  std::vector<ExtraVertex> extra;
  std::list<Vertex>        verts; // vertices on contour
  std::list<Face>          faces; // faces with two vertices from verts list
  //Site *s;
  Contour(void);
  ~Contour(void);
  void setNum (int i)
  {
    num = i;
  }
  int getNum (void)
  {
    return num;
  }
  void initSites (void);
  void incrementSiteVertCount (int i)
  {
    s[i].incrementN();
  }
  int getSiteVertCount (int i)
  {
    return s[i].getVertCount();
  }
  void initSiteExtraVertex (int i)
  {
    s[i].initExtraVertex();
  }
private:
  std::vector<Site> s;
};

#endif
