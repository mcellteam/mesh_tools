#include "contour.h"
#include "meshstitch.h"

Contour::Contour(void)
:num(0),extra(),verts(),faces(),s(NULL)
{
};

Contour::~Contour(void)
{
  //delete[] s;
};

void Contour::initSites (void)
{
  s.reserve(num);
  // for each site in contour
  for (int j=0;j<num;j++)
  {
    // for each extra vertex
    for (iter_vec_ev p=extra.begin();p!=extra.end();p++)
    {
      ExtraVertex *ev=&(*p);
      // if extra vertex group number matches site group
      if (ev->g==j)
      {
        //c[i].s[j].n++;
        s[j].incrementN();
      }
    }
    // allocate space for ExtraVertex*
    //c[i].s[j].ev = new ExtraVertex*[c[i].getSiteVertCount(j)];
    s[j].initExtraVertex();
  }
}
