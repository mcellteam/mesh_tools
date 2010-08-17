#include "contour.h"

#include "control_points.h"

Control_Points::Control_Points ()
:index()
{
}

/** Describe the control points for each spline as collections
  of raw contour points.
 * \param[in] closed_path Sequence of raw contour points to use as control points.
 */

void Control_Points::loadMatrix (std::list<int> & closed_path)
{
  std::list<int>::iterator p=closed_path.begin();
  std::list<int>::iterator q=p;q++;
  std::list<int>::iterator r=q;r++;
  std::list<int>::iterator s=r;s++;
  while (s!=closed_path.end())
  {
    index.push_back(*p);
    index.push_back(*q);
    index.push_back(*r);
    index.push_back(*s);
    p++;
    q++;
    r++;
    s++;
  }
  assert(index.size()>0);
}
