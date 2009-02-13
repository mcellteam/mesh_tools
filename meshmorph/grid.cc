// Author: Justin Kinney
// Date: Sep 2008

#include "grid.h"

#include <cmath>

#include "face.h"


/** Compute barycenter of face tile.
 * \param[out] p Barycenter of face tile.
 * \param[in] index Index of face tile.
 * \param[in] n Not sure what this is.
 * \param[in] f Face of interes.
 */

void Grid::computeBarycenter (vector3 & p,int index,int n,Face * f)
{
  int root = static_cast<int>( sqrt(static_cast<double>(index)) );
  int rootrem = index - root*root;
  int k = n - root - 1;
  int j = rootrem/2;
  int i = rootrem - 2*j;

  double over3n = 1.0 / static_cast<double>(3*n);

  double ucoef = (static_cast<double>(3*j+i+1))*over3n*uv_vert1_u +
                 (static_cast<double>(3*k+i+1))*over3n*uv_vert2[0];
  double vcoef = (static_cast<double>(3*k+i+1))*over3n*uv_vert2[1];

  p.p[0] = ucoef*unit_u.p[0] + vcoef*unit_v.p[0] + *f->getVertex(0)->getCoord(0);
  p.p[1] = ucoef*unit_u.p[1] + vcoef*unit_v.p[1] + *f->getVertex(0)->getCoord(1);
  p.p[2] = ucoef*unit_u.p[2] + vcoef*unit_v.p[2] + *f->getVertex(0)->getCoord(2);
}

/** Create instance of Grid.
 * \param[in] f Face of interes.
 */

Grid::Grid (Face * f)
  :unit_u(),unit_v(),uv_vert1_u(0),uv_vert2()
{
  vector3 fv(*f->getVertex(1)->getPos() - *f->getVertex(0)->getPos());
  double ff = 1.0 / sqrt( fv.dot(fv) );
  unit_u = fv * ff;

  fv = *f->getVertex(2)->getPos() - *f->getVertex(0)->getPos();
  vector3 normal(unit_u.cross(fv));
  ff = 1.0 / sqrt( normal.dot(normal));
  normal *= ff;
  unit_v = normal.cross(unit_u);

  uv_vert1_u  = unit_u.dot(*f->getVertex(1)->getPos() - *f->getVertex(0)->getPos());
  uv_vert2[0] = unit_u.dot(*f->getVertex(2)->getPos() - *f->getVertex(0)->getPos());
  uv_vert2[1] = unit_v.dot(*f->getVertex(2)->getPos() - *f->getVertex(0)->getPos());
}

