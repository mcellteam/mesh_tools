#include <math.h> 
#include "tetgen2FV.h"
#include "polyarea.h"

double	polyarea(struct tet_mesh *tet_mesh, struct voronoi_facet *vfp)

{


  double area;
  u_int j;
  u_int k1;
  u_int k2;
  double v[3];

  area = 0.0;

    v[0]=0.0;
    v[1]=0.0;
    v[2]=0.0;

 /* For each triangle in the face, compute the normal vector*/

     for (j=0; j<vfp->nnodes; j++)
     {

        k1 = vfp->node_index[j];

        if ( j < vfp->nnodes-1) 
        {
           k2 = vfp->node_index[j+1];
        }
        else
        {
          k2 = vfp->node_index[0];
        }

 /* Compute the cross product */

      v[0] = v[0]+  tet_mesh->voronoi_nodes[k1]->y * tet_mesh->voronoi_nodes[k2]->z - 
           tet_mesh->voronoi_nodes[k1]->z*tet_mesh->voronoi_nodes[k2]->y; 
      v[1] = v[1] + tet_mesh->voronoi_nodes[k1]->z* tet_mesh->voronoi_nodes[k2]->x - 
           tet_mesh->voronoi_nodes[k1]->x*tet_mesh->voronoi_nodes[k2]->z;
      v[2] = v[2] + tet_mesh->voronoi_nodes[k1]->x*tet_mesh->voronoi_nodes[k2]->y - 
           tet_mesh->voronoi_nodes[k1]->y * tet_mesh->voronoi_nodes[k2]->x;
     }

 /*  Add the magnitude of the normal vector to the sum.  */

    area = 0.5*sqrt (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    return(area);
}


