#include <stdio.h>
#include "tetgen2FV.h"
#include "vector.h"
#include "polycent.h"

	
/* compute circumcenters of boundary face polygons in tet_mesh */
/* do this by finding the point of intersection of 3 planes. */
/* plane1: plane whose normal is one edge of boundary face and */
/*         bisects the edge. */
/* plane2: plane whose normal is another edge of boundary face and */
/*         bisects the edge. */
/* plane3: the plane of the boundary face itself */
/* the point of intersection of these 3 planes is the circumcenter */

void polycent(struct tet_mesh *tet_mesh)

{
  struct polygon_list *plp;
  struct polygon *pop;
  struct vector3 p1,p2,mid1,mid2;
  struct vector3 n1,n2,n3;
  struct vector3 xp1,xp2,xp3;
  double d1,d2,d3;
  double dot1;

  for (plp=tet_mesh->boundary_head; plp!=NULL; plp=plp->next)
  {
    pop=plp->polygon;

    p1.x=tet_mesh->nodes[pop->node_index[0]]->x;
    p1.y=tet_mesh->nodes[pop->node_index[0]]->y;
    p1.z=tet_mesh->nodes[pop->node_index[0]]->z;
    p2.x=tet_mesh->nodes[pop->node_index[1]]->x;
    p2.y=tet_mesh->nodes[pop->node_index[1]]->y;
    p2.z=tet_mesh->nodes[pop->node_index[1]]->z;
    mid1.x=0.5*(p1.x+p2.x);
    mid1.y=0.5*(p1.y+p2.y);
    mid1.z=0.5*(p1.z+p2.z);
    vectorize(&p1,&p2,&n1);
    normalize(&n1);
    d1=dot_prod(&mid1,&n1);

    p2.x=tet_mesh->nodes[pop->node_index[2]]->x;
    p2.y=tet_mesh->nodes[pop->node_index[2]]->y;
    p2.z=tet_mesh->nodes[pop->node_index[2]]->z;
    mid2.x=0.5*(p1.x+p2.x);
    mid2.y=0.5*(p1.y+p2.y);
    mid2.z=0.5*(p1.z+p2.z);
    vectorize(&p1,&p2,&n2);
    normalize(&n2);
    d2=dot_prod(&mid2,&n2);

    cross_prod(&n1,&n2,&n3);
    normalize(&n3);
    d3=dot_prod(&p1,&n3);

    cross_prod(&n2,&n3,&xp1);
    cross_prod(&n3,&n1,&xp2);
    cross_prod(&n1,&n2,&xp3);
    dot1=dot_prod(&n1,&xp1);

    xp1.x=d1*xp1.x;
    xp1.y=d1*xp1.y;
    xp1.z=d1*xp1.z;
    
    xp2.x=d2*xp2.x;
    xp2.y=d2*xp2.y;
    xp2.z=d2*xp2.z;

    xp3.x=d3*xp3.x;
    xp3.y=d3*xp3.y;
    xp3.z=d3*xp3.z;

    pop->cent.x=(xp1.x+xp2.x+xp3.x)/dot1;
    pop->cent.y=(xp1.y+xp2.y+xp3.y)/dot1;
    pop->cent.z=(xp1.z+xp2.z+xp3.z)/dot1;

  }

  return;
}

