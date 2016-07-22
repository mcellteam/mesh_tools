#include <stdio.h>
#include <math.h>
#include "meshproc.h"
#include "vector.h"
#include "polyvol.h"



/*
  Calculate volume of a material-specific Voronoi polyhedron
  by summing over all surface triangles of the polyhedron, the quantity:

              |  x1   x2   x3  |
    (1/6) det |  y1   y2   y3  |
              |  z1   z2   z3  |

  where (xn,yn,zn) is vertex n of a surface triangle.  The surface normals of
  all triangles of the polyhedron must point toward the outside of the
  polyhedron.
*/
double output_voronoi_element_by_material(struct tet_mesh *tet_mesh,
                                          struct node *node,
                                          u_int node_matindex)
{
  struct node **nodes;
  struct edge_list *elp;
  struct edge *ep;
  struct voronoi_facet *facet_head,*vfp;
  struct vector3 p1,p2,pint,v1,v2,v_node,v_norm;
  struct vector3 **vnp;
  double x1,x2,x3,y1,y2,y3,z1,z2,z3;
  double dot,det,volume;
  u_int n1;
  u_int n2;
  u_int n3;
  u_int v,nnodes,node_count;
  u_int edge_matindex;
  u_int i;
  byte done;

  nodes=tet_mesh->nodes;
  vnp=tet_mesh->voronoi_nodes;

  /* first find a point on the interior of the material-specific Voronoi cell */
  /* this is done by finding a material specific interior Delaunay edge  */
  /* and placing a point 1/4 of the way along this edge away from node */
  done=0;
  for (elp=node->node_material[node_matindex].edge_head; elp!=NULL && !done; elp=elp->next)
  {
    ep=elp->edge;

    if (ep->edge_type == INTERIOR)
    {
      if (node->node_index == ep->node_index[0])
      {
        p1.x=nodes[ep->node_index[0]]->x;
        p1.y=nodes[ep->node_index[0]]->y;
        p1.z=nodes[ep->node_index[0]]->z;
        p2.x=nodes[ep->node_index[1]]->x;
        p2.y=nodes[ep->node_index[1]]->y;
        p2.z=nodes[ep->node_index[1]]->z;
      } 
      else
      {
        p1.x=nodes[ep->node_index[1]]->x;
        p1.y=nodes[ep->node_index[1]]->y;
        p1.z=nodes[ep->node_index[1]]->z;
        p2.x=nodes[ep->node_index[0]]->x;
        p2.y=nodes[ep->node_index[0]]->y;
        p2.z=nodes[ep->node_index[0]]->z;
      }
      vectorize(&p1,&p2,&v1);
      pint.x=p1.x+0.25*v1.x;
      pint.y=p1.y+0.25*v1.y;
      pint.z=p1.z+0.25*v1.z;
      done=1;
    }
  }

  /* now calculate volume of Voronoi cell */
  volume = 0.0;
  for (elp=node->node_material[node_matindex].edge_head; elp!=NULL; elp=elp->next)
  {
    ep=elp->edge;

    if (ep->material_count>1)
    {
      edge_matindex=0;
      for (i=0;i<ep->material_count;i++)
      {
        if (node->node_material[node_matindex].matnum==ep->edge_material[i].matnum)
        {
          edge_matindex=i;
        }
      }
      facet_head=ep->edge_material[edge_matindex].facet;
    }
    else
    {
      facet_head=ep->facet;
    }

    for (vfp=facet_head; vfp!=NULL; vfp=vfp->next)
    {

    /* orient and triangulate each face.*/

      n3 = vfp->node_index[vfp->nnodes-1]; 

    /* determine orientation of Voronoi facet */
      /* construct vector from interior node to first vertex of facet */
      p1.x = vnp[vfp->node_index[0]]->x; 
      p1.y = vnp[vfp->node_index[0]]->y; 
      p1.z = vnp[vfp->node_index[0]]->z; 
      vectorize(&pint,&p1,&v_node);
      normalize(&v_node);

      /* construct normal vector of facet */
      p2.x = vnp[vfp->node_index[1]]->x; 
      p2.y = vnp[vfp->node_index[1]]->y; 
      p2.z = vnp[vfp->node_index[1]]->z; 
      vectorize(&p1,&p2,&v1);
      p2.x = vnp[n3]->x; 
      p2.y = vnp[n3]->y; 
      p2.z = vnp[n3]->z; 
      vectorize(&p1,&p2,&v2);
      cross_prod(&v1,&v2,&v_norm);
      normalize(&v_norm);

      dot=dot_prod(&v_node,&v_norm);

      nnodes=vfp->nnodes;
      if (dot>=0)
      {
        /* facet normal already points outward */
        n3 = vfp->node_index[nnodes-1]; 
        for (v=0; v<nnodes-2; v++)
        {

          n1 = vfp->node_index[v];
          n2 = vfp->node_index[v+1];
      
          x1=vnp[n1]->x;
          y1=vnp[n1]->y;
          z1=vnp[n1]->z;
          x2=vnp[n2]->x;
          y2=vnp[n2]->y;
          z2=vnp[n2]->z;
          x3=vnp[n3]->x;
          y3=vnp[n3]->y;
          z3=vnp[n3]->z;

          /* compute determinant of oriented triangle */
          det=x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1);

          volume+=det;
        }
      }
      else
      {
        /* facet normal points inward -- need to reverse it */
        n3 = vfp->node_index[1]; 
        for (v=0; v<nnodes-2; v++)
        {

          if (v==0)
          {
            n1 = vfp->node_index[0];
          }
          else
          {
            n1 = vfp->node_index[nnodes-v];
          }
          n2 = vfp->node_index[nnodes-v-1];
      
          x1=vnp[n1]->x;
          y1=vnp[n1]->y;
          z1=vnp[n1]->z;
          x2=vnp[n2]->x;
          y2=vnp[n2]->y;
          z2=vnp[n2]->z;
          x3=vnp[n3]->x;
          y3=vnp[n3]->y;
          z3=vnp[n3]->z;

          /* compute determinant of oriented triangle */
          det=x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1);

          volume+=det;
        }
      }
    }
  }
  volume = volume/6.0;

  return(volume); 
}



