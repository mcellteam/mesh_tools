#include <stdio.h>
#include <math.h>
#include "tetgen2FV.h"
#include "vector.h"
#include "polyvol.h"


/*
  Calculate volume of a Voronoi polyhedron
  This is done by creating a tetrahedralization of the Voronoi cell
  and then summing the volumes of all tets.
  The volume of each tet is given by:

    (1/6) * abs(dot_prod(a,cross_prod(b,c)))

  where a, b, and c are three edges of the tet that share a common
  vertex of the tet.
*/
double polyvol_1(struct tet_mesh *tet_mesh, struct node *node)
{
  struct node **nodes;
  struct edge_list *elp;
  struct edge *ep;
  struct voronoi_facet *vfp;
  struct vector3 p1,p2,a,b,c,bxc;
  struct vector3 **vnp;
  double volume,facet_volume;
  u_int i;

  nodes = tet_mesh->nodes;
  vnp=tet_mesh->voronoi_nodes;

  volume = 0.0;
  for (elp=node->edge_head; elp!=NULL; elp=elp->next)
  {
    ep=elp->edge;

    if ((node->node_type!=INTERFACE)
        || (node->node_type==INTERFACE && ep->edge_type!=VIRTUAL))
    {
      /* use the central node of the Voronoi cell as the apex of each tet */
      p1.x=node->x;
      p1.y=node->y;
      p1.z=node->z;

      /* Triangulate each Voronoi facet and use each triangle as */
      /* the base of each tet */
      for (vfp=ep->facet; vfp!=NULL; vfp=vfp->next)
      {
        facet_volume=0;

        for (i=2; i<vfp->nnodes; i++)
        {
          p2.x=vnp[vfp->node_index[0]]->x;
          p2.y=vnp[vfp->node_index[0]]->y;
          p2.z=vnp[vfp->node_index[0]]->z;
          vectorize(&p1,&p2,&a);

          p2.x=vnp[vfp->node_index[i-1]]->x;
          p2.y=vnp[vfp->node_index[i-1]]->y;
          p2.z=vnp[vfp->node_index[i-1]]->z;
          vectorize(&p1,&p2,&b);

          p2.x=vnp[vfp->node_index[i]]->x;
          p2.y=vnp[vfp->node_index[i]]->y;
          p2.z=vnp[vfp->node_index[i]]->z;
          vectorize(&p1,&p2,&c);

          cross_prod(&b,&c,&bxc);

          facet_volume += dot_prod(&a,&bxc);
        }
      }
      volume += fabs(facet_volume);
    }
  }

  volume = volume/6.0;
  
  return(volume);

}



/*
  Calculate the material-specific volume of a Voronoi polyhedron
  This is done by creating a tetrahedralization of the Voronoi cell
  and then summing the volumes of all tets.
  The volume of each tet is given by:

    (1/6) * abs(dot_prod(a,cross_prod(b,c)))

  where a, b, and c are three edges of the tet that share a common
  vertex of the tet.
*/
double polyvol_by_material_1(struct tet_mesh *tet_mesh,
                             struct node *node,
                             u_int node_matindex)
{
  struct node **nodes;
  struct edge_list *elp;
  struct edge *ep;
  struct voronoi_facet *facet_head,*vfp;
  struct vector3 p1,p2,a,b,c,bxc;
  struct vector3 **vnp;
  double volume,facet_volume;
  u_int edge_matindex;
  u_int i;

  nodes = tet_mesh->nodes;
  vnp=tet_mesh->voronoi_nodes;

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

    /* use the central node of the Voronoi cell as the apex of each tet */
    p1.x=node->x;
    p1.y=node->y;
    p1.z=node->z;

    /* Triangulate each Voronoi facet and use each triangle as */
    /* the base of each tet */
    for (vfp=facet_head; vfp!=NULL; vfp=vfp->next)
    {
      facet_volume=0;

      for (i=2; i<vfp->nnodes; i++)
      {
        p2.x=vnp[vfp->node_index[0]]->x;
        p2.y=vnp[vfp->node_index[0]]->y;
        p2.z=vnp[vfp->node_index[0]]->z;
        vectorize(&p1,&p2,&a);

        p2.x=vnp[vfp->node_index[i-1]]->x;
        p2.y=vnp[vfp->node_index[i-1]]->y;
        p2.z=vnp[vfp->node_index[i-1]]->z;
        vectorize(&p1,&p2,&b);

        p2.x=vnp[vfp->node_index[i]]->x;
        p2.y=vnp[vfp->node_index[i]]->y;
        p2.z=vnp[vfp->node_index[i]]->z;
        vectorize(&p1,&p2,&c);

        cross_prod(&b,&c,&bxc);

        facet_volume += dot_prod(&a,&bxc);
      }
    }
    volume += fabs(facet_volume);
  }

  volume = volume/6.0;
  
  return(volume);

}



/*
  Calculate volume of a Voronoi polyhedron
  by summing over all surface triangles of the polyhedron, the quantity:

              |  x1   x2   x3  |
    (1/6) det |  y1   y2   y3  |
              |  z1   z2   z3  |

  where (xn,yn,zn) is vertex n of a surface triangle.  The surface normals of
  all triangles of the polyhedron must point toward the outside of the
  polyhedron.
*/
double polyvol_2(struct tet_mesh *tet_mesh, struct node *node)
{
  struct node **nodes;
  struct edge_list *elp;
  struct edge *ep;
  struct voronoi_facet *vfp;
  struct vector3 p1,p2,pint,v1,v2,v_node,v_norm;
  struct vector3 **vnp;
  double x1,x2,x3,y1,y2,y3,z1,z2,z3;
  double dot,det,volume;
  double edge_ratio,l1,l2;
  u_int n1;
  u_int n2;
  u_int n3;
  u_int v,nnodes;
  u_int edge_node_index;
  u_int i;
  byte found;

  nodes = tet_mesh->nodes;
  vnp=tet_mesh->voronoi_nodes;

  /* first find a point on the interior of the Voronoi cell */
  /* this is done by finding an interior Delaunay edge  */
  /* and placing a point 1/4 of the way along this edge away from node */
  i=0;
  found=0;
  for (elp=node->edge_head; elp!=NULL && !found; elp=elp->next)
  {
    ep=elp->edge;

    if (ep->edge_type == INTERIOR)
    {
      if (node->node_index == ep->node_index[0])
      {
        edge_node_index=ep->node_index[1];
        p1.x=nodes[ep->node_index[0]]->x;
        p1.y=nodes[ep->node_index[0]]->y;
        p1.z=nodes[ep->node_index[0]]->z;
        p2.x=nodes[ep->node_index[1]]->x;
        p2.y=nodes[ep->node_index[1]]->y;
        p2.z=nodes[ep->node_index[1]]->z;
      } 
      else
      {
        edge_node_index=ep->node_index[0];
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
      found=1;
      if (node->node_index==1791)
      {
        printf("polyvol: node %u  using interior edge %u to %u\n",node->node_index,i,edge_node_index);
      }
    }
    i++;
  }

/*
  if (!found)
  {
    fprintf(stderr,"polyvol: warning: no interior point found for node %u\n",node->node_index);
  }
*/

  /* Set interior point to node if no suitable interior edge is found */
  if (!found)
  {
    pint.x=node->x;
    pint.y=node->y;
    pint.z=node->z;
  }


  i=0;
  volume = 0.0;

  for (elp=node->edge_head; elp!=NULL; elp=elp->next)
  {
    ep=elp->edge;

    if ((node->node_type!=INTERFACE)
        || (node->node_type==INTERFACE && ep->edge_type!=VIRTUAL))
    {
      for (vfp=ep->facet; vfp!=NULL; vfp=vfp->next)
      {

        /* orient and triangulate each facet.*/

        /* calculate facet edge length ratio of consecutive edges */
        /* find two consecutive edges with a reasonable ratio */
        /* this procedure avoids nearly coincident tet circumcenters */
        found=0;
        for (v=0; v<vfp->nnodes && !found; v++)
        {
          if (v==0)
          {
            n2 = vfp->node_index[v+1];
            n3 = vfp->node_index[vfp->nnodes-1];
          }
          else if (v==vfp->nnodes-1)
          {
            n2 = vfp->node_index[0];
            n3 = vfp->node_index[v-1];
          }
          else
          {
            n2 = vfp->node_index[v+1];
            n3 = vfp->node_index[v-1];
          }
          n1 = vfp->node_index[v]; 
          p1.x = vnp[n1]->x; 
          p1.y = vnp[n1]->y; 
          p1.z = vnp[n1]->z; 
          p2.x = vnp[n2]->x; 
          p2.y = vnp[n2]->y; 
          p2.z = vnp[n2]->z; 
          vectorize(&p1,&p2,&v1);
          p2.x = vnp[n3]->x; 
          p2.y = vnp[n3]->y; 
          p2.z = vnp[n3]->z; 
          vectorize(&p1,&p2,&v2);
          l1=vect_length(&v1);
          l2=vect_length(&v2);
          if (l1<l2)
          {
            edge_ratio=l1/l2;
          }
          else
          {
            edge_ratio=l2/l1;
          }
          if (edge_ratio>1e-9)
          {
            found=1;
          }
        }

/*
        if (!found)
        {
          fprintf(stderr,"polyvol: warning: no suitable edge length ratio found %u\n",node->node_index);
        }
*/

      /* determine orientation of Voronoi facet */
        /* construct vector from interior node to first vertex of facet */
        p1.x = vnp[n1]->x; 
        p1.y = vnp[n1]->y; 
        p1.z = vnp[n1]->z; 
/*
        p2.x=node->x;
        p2.y=node->y;
        p2.z=node->z;
*/
        p2.x=pint.x;
        p2.y=pint.y;
        p2.z=pint.z;
        vectorize(&p2,&p1,&v_node);
        normalize(&v_node);

        /* construct normal vector of facet */
        p2.x = vnp[n2]->x; 
        p2.y = vnp[n2]->y; 
        p2.z = vnp[n2]->z; 
        vectorize(&p1,&p2,&v1);
        p2.x = vnp[n3]->x; 
        p2.y = vnp[n3]->y; 
        p2.z = vnp[n3]->z; 
        vectorize(&p1,&p2,&v2);
        cross_prod(&v1,&v2,&v_norm);
        normalize(&v_norm);

        dot=dot_prod(&v_node,&v_norm);

        l1=vect_length(&v1);
        l2=vect_length(&v2);
        if (node->node_index==1791)
        {
          printf("polyvol: node %u facet %u n_voronoi_nodes %u dot %.9g l1 %.9g l2 %.9g\n",node->node_index,i,vfp->nnodes,dot,l1,l2);
        }

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
    i++;
  }
  volume = volume/6.0;

  return(volume); 
}



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
double polyvol_by_material_2(struct tet_mesh *tet_mesh,
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
  double l1,l2,edge_ratio;
  u_int n1;
  u_int n2;
  u_int n3;
  u_int v,nnodes,node_count;
  u_int edge_matindex;
  u_int i;
  byte found;

  nodes=tet_mesh->nodes;
  vnp=tet_mesh->voronoi_nodes;

  /* first find a point on the interior of the material-specific Voronoi cell */
  /* this is done by finding a material specific interior Delaunay edge  */
  /* and placing a point 1/4 of the way along this edge away from node */
  found=0;
  for (elp=node->node_material[node_matindex].edge_head; elp!=NULL && !found; elp=elp->next)
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
      found=1;
    }
  }
/*
  if (!found)
  {
    fprintf(stderr,"polyvol_by_material: warning: no interior point found for node %u\n",node->node_index);
  }
*/

  /* NOTE: the following is disabled in favor of the above method */
  /* Use the following only if no suitable interior edge is found */
  /* first find a point on the interior of the material-specific Voronoi cell */
  /* this is done by calculating the average value of the Voronoi vertices */
  if (!found)
  {
    pint.x=0;
    pint.y=0;
    pint.z=0;
    node_count=0;
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
        nnodes=vfp->nnodes;
        for (v=0; v<nnodes; v++)
        {
          n1 = vfp->node_index[v];
          pint.x+=vnp[n1]->x;
          pint.y+=vnp[n1]->y;
          pint.z+=vnp[n1]->z;
          node_count++;
        }
  
      }
    }
    pint.x=pint.x/node_count;
    pint.y=pint.y/node_count;
    pint.z=pint.z/node_count;
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
        /* calculate facet edge length ratio of consecutive edges */
        /* find two consecutive edges with a reasonable ratio */
        /* this procedure avoids nearly coincident tet circumcenters */
        found=0;
        for (v=0; v<vfp->nnodes && !found; v++)
        {
          if (v==0)
          {
            n2 = vfp->node_index[v+1];
            n3 = vfp->node_index[vfp->nnodes-1];
          }
          else if (v==vfp->nnodes-1)
          {
            n2 = vfp->node_index[0];
            n3 = vfp->node_index[v-1];
          }
          else
          {
            n2 = vfp->node_index[v+1];
            n3 = vfp->node_index[v-1];
          }
          n1 = vfp->node_index[v]; 
          p1.x = vnp[n1]->x; 
          p1.y = vnp[n1]->y; 
          p1.z = vnp[n1]->z; 
          p2.x = vnp[n2]->x; 
          p2.y = vnp[n2]->y; 
          p2.z = vnp[n2]->z; 
          vectorize(&p1,&p2,&v1);
          p2.x = vnp[n3]->x; 
          p2.y = vnp[n3]->y; 
          p2.z = vnp[n3]->z; 
          vectorize(&p1,&p2,&v2);
          l1=vect_length(&v1);
          l2=vect_length(&v2);
          if (l1<l2)
          {
            edge_ratio=l1/l2;
          }
          else
          {
            edge_ratio=l2/l1;
          }
          if (edge_ratio>1e-9)
          {
            found=1;
          }
        }

    /* determine orientation of Voronoi facet */
      /* construct vector from interior node to first vertex of facet */
      p1.x = vnp[n1]->x; 
      p1.y = vnp[n1]->y; 
      p1.z = vnp[n1]->z; 
      vectorize(&pint,&p1,&v_node);
      normalize(&v_node);

      /* construct normal vector of facet */
      p2.x = vnp[n2]->x; 
      p2.y = vnp[n2]->y; 
      p2.z = vnp[n2]->z; 
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



/*
  Calculate volume of a tet_mesh by summing the volumes of all tets.
  Volume of each tet is given by:

    (1/6) * abs(dot_prod(a,cross_prod(b,c)))

  where a, b, and c are three edges of the tet that share a common
  vertex of the tet.
*/
double meshvol(struct tet_mesh *tet_mesh)
{


  struct node **nodes;
  struct tet_list *tlp;
  struct tet *tp;
  struct vector3 p1,p2,a,b,c,bxc;
  double volume;

  nodes=tet_mesh->nodes;

  volume = 0.0;

  for (tlp=tet_mesh->tet_head; tlp!=NULL; tlp=tlp->next)
  {
    tp=tlp->tet;
    p1.x=nodes[tp->node_index[0]]->x;
    p1.y=nodes[tp->node_index[0]]->y;
    p1.z=nodes[tp->node_index[0]]->z;

    p2.x=nodes[tp->node_index[1]]->x;
    p2.y=nodes[tp->node_index[1]]->y;
    p2.z=nodes[tp->node_index[1]]->z;
    vectorize(&p1,&p2,&a);

    p2.x=nodes[tp->node_index[2]]->x;
    p2.y=nodes[tp->node_index[2]]->y;
    p2.z=nodes[tp->node_index[2]]->z;
    vectorize(&p1,&p2,&b);

    p2.x=nodes[tp->node_index[3]]->x;
    p2.y=nodes[tp->node_index[3]]->y;
    p2.z=nodes[tp->node_index[3]]->z;
    vectorize(&p1,&p2,&c);

    cross_prod(&b,&c,&bxc);

    volume += fabs(dot_prod(&a,&bxc));
  }

  volume = volume/6.0;

  return(volume); 
}



#if 0

/* find a point inside a possibly concave voronoi polyhedron */
/* this is done by ray-tracing from a point to edge of the polyhedron */
/* and counting intersections with the polygons of the polyhedron */
int point_in_polyhedron(struct tet_mesh *tet_mesh, struct node *node, int node_matindex, struct vector3 *p)
{
}



/* trace a ray from pnt1 to pnt2 */
/* and check for intersection with a single piece of a voronoi facet */
/* return 0 on no intersection and 1 on intersection */
int ray_trace(struct tet_mesh *tet_mesh,
              struct voronoi_facet *vfp,
              struct vector3 *pnt1,
              struct vector3 *pnt2)
{
  struct vector3 **vnp;

  vnp=tet_mesh->voronoi_nodes;

  /* calculate facet edge length ratio of consecutive edges */
  /* find two consecutive edges with a reasonable ratio */
  /* this procedure avoids nearly coincident facet vertices */
  found=0;
  for (v=0; v<vfp->nnodes && !found; v++)
  {
    if (v==0)
    {
      n2 = vfp->node_index[v+1];
      n3 = vfp->node_index[vfp->nnodes-1];
    }
    else if (v==vfp->nnodes-1)
    {
      n2 = vfp->node_index[0];
      n3 = vfp->node_index[v-1];
    }
    else
    {
      n2 = vfp->node_index[v+1];
      n3 = vfp->node_index[v-1];
    }
    n1 = vfp->node_index[v];
    p1.x = vnp[n1]->x;
    p1.y = vnp[n1]->y;
    p1.z = vnp[n1]->z;
    p2.x = vnp[n2]->x;
    p2.y = vnp[n2]->y;
    p2.z = vnp[n2]->z;
    vectorize(&p1,&p2,&v1);
    p2.x = vnp[n3]->x;
    p2.y = vnp[n3]->y;
    p2.z = vnp[n3]->z;
    vectorize(&p1,&p2,&v2);
    l1=vect_length(&v1);
    l2=vect_length(&v2);
    if (l1<l2)
    {
      edge_ratio=l1/l2;
    }
    else
    {
      edge_ratio=l2/l1;
    }
    if (edge_ratio>1e-9)
    {
      found=1;
    }
  }

  /* construct normal vector of facet */
  p1.x = vnp[n1]->x;
  p1.y = vnp[n1]->y;
  p1.z = vnp[n1]->z;
  p2.x = vnp[n2]->x;
  p2.y = vnp[n2]->y;
  p2.z = vnp[n2]->z;
  vectorize(&p1,&p2,&v1);
  p2.x = vnp[n3]->x;
  p2.y = vnp[n3]->y;
  p2.z = vnp[n3]->z;
  vectorize(&p1,&p2,&v2);
  cross_prod(&v1,&v2,&v_norm);
  normalize(&v_norm);

  a=v_norm.x;
  b=v_norm.y;
  c=v_norm.z;
  d=(v_norm.x*p1.x)+(v_norm.y*p1.y)+(v_norm.z*p1.z);
  projection=0;
  tmp=v_norm.x;
  if (tmp*tmp < v_norm.y*v_norm.y)
  {
    projection=1;
    tmp=v_norm.y;
  }
  if (tmp*tmp < v_norm.z*v_norm.z)
  {
    projection=2;
  }


  n_dot_delta=a*dx+b*dy+c*dz;
  if (n_dot_delta==0)
  {
    t=-1000;
  }
  else
  {
    t=(d-a*x-b*y-c*z)/n_dot_delta;
  }

  if (t>=0.0)
  {
    /* take 2D projection of polygon and hit point
    * and determine if hit point is within the projected polygon boundaries */
    x=pnt1->x;
    y=pnt1->y;
    z=pnt1->z;
    dx=pnt2->x-pnt1->x;
    dy=pnt2->y-pnt1->y;
    dz=pnt2->z-pnt1->z;
    vert=vnp[vfp->node_index[0]];
    if (projection == 0)
    {
      h1=y+t*dy;
      h2=z+t*dz;
      a1=vert->y-h1;
      a2=vert->z-h2;
    }
    else if (projection == 1)
    {
      h1=x+t*dx;
      h2=z+t*dz;
      a1=vert->x-h1;
      a2=vert->z-h2;
    }
    else
    {
      h1=x+t*dx;
      h2=y+t*dy;
      a1=vert->x-h1;
      a2=vert->y-h2;
    }
    j=1;
    prev_det=0;
    inside=1;
    n_vert = wp->n_vert;
    while (inside && j<vfp->nnodes;)
    {
      vert=vnp[vfp->node_index[j]];
      if (projection == 0)
      {
        b1=vert->y-h1;
        b2=vert->z-h2;
      }
      else if (projection == 1)
      {
        b1=vert->x-h1;
        b2=vert->z-h2;
      }
      else
      {
        b1=vert->x-h1;
        b2=vert->y-h2;
      }
      det=a1*b2-a2*b1;
      a1=b1;
      a2=b2;
      if (j>1)
      {
        inside=(!((det<0 && prev_det>=0) || (det>=0 && prev_det<0)));
      }
      prev_det=det;
      j++;
    }
    if (inside)
    {
      vert=vnp[vfp->node_index[0]];
      if (projection == 0)
      {
        b1=vert->y-h1;
        b2=vert->z-h2;
      }
      else if (projection == 1)
      {
        b1=vert->x-h1;
        b2=vert->z-h2;
      }
      else
      {
        b1=vert->x-h1;
        b2=vert->y-h2;
      }
      det=a1*b2-a2*b1;
      inside=(!((det<0 && prev_det>=0) || (det>=0 && prev_det<0)));
      prev_det=det;
    }

    if (inside)
    {
      return(1);
    }

  }

  return(0);
}

#endif /* 0 */
