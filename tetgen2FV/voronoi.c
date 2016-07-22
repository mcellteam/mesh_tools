#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "tetgen2FV.h"
#include "polycent.h"
#include "c_3d.h"
#include "polyarea.h"
#include "polyvol.h"
#include "strfunc.h"
#include "voronoi.h"


/* given 3D Delaunay (tet) mesh construct 3D Voronoi dual mesh */

int construct_voronoi_mesh(struct tet_mesh *tet_mesh)
{
  struct tet *tp;
  struct tet_list *tlp;
  struct node *np;
  struct edge *ep,*right_ep,*left_ep;
  struct edge_list *elp;
  struct hash_table *hp;
  struct voronoi_volume_element *vvep;
  struct voronoi_facet *vfp;
  struct vector3 **voronoi_nodes;
  struct vector3 v1,p1,p2;
  struct polygon_list *plp;
  struct polygon *pop;
  double tot_mesh_volume,tot_voronoi_volume,boundary_area;
  int found_node;
  u_int rn,ln,tn1,tn2,itmp;
  u_int n_voronoi_nodes,voronoi_node_index_offset;
  u_int n_voronoi_facet_nodes;
  u_int i,j,k;
  char ckey[17];
  byte prev_tet_type,not_found,boundary_tet,boundary_tet_next;

  fprintf(stderr,"tetgen2FV: constructing Voronoi mesh:\n");

  tot_mesh_volume=0;
  tot_voronoi_volume=0;
  boundary_area=0;

  for (i=0; i<17; i++)
  {
    ckey[i]=0x00;
  }

  /* calculate circumcenters of boundary faces in tet_mesh */
  fprintf(stderr,"          calculating boundary face circumcenters...\n");
  polycent(tet_mesh);

  /* now calculate circumcenters of tets in tet_mesh */
  fprintf(stderr,"          calculating tet circumcenters...\n");
  if (circumcent(tet_mesh))
  {
    fprintf(stderr,
      "tetgen2FV: fatal error while computing tet_mesh circumcenters\n");
    return(1);
  }

  /* next, construct edge neighbor list of Delaunay mesh
     from which Voronoi mesh will be built */
  fprintf(stderr,"          constructing edge neighbor list...\n");
  if (construct_edge_list(tet_mesh))
  {
    fprintf(stderr,"tetgen2FV: fatal error while constructing edge list\n");
    return(1);
  }

  /* associate boundary faces with their parent nodes and edges */
  fprintf(stderr,"          constructing node boundary face list...\n");
  if (construct_node_boundary_face_list(tet_mesh))
  {
    fprintf(stderr,
      "tetgen2FV: fatal error while constructing node boundary face list\n");
    return(1);
  }

  /* allocate array of pointers to circumcenters for Voronoi nodes */
  n_voronoi_nodes=tet_mesh->ntets+tet_mesh->nboundaryfaces+tet_mesh->nboundaryedges+tet_mesh->nboundarynodes;
  tet_mesh->nvoronoinodes=n_voronoi_nodes;

  if ((voronoi_nodes=(struct vector3 **)malloc
     (n_voronoi_nodes*sizeof(struct vector3 *)))==NULL)
  {
    fprintf(stderr,
      "tetgen2FV: out of memory while storing Voronoi mesh nodes\n");
    return(1);
  }
  tet_mesh->voronoi_nodes=voronoi_nodes;

  /* map tet circumcenters to Voronoi_nodes */
  for (tlp=tet_mesh->tet_head; tlp!=NULL; tlp=tlp->next)
  {
    tp=tlp->tet;
    voronoi_nodes[tp->tet_index]=&(tp->cent);
  }
  voronoi_node_index_offset=tet_mesh->ntets;

  /* map boundary face circumcenters to Voronoi_nodes */
  for (plp=tet_mesh->boundary_head; plp!=NULL; plp=plp->next)
  {
    pop=plp->polygon;
    
    pop->voronoi_index=voronoi_node_index_offset;
    voronoi_nodes[pop->voronoi_index]=&(pop->cent);
    voronoi_node_index_offset++;
  }

  /* map boundary edge circumcenters to Voronoi_nodes */
  for (i=0; i<tet_mesh->hashsize; i++)
  {
    for (hp=tet_mesh->edge_hashtab[i]; hp!=NULL; hp=hp->next)
    {
      ep=(struct edge *)hp->contents;
      if (ep->edge_type==BOUNDARY || ep->edge_type==INTERFACE)
      {
        ep->voronoi_index=voronoi_node_index_offset;
        voronoi_nodes[ep->voronoi_index]=ep->cent;
        voronoi_node_index_offset++;
      }
    }
  }

  /* map boundary nodes to Voronoi_nodes */
  for (i=0; i<tet_mesh->nnodes; i++)
  {
    np=tet_mesh->nodes[i];
    if (np->node_type==BOUNDARY || np->node_type==INTERFACE)
    {
      np->voronoi_index=voronoi_node_index_offset;
      voronoi_nodes[np->voronoi_index]=(struct vector3 *)np;
      voronoi_node_index_offset++;
    }
  }


  /* now construct Voronoi mesh */
  /* by traversing each node in Delaunay tet_mesh */
  fprintf(stderr,"          constructing Voronoi cells...\n");
  for (i=0; i<tet_mesh->nnodes; i++)
  {
    np=tet_mesh->nodes[i];

    /* allocate memory for Voronoi cell associated with this node */
    if ((vvep=(struct voronoi_volume_element *)malloc
       (sizeof(struct voronoi_volume_element)))==NULL)
    {
      fprintf(stderr,
        "tetgen2FV: out of memory while constructing Voronoi volume element\n");
      return(1);
    }
    vvep->volume=0;
    vvep->nedges=0;
    vvep->nedges_tot=0;
    np->volume_element=vvep;

    /* given a list of Delaunay edges connected to this Delaunay node
       construct Voronoi cell around Delaunay node */
    for (elp=np->edge_head; elp!=NULL; elp=elp->next)
    {
      ep=elp->edge;

      if (ep->edge_type==INTERIOR || ep->edge_type==INTERFACE)
      {
        n_voronoi_facet_nodes=ep->ntets;
      }
      else if (ep->edge_type==BOUNDARY)
      {
        n_voronoi_facet_nodes=ep->ntets+3;
      }
      else if (ep->edge_type==VIRTUAL)
      {
        n_voronoi_facet_nodes=(2*np->nboundaryfaces)+1;
      }

      /* construct Voronoi facet for this edge
         only if this has not already been done
         and only if facet is not degenerate */
      if (n_voronoi_facet_nodes < 3)
      {
        fprintf(stderr,"tetgen2FV: oops! should never happen: degenerate Voronoi facet found at node %d: node_type %d, edge_type %d, ntets %d\n",np->node_index,np->node_type,ep->edge_type,ep->ntets);
      }
      if (ep->facet==NULL && n_voronoi_facet_nodes > 2)
      {
        if (ep->edge_type!=VIRTUAL)
        {
          /* make the one Voronoi facet associated with this edge --
             the facet is a planar-polygon whose vertices are the
             tet circumcenter points associated with this edge */

          /* allocate memory for the Voronoi facet associated with this edge */
          if ((vfp=(struct voronoi_facet *)malloc
             (sizeof(struct voronoi_facet)))==NULL)
          {
            fprintf(stderr,
             "tetgen2FV: out of memory while constructing Voronoi volume element\n");
            return(1);
          }
          vfp->area=0;
          vfp->nnodes=n_voronoi_facet_nodes;
          tet_mesh->voronoi_facet_ntriangles_tot+=(vfp->nnodes-2);
          vvep->nedges_tot+=vfp->nnodes;
          tet_mesh->voronoi_facet_nedges_tot+=vfp->nnodes;
          vfp->node_index=NULL;
          vfp->next=NULL;
          ep->facet=vfp;
          np->nfacets++;
  
          if ((vfp->node_index=(u_int *)malloc
             (vfp->nnodes*sizeof(u_int)))==NULL)
          {
            fprintf(stderr,
             "tetgen2FV: out of memory while constructing Voronoi volume element\n");
            return(1);
          }
  
          /* given perimeter-ordered list of tets that share this Delaunay edge
             construct Voronoi cell facet intersected by this Delaunay edge */
          if (ep->edge_type==INTERIOR || ep->edge_type==INTERFACE)
          {
            /* circumcenter points form complete perimeter around edge */
            /* construct perimeter of Voronoi facet
               one circumcenter point at a time */
            j=0;
            for (tlp=ep->tet_head; tlp!=NULL; tlp=tlp->next)
            {
              tp=tlp->tet;

              /* include tet circumcenter */
              vfp->node_index[j++]=tp->tet_index;
            }
          }
          else if (ep->edge_type==BOUNDARY)
          {
            if (ep->ntets==1)
            {
              /* circumcenter points form partial perimeter around edge --
                 close the perimeter by using additional points lying on the
                 boundary faces sharing this edge and the edge mid-point */
              j=0;
              vfp->node_index[j++]=ep->boundary_face[0]->voronoi_index;
              vfp->node_index[j++]=ep->tet_head->tet->tet_index;
              vfp->node_index[j++]=ep->boundary_face[1]->voronoi_index;
              vfp->node_index[j++]=ep->voronoi_index;
            }
            else if (ep->ntets==2)
            {
              /* circumcenter points form partial perimeter around edge --
                 close the perimeter by using additional points lying on the
                 boundary faces sharing this edge and the edge mid-point */
              j=0;
              for (tlp=ep->tet_head; tlp!=NULL; tlp=tlp->next)
              {
                tp=tlp->tet;

                if (tlp==ep->tet_head)
                {
                  /* include tet circumcenter before additional points */
                  vfp->node_index[j++]=tp->tet_index;

                  /* determine perimeter order of boundary faces and
                     include additional points lying on the boundary faces
                     sharing this edge and the edge mid_point */
                  /* find boundary face of this edge
                     which is also shared by this tet */
                  not_found=1;
                  for (k=0; k<2 && not_found; k++)
                  {
                    pop=ep->boundary_face[k];
                    if (pop->shared_tet[0]==tp || pop->shared_tet[1]==tp)
                    {
                      not_found=0;
                      if (k==0)
                      {
                        vfp->node_index[j++]=ep->boundary_face[0]->voronoi_index;
                        vfp->node_index[j++]=ep->voronoi_index;
                        vfp->node_index[j++]=ep->boundary_face[1]->voronoi_index;
                      }
                      else
                      {
                        vfp->node_index[j++]=ep->boundary_face[1]->voronoi_index;
                        vfp->node_index[j++]=ep->voronoi_index;
                        vfp->node_index[j++]=ep->boundary_face[0]->voronoi_index;
                      }
                    }
                  }
                }
                else
                {
                  /* include tet circumcenter */
                  vfp->node_index[j++]=tp->tet_index;
                }
              }
            }
            else if (ep->ntets > 2)
            {

              /* circumcenter points form partial perimeter around edge --
                 close the perimeter by using additional points lying on the
                 boundary faces sharing this edge and the edge mid-point */
              j=0;
              prev_tet_type=INTERIOR;
              for (tlp=ep->tet_head; tlp!=NULL; tlp=tlp->next)
              {
                tp=tlp->tet;

                boundary_tet=0;
                if (tp->tet_type==BOUNDARY)
                {
                  /* boundary status only valid
                     if tet and edge share boundary face */
                  /* find boundary face of this edge
                     which is also shared by this tet, if there is one */
                  not_found=1;
                  for (k=0; k<2 && not_found; k++)
                  {
                    pop=ep->boundary_face[k];
                    if (pop->shared_tet[0]==tp || pop->shared_tet[1]==tp)
                    {
                      not_found=0;
                      boundary_tet=1;
                    }
                  }
                }

                if (!boundary_tet)
                {
                  /* include tet circumcenter */
                  vfp->node_index[j++]=tp->tet_index;
                }
                if (boundary_tet && prev_tet_type==BOUNDARY)
                {
                  /* include tet circumcenter */
                  vfp->node_index[j++]=tp->tet_index;
                }
                if (boundary_tet && tlp==ep->tet_head)
                {
                  /* include tet circumcenter before additional points */
                  vfp->node_index[j++]=tp->tet_index;

                  boundary_tet_next=0;
                  if (tlp->next->tet->tet_type==BOUNDARY)
                  {
                    /* boundary status only valid
                       if next tet and edge share boundary face */
                    /* find boundary face of this edge
                       which is also shared by next tet, if there is one */
                    not_found=1;
                    for (k=0; k<2 && not_found; k++)
                    {
                      pop=ep->boundary_face[k];
                      if (pop->shared_tet[0]==tlp->next->tet
                          || pop->shared_tet[1]==tlp->next->tet)
                      {
                        not_found=0;
                        boundary_tet_next=1;
                      }
                    }
                  }

                  if (boundary_tet_next)
                  {
                    /* determine perimeter order of boundary faces and
                       include additional points lying on the boundary faces
                       sharing this edge and the edge mid_point */
                    /* find boundary face of this edge
                       which is also shared by this tet */
                    not_found=1;
                    for (k=0; k<2 && not_found; k++)
                    {
                      pop=ep->boundary_face[k];
                      if (pop->shared_tet[0]==tp || pop->shared_tet[1]==tp)
                      {
                        not_found=0;
                        if (k==0)
                        {
                          vfp->node_index[j++]=ep->boundary_face[0]->voronoi_index;
                          vfp->node_index[j++]=ep->voronoi_index;
                          vfp->node_index[j++]=ep->boundary_face[1]->voronoi_index;
                        }
                        else
                        {
                          vfp->node_index[j++]=ep->boundary_face[1]->voronoi_index;
                          vfp->node_index[j++]=ep->voronoi_index;
                          vfp->node_index[j++]=ep->boundary_face[0]->voronoi_index;
                        }
                      }
                    }
                  }
                }
                if (boundary_tet
                    && tlp!=ep->tet_head
                    && prev_tet_type!=BOUNDARY)
                {
                  /* include tet circumcenter before additional points */
                  vfp->node_index[j++]=tp->tet_index;

                  /* determine perimeter order of boundary faces and
                     include additional points lying on the boundary faces
                     sharing this edge and the edge mid_point */
                  /* find boundary face of this edge
                     which is also shared by this tet */
                  not_found=1;
                  for (k=0; k<2 && not_found; k++)
                  {
                    pop=ep->boundary_face[k];
                    if (pop->shared_tet[0]==tp || pop->shared_tet[1]==tp)
                    {
                      not_found=0;
                      if (k==0)
                      {
                        vfp->node_index[j++]=ep->boundary_face[0]->voronoi_index;
                        vfp->node_index[j++]=ep->voronoi_index;
                        vfp->node_index[j++]=ep->boundary_face[1]->voronoi_index;
                      }
                      else
                      {
                        vfp->node_index[j++]=ep->boundary_face[1]->voronoi_index;
                        vfp->node_index[j++]=ep->voronoi_index;
                        vfp->node_index[j++]=ep->boundary_face[0]->voronoi_index;
                      }
                    }
                  }
                }
                if (boundary_tet)
                {
                  prev_tet_type=BOUNDARY;
                }
                else
                {
                  prev_tet_type=INTERIOR;
                }

              }

            }
          }

          /* calculate area of Voronoi facet (dA) */
	  vfp->area=polyarea(tet_mesh,vfp);
          if (vfp->area==0)
          {
            fprintf(stderr,"warning: zero Voronoi facet area at node: %d %g\n",i,vfp->area);
          }
/*
          printf("%d %g\n",i,vfp->area);
*/
        }
        else if (ep->edge_type==VIRTUAL)
        {
          /* make a group of Voronoi facets to divide an interface cell */
          /* or to cap off a boundary cell */
          /* each facet in the group is a kite-shaped polygon with 4 vertices */
 
          /* make one Voronoi facet for each boundary face at this node */
          for (plp=np->boundary_head; plp!=NULL; plp=plp->next)
          {
            pop=plp->polygon;

            /* allocate memory for a Voronoi facet associated with this edge */
            if ((vfp=(struct voronoi_facet *)malloc
               (sizeof(struct voronoi_facet)))==NULL)
            {
              fprintf(stderr,
               "tetgen2FV: out of memory while constructing Voronoi volume element\n");
              return(1);
            }
            vfp->area=0;
            vfp->nnodes=4;
            tet_mesh->voronoi_cap_ntriangles_tot+=2;
            vvep->nedges_tot+=vfp->nnodes;
            tet_mesh->voronoi_cap_nedges_tot+=vfp->nnodes;
            vfp->node_index=NULL;
            vfp->next=ep->facet;
            ep->facet=vfp;
            np->nfacets++;
  
            if ((vfp->node_index=(u_int *)malloc
               (vfp->nnodes*sizeof(u_int)))==NULL)
            {
              fprintf(stderr,
               "tetgen2FV: out of memory while constructing Voronoi volume element\n");
              return(1);
            }

            /* make kite-shaped facet */
            /* facet is composed of 4 vertices:
                1) this node
                2) mid-point of right-edge of boundary face
                3) centroid of boundary face
                4) mid-point of left-edge of boundary face */

            /* find right and left edges of boundary face: */
            /* first find index of this node in boundary face node_index */
            found_node=-1;
            for (j=0; (j<3)&&(found_node<0); j++)
            {
              if (pop->node_index[j]==np->node_index)
              {
                found_node=j;
              }
            }

            if (found_node==0)
            {
              rn=1;
              ln=2;
            }
            else if (found_node==1)
            {
              rn=2;
              ln=0;
            }
            else if (found_node==2)
            {
              rn=0;
              ln=1;
            }

            /* find right edge: */
            /* sort edge nodes into increasing order */
            tn1 = pop->node_index[found_node];
            tn2 = pop->node_index[rn];
            if (tn1>tn2)
            {
              itmp=tn1;
              tn1=tn2;
              tn2=itmp;
            }

            /* construct edge key for this edge */
            sprintf(ckey,"%08x%08x",tn1,tn2);
     
            /* look-up edge key in edge hash table */
            hp=retrieve_key(ckey,tet_mesh->hashmask,tet_mesh->edge_hashtab);
            right_ep=(struct edge *)hp->contents;

            /* find left edge: */
            /* sort edge nodes into increasing order */
            tn1 = pop->node_index[found_node];
            tn2 = pop->node_index[ln];
            if (tn1>tn2)
            {
              itmp=tn1;
              tn1=tn2;
              tn2=itmp;
            }

            /* construct edge key for this edge */
            sprintf(ckey,"%08x%08x",tn1,tn2);
     
            /* look-up edge key in edge hash table */
            hp=retrieve_key(ckey,tet_mesh->hashmask,tet_mesh->edge_hashtab);
            left_ep=(struct edge *)hp->contents;

            /* now make kite-shaped facet: */
            vfp->node_index[0]=np->voronoi_index;
            vfp->node_index[1]=right_ep->voronoi_index;
            vfp->node_index[2]=pop->voronoi_index;
            vfp->node_index[3]=left_ep->voronoi_index;

            /* calculate area of voronoi facet (dA) */
	    vfp->area=polyarea(tet_mesh,vfp);
            boundary_area += vfp->area;
            if (vfp->area==0)
            {
              fprintf(stderr,"warning: zero Voronoi facet area at node: %d %g\n",i,vfp->area);
            }
/*
            printf("%d %g\n",i,vfp->area);
*/
          }
        }
  

        /* calculate length of Delaunay edge normal to the Voronoi facet (dS) */
        p1.x=tet_mesh->nodes[ep->node_index[0]]->x;
        p1.y=tet_mesh->nodes[ep->node_index[0]]->y;
        p1.z=tet_mesh->nodes[ep->node_index[0]]->z;
        p2.x=tet_mesh->nodes[ep->node_index[1]]->x;
        p2.y=tet_mesh->nodes[ep->node_index[1]]->y;
        p2.z=tet_mesh->nodes[ep->node_index[1]]->z;
        vectorize(&p1,&p2,&v1);
        ep->length=vect_length(&v1);
/*
        if (ep->length==0)
        {
          fprintf(stderr,"warning: zero length of Delaunay edge at node: %d %g\n",i,ep->length);
        }
*/

      }

    }

    tet_mesh->voronoi_nfacets_tot+=np->nfacets;
 
    /* calculate volume of Voronoi cell (dV) */
    vvep->volume=polyvol_1(tet_mesh,np);
    if (vvep->volume<=0)
    {
      fprintf(stderr,"warning: negative or zero volume for Voronoi element at node: %d %g\n",i,vvep->volume);
    }
    tot_voronoi_volume += vvep->volume;
/*
    printf("%d %g\n",i,vvep->volume);
*/

  }

  /* now construct material-specific Voronoi meshes */
  if (construct_material_info(tet_mesh))
  {
    fprintf(stderr,"tetgen2FV: fatal error while constructing material info\n");
    return(1);
  }

/*
  check_normals(tet_mesh);
*/

  tot_mesh_volume=meshvol(tet_mesh);
  printf("total Voronoi volume = %.17g\n",tot_voronoi_volume);
  printf("total tet volume = %.17g\n",tot_mesh_volume);
  printf("boundary surface area = %.17g\n",boundary_area);

/*
  output_voronoi_element(tet_mesh,tet_mesh->nodes[77551]);
  output_voronoi_element(tet_mesh,tet_mesh->nodes[354346]);
  output_voronoi_element_by_material(tet_mesh,tet_mesh->nodes[354346],1);
  output_voronoi_element_dx(tet_mesh,tet_mesh->nodes[77551],"voronoi_element_77551.dx");
  output_voronoi_element_dx(tet_mesh,tet_mesh->nodes[354350],"voronoi_element_354350.dx");
  output_voronoi_element_by_material(tet_mesh,tet_mesh->nodes[1791],0);
  output_voronoi_element_dx(tet_mesh,tet_mesh->nodes[1791],"voronoi_element_1791.dx");
  output_voronoi_element_by_material(tet_mesh,tet_mesh->nodes[26492],1);
  output_voronoi_element_dx(tet_mesh,tet_mesh->nodes[26492],"voronoi_element_26492.dx");
  output_voronoi_element_dx(tet_mesh,tet_mesh->nodes[460594],"voronoi_element_460594.dx");
*/

  return(0);
}



/* associate boundary faces with their parent nodes, edges, and tets */
int construct_node_boundary_face_list(struct tet_mesh *tet_mesh)
{
  struct polygon_list *plp,*new_plp;
  struct polygon_list *new_boundary_head,*new_boundary_tail;
  struct polygon_list *found_plp,*prev_plp,*chk_plp;
  struct polygon *pop,*chk_pop;
  struct node *np,**nodes;
  struct edge *ep;
  struct tet *tp;
  struct tet_list *tlp;
  struct hash_table *hp;
  struct vector3 p1,p2;
  u_int node_index,ref_node_index,ref_node_pos,found_ref_node_pos;
  u_int ni1,ni2,tn1,tn2,itmp;
  u_int cw_target_node,ccw_target_node;
  u_int cw_target_found,ccw_target_found;
  u_int n_mismatch,tet_node_index;
  u_int i,j;
  char ckey[17];

  nodes=tet_mesh->nodes;

  for (i=0; i<17; i++)
  {
    ckey[i]=0x00;
  }

  /* Step 1: build unsorted list of boundary faces */
  /*         associated with each parent boundary node */
  /*         for each boundary face */
  for (plp=tet_mesh->boundary_head; plp!=NULL; plp=plp->next)
  {
    pop=plp->polygon;

    /* for each node in boundary face */
    for (i=0; i<3; i++)
    {
      node_index=pop->node_index[i];
      np=nodes[node_index];

      /* allocate a new polygon_list element to hold this boundary face */
      /* and append it to the boundary face list for this node */
      /* later, we'll sort this list into ccw perimeter order */
      if ((new_plp=(struct polygon_list *)malloc
         (sizeof(struct polygon_list)))==NULL)
      {
        fprintf(stderr,
          "tetgen2FV: out of memory while constructing node boundary face list\n");
        return(1);
      }
      new_plp->polygon=pop;
      new_plp->next=np->boundary_head;
      np->boundary_head=new_plp;
      np->nboundaryfaces++;
    }
  }

  /* Step 2: for each boundary node in tet mesh */
  /*         sort boundary face list into ccw perimeter order, */
  /*         associate boundary faces with their parent edge, */
  /*         and associate boundary faces with their parent tets */
  for (ref_node_index=0; ref_node_index<tet_mesh->nnodes; ref_node_index++)
  {
    np=nodes[ref_node_index];
    if (np->node_type==BOUNDARY || np->node_type==INTERFACE)
    {

      tet_mesh->voronoi_boundary_nedges_tot+=np->nboundaryfaces;

      /* associate parent tet with boundary faces of this node */
      for (plp=np->boundary_head; plp!=NULL; plp=plp->next)
      {
        pop=plp->polygon;

        /* find position of ref_node_index in pop */
        ref_node_pos=0;
        for (i=0; i<3; i++)
        {
          if (pop->node_index[i]==ref_node_index)
          {
            ref_node_pos=i;
          }
        }
        if (ref_node_pos==0)
        {
          ni1 = pop->node_index[1];
          ni2 = pop->node_index[2];
        }
        else if (ref_node_pos==1)
        {
          ni1 = pop->node_index[2];
          ni2 = pop->node_index[0];
        }
        else if (ref_node_pos==2)
        {
          ni1 = pop->node_index[0];
          ni2 = pop->node_index[1];
        }
        /* retrieve one boundary face edge connected to ref_node_index */

        /* sort edge nodes into increasing order */
        tn1 = ref_node_index;
        tn2 = ni1;
        if (tn1>tn2)
        {
          itmp=tn1;
          tn1=tn2;
          tn2=itmp;
        }

        /* construct edge key for this edge */
        sprintf(ckey,"%08x%08x",tn1,tn2);
     
        /* look-up edge key in edge hash table */
        hp=retrieve_key(ckey,tet_mesh->hashmask,tet_mesh->edge_hashtab);
        ep=(struct edge *)hp->contents;

        /* find tets in tet_list of edge sharing this boundary face */
        /* which also share the boundary face */
        for (tlp=ep->tet_head; tlp!=NULL; tlp=tlp->next)
        {
          tp=tlp->tet;

          /* find tet having 3 nodes in common with boundary face */
          /* and only one non-matching node */
          n_mismatch=0;
          for (i=0; i<4; i++)
          {
            tet_node_index=tp->node_index[i];
            if (!(tet_node_index==ref_node_index
                || tet_node_index==ni1
                || tet_node_index==ni2))
            {
              n_mismatch++;
            }
          }
          /* tet found */
          if (n_mismatch==1)
          {

            /* set tet_type */
            tp->tet_type=nodes[pop->node_index[0]]->node_type;

#if 0
            /* mark tet with srfnum of boundary face */
            if (tp->srfnum==-1)
            {
              tp->srfnum=pop->srfnum;
            }
            else if (tp->srfnum!=pop->srfnum)
            {
              fprintf(stderr,
                "tetgen2FV: found tet which spans multiple boundary types\n");
              return(1);
            }
#endif

            /* associate tet with this boundary face */
            /* add tet to shared tet array of this boundary face only once */
            if (pop->shared_tet[0]==NULL)
            {
              pop->shared_tet[0]=tp;
            }
            else if (tp!=pop->shared_tet[0])
            {
              if (tp!=pop->shared_tet[1])
              {
                pop->shared_tet[1]=tp;
              }
            }

#if 0
            /* add boundary face to boundary face list of this tet only once */
            found_pop=0;
            for (plp2=tp->boundary_head;
                 plp2!=NULL && !found_pop;
                 plp2=plp2->next)
            {
              pop2=plp2->polygon;
              if (pop2==pop)
              {
                found_pop=1;
              }
            }
            if (!found_pop)
            {
              if ((new_plp=(struct polygon_list *)malloc
                 (sizeof(struct polygon_list)))==NULL)
              {
                fprintf(stderr,
                  "tetgen2FV: out of memory while constructing node boundary face list\n");
                return(1);
              }
              new_plp->polygon=pop;
              new_plp->next=tp->boundary_head;
              tp->boundary_head=new_plp;
              tp->nboundaryfaces++;
            }
#endif
          }
        }
      }
      

      /* for the first face of boundary face list */
      /* find index of face node that matches ref_node_index */
      /* then find node_index in cw direction */
      /* and call this the "cw_target_node" */
      /* then find node_index in ccw direction */
      /* and call this the "ccw_target_node" */
      pop=np->boundary_head->polygon;

      /* find position of ref_node_index in pop */
      ref_node_pos=0;
      for (i=0; i<3; i++)
      {
        if (pop->node_index[i]==ref_node_index)
        {
          ref_node_pos=i;
        }
      }
      /* set cw and ccw target nodes to be used in sorting the list */
      if (ref_node_pos==0)
      {
        cw_target_node = pop->node_index[1];
        ccw_target_node = pop->node_index[2];
      }
      else if (ref_node_pos==1)
      {
        cw_target_node = pop->node_index[2];
        ccw_target_node = pop->node_index[0];
      }
      else if (ref_node_pos==2)
      {
        cw_target_node = pop->node_index[0];
        ccw_target_node = pop->node_index[1];
      }
      
      /* build a new boundary face list in which to sort the original list */
      /* detach original first face from original list */
      /* and attach the original first face to the top of the new list */
      /* note that the new list now contains just one face */
      /* so the cw (head) and ccw (tail) ends of the list start off */
      /* as coincident */
      new_boundary_head=np->boundary_head;
      np->boundary_head=np->boundary_head->next;
      new_boundary_head->next=NULL;
      new_boundary_tail=new_boundary_head;

      /* now sort remaining faces into new boundary face list */
      for (j=0; j<np->nboundaryfaces-1; j++)
      {
        /* find face that mates with either the cw_target or ccw_target */
        cw_target_found=0;
        ccw_target_found=0;
        found_ref_node_pos=0;
        found_plp=NULL;
        prev_plp=NULL;
        for (chk_plp=np->boundary_head;
             chk_plp!=NULL && (!cw_target_found && !ccw_target_found);
             chk_plp=chk_plp->next)
        {
          chk_pop=chk_plp->polygon;

          /* find position of ref_node_index in chk_pop */
          ref_node_pos=0;
          for (i=0; i<3; i++)
          {
            if (chk_pop->node_index[i]==ref_node_index)
            {
              ref_node_pos=i;
            }
          }

          /* try to match edges of chk_pop with targets */
          if (ref_node_pos==0)
          {
            /* match ccw edge of chk_pop with cw edge of target */
            if (chk_pop->node_index[2]==cw_target_node)
            {
              found_plp=chk_plp;
              cw_target_found=1;
            }
            /* else match cw edge of chk_pop with ccw edge of target */
            else if (chk_pop->node_index[1]==ccw_target_node)
            {
              found_plp=chk_plp;
              ccw_target_found=1;
            }
          }
          else if (ref_node_pos==1)
          {
            /* match ccw edge of chk_pop with cw edge of target */
            if (chk_pop->node_index[0]==cw_target_node)
            {
              found_plp=chk_plp;
              cw_target_found=1;
            }
            /* else match cw edge of chk_pop with ccw edge of target */
            else if (chk_pop->node_index[2]==ccw_target_node)
            {
              found_plp=chk_plp;
              ccw_target_found=1;
            }
          }
          else if (ref_node_pos==2)
          {
            /* match ccw edge of chk_pop with cw edge of target */
            if (chk_pop->node_index[1]==cw_target_node)
            {
              found_plp=chk_plp;
              cw_target_found=1;
            }
            /* else match cw edge of chk_pop with ccw edge of target */
            else if (chk_pop->node_index[0]==ccw_target_node)
            {
              found_plp=chk_plp;
              ccw_target_found=1;
            }
          }
          if (!cw_target_found && !ccw_target_found)
          {
            prev_plp=chk_plp;
          }
          else {
            found_ref_node_pos=ref_node_pos;
          }
        }

        /* detach the found face from the original face list, */
        /* attach it to the head or tail of the new face list, */
        /* and associate the adjacent faces with the boundary edge they share */
        if (cw_target_found)
        {
          /* associate cw face and head face with their shared edge */
          ni1 = ref_node_index;
          ni2 = cw_target_node;

          if (ni1>ni2)
          {
            itmp=ni1;
            ni1=ni2;
            ni2=itmp;
          }

          /* construct edge key for this edge */
          sprintf(ckey,"%08x%08x",ni1,ni2);
     
          /* look-up edge key in edge hash table */
          hp=retrieve_key(ckey,tet_mesh->hashmask,tet_mesh->edge_hashtab);
          ep=(struct edge *)hp->contents;

          if (ep->boundary_face[0]==NULL)
          {
            ep->boundary_face[0]=found_plp->polygon;
            ep->boundary_face[1]=new_boundary_head->polygon;

            ep->edge_type=tet_mesh->nodes[ni1]->node_type;

            tet_mesh->nboundaryedges++;
            if ((ep->cent=(struct vector3 *)malloc
               (sizeof(struct vector3)))==NULL)
            {
              fprintf(stderr,
                "tetgen2FV: out of memory while creating edge circumcenter\n");
              return(1);
            }
            p1.x=tet_mesh->nodes[ni1]->x;
            p1.y=tet_mesh->nodes[ni1]->y;
            p1.z=tet_mesh->nodes[ni1]->z;
            p2.x=tet_mesh->nodes[ni2]->x;
            p2.y=tet_mesh->nodes[ni2]->y;
            p2.z=tet_mesh->nodes[ni2]->z;
            ep->cent->x=0.5*(p1.x+p2.x);
            ep->cent->y=0.5*(p1.y+p2.y);
            ep->cent->z=0.5*(p1.z+p2.z);
          }

          /* attach found face to head of new face list */
          if (found_plp==np->boundary_head)
          {
            np->boundary_head=np->boundary_head->next;
          }
          else {
            prev_plp->next=found_plp->next;
          }
          found_plp->next=new_boundary_head;
          new_boundary_head=found_plp;

          /* now set cw_target_node for the next sorting pass */
          if (found_ref_node_pos==0)
          {
            cw_target_node=found_plp->polygon->node_index[1];
          }
          else if (found_ref_node_pos==1)
          {
            cw_target_node=found_plp->polygon->node_index[2];
          }
          else if (found_ref_node_pos==2)
          {
            cw_target_node=found_plp->polygon->node_index[0];
          }
        }
        else /* ccw_target_found */
        {
          /* associate ccw face and tail face with their shared edge */
          ni1 = ref_node_index;
          ni2 = ccw_target_node;

          if (ni1>ni2)
          {
            itmp=ni1;
            ni1=ni2;
            ni2=itmp;
          }

          /* construct edge key for this edge */
          sprintf(ckey,"%08x%08x",ni1,ni2);
     
          /* look-up edge key in edge hash table */
          hp=retrieve_key(ckey,tet_mesh->hashmask,tet_mesh->edge_hashtab);
          ep=(struct edge *)hp->contents;

          if (ep->boundary_face[0]==NULL)
          {
            ep->boundary_face[0]=found_plp->polygon;
            ep->boundary_face[1]=new_boundary_tail->polygon;

            ep->edge_type=tet_mesh->nodes[ni1]->node_type;

            tet_mesh->nboundaryedges++;
            if ((ep->cent=(struct vector3 *)malloc
               (sizeof(struct vector3)))==NULL)
            {
              fprintf(stderr,
                "tetgen2FV: out of memory while creating edge circumcenter\n");
              return(1);
            }
            p1.x=tet_mesh->nodes[ni1]->x;
            p1.y=tet_mesh->nodes[ni1]->y;
            p1.z=tet_mesh->nodes[ni1]->z;
            p2.x=tet_mesh->nodes[ni2]->x;
            p2.y=tet_mesh->nodes[ni2]->y;
            p2.z=tet_mesh->nodes[ni2]->z;
            ep->cent->x=0.5*(p1.x+p2.x);
            ep->cent->y=0.5*(p1.y+p2.y);
            ep->cent->z=0.5*(p1.z+p2.z);
          }

          /* attach found face to tail of new face list */
          if (found_plp==np->boundary_head)
          {
            np->boundary_head=np->boundary_head->next;
          }
          else {
            prev_plp->next=found_plp->next;
          }
          new_boundary_tail->next=found_plp;
          new_boundary_tail=found_plp;
          new_boundary_tail->next=NULL;

          /* now set ccw_target_node for the next sorting pass */
          if (found_ref_node_pos==0)
          {
            ccw_target_node=found_plp->polygon->node_index[2];
          }
          else if (found_ref_node_pos==1)
          {
            ccw_target_node=found_plp->polygon->node_index[0];
          }
          else if (found_ref_node_pos==2)
          {
            ccw_target_node=found_plp->polygon->node_index[1];
          }
        }
      }

      /* now attach newly sorted list to this node */
      np->boundary_head=new_boundary_head;

    }
  }

  return(0);
}



/* given a 3D Delaunay (tet) mesh
   construct the list of edges that share each Delaunay node
   and the list of tets that share each edge */
int construct_edge_list(struct tet_mesh *tet_mesh)
{
  struct tet_list *tlp,*found_tlp,*prev_tlp;
  struct tet_list *etlp;
  struct tet_list *new_tet_head,*new_tet_tail;
  struct tet *tp,*tmp_tp;
  struct hash_table *hp;
  struct edge *ep;
  struct edge_list *elp;
  double x,y,z,cx,cy,cz,rc,rt;
  u_int i1[6],i2[6];
  u_int ni1,ni2,itmp;
  u_int head_target_node;
  u_int tail_target_node;
  u_int i,j,k;
  byte head_target_found;
  byte tail_target_found;
  char ckey[17];
  char *key;

  tet_mesh->hashsize=0x1000000;
  tet_mesh->hashmask=tet_mesh->hashsize-1;
  if ((tet_mesh->edge_hashtab=init_hashtab(tet_mesh->hashsize))==NULL)
  {
    fprintf(stderr,
      "tetgen2FV: out of memory while initializing edge hash table\n");
    return(1);
  }

  for (i=0; i<17; i++)
  {
    ckey[i]=0x00;
  }

  /* tet node indices are 0, 1, 2, 3 */
  /* tet edges are the 6 sets of node-to-node connections:
     0-1, 1-2, 2-0, 0-3, 1-3, 2-3 */
  /* set up array for scanning node-to-node connections */
  i1[0]=0;
  i2[0]=1;
  
  i1[1]=1;
  i2[1]=2;
  
  i1[2]=2;
  i2[2]=0;
  
  i1[3]=0;
  i2[3]=3;
  
  i1[4]=1;
  i2[4]=3;
  
  i1[5]=2;
  i2[5]=3;

  /* for each tet in the tet_mesh */
  for (tlp=tet_mesh->tet_head; tlp!=NULL; tlp=tlp->next)
  {
    tp=tlp->tet;

    for (i=0; i<4; i++)
    {
      tet_mesh->nodes[tp->node_index[i]]->ntets++;
    }

    /* scan each of the 6 edges of the tet */
    for (i=0; i<6; i++)
    {
      ni1=tp->node_index[i1[i]];
      ni2=tp->node_index[i2[i]];
      tet_mesh->nodes[ni1]->ntets++;
      
      /* sort n1,n2 into increasing order */
      if (ni1>ni2)
      {
        itmp=ni1;
        ni1=ni2;
        ni2=itmp;
      }

      /* construct edge key for this edge */
      sprintf(ckey,"%08x%08x",ni1,ni2);
/*
      printf("%d %d: %08x%08x\n",ni1,ni2,ni1,ni2);
*/
     
      /* look-up edge key in edge hash table */
      hp=retrieve_key(ckey,tet_mesh->hashmask,tet_mesh->edge_hashtab);

      /* if this edge key has not already been inserted into edge hash table */
      if (hp==NULL)
      {
        tet_mesh->nedges++;
        tet_mesh->nodes[ni1]->nedges++;
        tet_mesh->nodes[ni2]->nedges++;

        /* make permanent copy of key */
        if ((key=my_strdup(ckey))==NULL)
        {
          fprintf(stderr,
            "tetgen2FV: out of memory while duplicating edge key\n");
          return(1);
        }

        /* store edge key in edge hash table */
        if ((hp=store_key(key,tet_mesh->hashmask,tet_mesh->edge_hashtab))
          ==NULL)
        {
          fprintf(stderr,
            "tetgen2FV: out of memory while storing key in edge hash table\n");
          return(1);
        }

        /* allocate new edge and insert into contents of edge hash table */
        if ((ep=(struct edge *)malloc(sizeof(struct edge)))==NULL)
        {
          fprintf(stderr,
            "tetgen2FV: out of memory while creating edge\n");
          return(1);
        }
        hp->contents=(void *)ep;
        ep->length=0;
        ep->ntets=0;
        ep->node_index[0]=ni1;
        ep->node_index[1]=ni2;
        ep->voronoi_index=0;
        ep->material_count=0;
        ep->matnum=0;
        ep->edge_type=INTERIOR;
        ep->facet=NULL;
        ep->edge_material=NULL;
        ep->tet_head=NULL;
        ep->boundary_face[0]=NULL;
        ep->boundary_face[1]=NULL;
        ep->cent=NULL;
/*
        if ((tet_mesh->nodes[ni1]->node_type==BOUNDARY
            && tet_mesh->nodes[ni2]->node_type==BOUNDARY)
            || (tet_mesh->nodes[ni1]->node_type==INTERFACE
            && tet_mesh->nodes[ni2]->node_type==INTERFACE))
        {
          ep->edge_type=tet_mesh->nodes[ni1]->node_type;
          tet_mesh->nboundaryedges++;
          if ((ep->cent=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL)
          {
            fprintf(stderr,
              "tetgen2FV: out of memory while creating edge circumcenter\n");
            return(1);
          }
          p1.x=tet_mesh->nodes[ni1]->x;
          p1.y=tet_mesh->nodes[ni1]->y;
          p1.z=tet_mesh->nodes[ni1]->z;
          p2.x=tet_mesh->nodes[ni2]->x;
          p2.y=tet_mesh->nodes[ni2]->y;
          p2.z=tet_mesh->nodes[ni2]->z;
          ep->cent->x=0.5*(p1.x+p2.x);
          ep->cent->y=0.5*(p1.y+p2.y);
          ep->cent->z=0.5*(p1.z+p2.z);
        }
*/

        /* insert new edge into edge_list of the two nodes sharing this edge */
        if ((elp=(struct edge_list *)malloc(sizeof(struct edge_list)))==NULL)
        {
          fprintf(stderr,
            "tetgen2FV: out of memory while creating edge_list\n");
          return(1);
        }
        elp->edge=ep;
        elp->next=tet_mesh->nodes[ni1]->edge_head;
        tet_mesh->nodes[ni1]->edge_head=elp;

        if ((elp=(struct edge_list *)malloc(sizeof(struct edge_list)))==NULL)
        {
          fprintf(stderr,
            "tetgen2FV: out of memory while creating edge_list\n");
          return(1);
        }
        elp->edge=ep;
        elp->next=tet_mesh->nodes[ni2]->edge_head;
        tet_mesh->nodes[ni2]->edge_head=elp;
      }

      /* insert tet into tet_list for this edge */
      /* later we'll sort this list into perimeter order */
      ep=(struct edge *)hp->contents;
      ep->ntets++;

      if ((etlp=(struct tet_list *)malloc(sizeof(struct tet_list)))==NULL)
      {
        fprintf(stderr,
          "tetgen2FV: out of memory while creating tet_list\n");
        return(1);
      }
      etlp->tet=tp;
      etlp->next=ep->tet_head;
      ep->tet_head=etlp;
      ep->matnum=ep->tet_head->tet->matnum;
      

      /* if ni2 is a BOUNDARY or INTERFACE node, create VIRTUAL edge */
      /* by connecting ni2 to itself. The virtual edge will be used to */
      /* hold the Voronoi facet group that divides an INTERFACE Voronoi cell */
      /* or caps off a boundary Voronoi cell */
      
      /* set ni2 back to it's unsorted value */
      /* using ni2 guarantees that we'll visit all 4 nodes of the tet */
      /* at least once (see indexing order above) */
      ni2=tp->node_index[i2[i]];
      
      if ((tet_mesh->nodes[ni2]->node_type==BOUNDARY)
           || (tet_mesh->nodes[ni2]->node_type==INTERFACE))
      {
        /* construct edge key for this edge */
        sprintf(ckey,"%08x%08x",ni2,ni2);
     
        /* look-up edge key in edge hash table */
        hp=retrieve_key(ckey,tet_mesh->hashmask,tet_mesh->edge_hashtab);

        /*if this edge key has not already been inserted into edge hash table*/
        if (hp==NULL)
        {

/* do not include virtual edges in count */
/*
          tet_mesh->nedges++;
          tet_mesh->nodes[ni2]->nedges++;
*/

          /* make permanent copy of key */
          if ((key=my_strdup(ckey))==NULL)
          {
            fprintf(stderr,
              "tetgen2FV: out of memory while duplicating edge key\n");
            return(1);
          }

          /* store edge key in edge hash table */
          if ((hp=store_key(key,tet_mesh->hashmask,tet_mesh->edge_hashtab))
            ==NULL)
          {
            fprintf(stderr,
              "tetgen2FV: out of memory while storing key in edge hash table\n");
            return(1);
          }

          /* allocate new edge and insert into contents of edge hash table */
          if ((ep=(struct edge *)malloc(sizeof(struct edge)))==NULL)
          {
            fprintf(stderr,
              "tetgen2FV: out of memory while creating edge\n");
            return(1);
          }
          hp->contents=(void *)ep;
          ep->length=0;
          ep->ntets=0;
          ep->node_index[0]=ni2;
          ep->node_index[1]=ni2;
          ep->voronoi_index=0;
          ep->material_count=0;
          ep->matnum=0;
          ep->edge_type=VIRTUAL;
          ep->facet=NULL;
          ep->tet_head=NULL;
          ep->boundary_face[0]=NULL;
          ep->boundary_face[1]=NULL;
          ep->cent=NULL;
/*
          tet_mesh->nboundaryedges++;
*/

          /* insert new edge into edge_list of the node sharing this edge */
          if ((elp=(struct edge_list *)malloc(sizeof(struct edge_list)))==NULL)
          {
            fprintf(stderr,
              "tetgen2FV: out of memory while creating edge_list\n");
            return(1);
          }
          elp->edge=ep;
          elp->next=tet_mesh->nodes[ni2]->edge_head;
          tet_mesh->nodes[ni2]->edge_head=elp;
        }
      }
    }
  }

  /* for each edge in edge hash table */
  /* sort edge tet_list into perimeter order */
  for (i=0; i<tet_mesh->hashsize; i++)
  {
    for (hp=tet_mesh->edge_hashtab[i]; hp!=NULL; hp=hp->next)
    {
      ep=(struct edge *)hp->contents;

      if (ep->edge_type!=VIRTUAL)
      {     

        /* for the first tet of edge tet_list */
        /* find first node of tet that is not a node of this edge */
        /* and call this the "head_target_node" */
        /* then find remaining node of tet that is not a node of this edge */
        /* and call this the "tail_target_node" */
        tp=ep->tet_head->tet;
        ni1=ep->node_index[0];
        ni2=ep->node_index[1];

        head_target_found=0;
        tail_target_found=0;
        head_target_node=0;
        tail_target_node=0;
        for (j=0; j<4 && (!head_target_found || !tail_target_found); j++)
        {

          /* look for head_target_node */
          if (!head_target_found
              && tp->node_index[j]!=ni1
              && tp->node_index[j]!=ni2)
          {
            head_target_node=tp->node_index[j];
            head_target_found=1;
          }
  
          /* look for tail_target_node */
          if (head_target_found)
          {
            if (tp->node_index[j]!=ni1
                && tp->node_index[j]!=ni2
                && tp->node_index[j]!=head_target_node)
            {
              tail_target_node=tp->node_index[j];
              tail_target_found=1;
            }
          }
        }
  
        /* build a new tet_list in which to sort the original tet_list */
        /* detach original first tet from original tet_list */
        /* and attach the original first tet to the top of the new tet_list */
        /* note that the new tet_list now contains just one tet */
        /* so the head and tail of the new tet_list start off as coincident */
        new_tet_head=ep->tet_head;
        ep->tet_head=ep->tet_head->next;
        new_tet_head->next=NULL;
        new_tet_tail=new_tet_head;
  
        /* now sort the remaining tets into the new tet_list */
        for (j=0; j<ep->ntets-1; j++)
        {
  
          /* find the first tet that contains either */
          /* the head_target_node or the tail_target_node */
          head_target_found=0;
          tail_target_found=0;
          found_tlp=NULL;
          prev_tlp=NULL;
          for (tlp=ep->tet_head;
               tlp!=NULL && (!head_target_found && !tail_target_found);
               tlp=tlp->next)
          {
            tp=tlp->tet; 
            for (k=0; k<4 && (!head_target_found && !tail_target_found); k++)
            {
              if (tp->node_index[k]==head_target_node)
              {
                found_tlp=tlp;
                head_target_found=1;
              }
              else if (tp->node_index[k]==tail_target_node)
              {
                found_tlp=tlp;
                tail_target_found=1;
              }
            }
            if (!head_target_found && !tail_target_found)
            {
              prev_tlp=tlp;
            }
          }
  
          /* detach the found tet from the original tet_list */
          /* and attach it to the head or tail of the new tet_list as needed */
          if (head_target_found)
          {
            /* attach found tet to the head new tet_list */
            if (found_tlp==ep->tet_head)
            {
              ep->tet_head=ep->tet_head->next;
            }
            else {
              prev_tlp->next=found_tlp->next;
            }
            rc=new_tet_head->tet->r;
            cx=new_tet_head->tet->cent.x;
            cy=new_tet_head->tet->cent.y;
            cz=new_tet_head->tet->cent.z;
            tmp_tp=new_tet_head->tet;
            
            found_tlp->next=new_tet_head;
            new_tet_head=found_tlp;
    
            /* now set the head_target_node for the next sorting pass */
            tp=new_tet_head->tet;
            head_target_found=0;
            for (k=0; k<4 && !head_target_found; k++)
            {
              if (tp->node_index[k]!=ni1
                  && tp->node_index[k]!=ni2
                  && tp->node_index[k]!=head_target_node)
              {
                head_target_node=tp->node_index[k];
                head_target_found=1;
              }
            }
            x=cx-tet_mesh->nodes[head_target_node]->x;
            y=cy-tet_mesh->nodes[head_target_node]->y;
            z=cz-tet_mesh->nodes[head_target_node]->z;
            rt=sqrt(x*x + y*y + z*z);
  /*
            if (rt < rc)
  */
            if (rc-rt > 1e-15)
            {
              tp->delaunay=0;
              fprintf(stderr,"Delaunay circumcenter violation at tets: %d %d; r: %.17g %.17g; edge: %d %d\n",new_tet_head->tet->tet_index,tmp_tp->tet_index,rc,rt,ni1,ni2);
            }
          }
          else /* tail_target_found */
          {
            /* attach found tet to the tail new tet_list */
            if (found_tlp==ep->tet_head)
            {
              ep->tet_head=ep->tet_head->next;
            }
            else {
              prev_tlp->next=found_tlp->next;
            }
            rc=new_tet_tail->tet->r;
            cx=new_tet_tail->tet->cent.x;
            cy=new_tet_tail->tet->cent.y;
            cz=new_tet_tail->tet->cent.z;
            tmp_tp=new_tet_tail->tet;
    
            new_tet_tail->next=found_tlp;
            new_tet_tail=found_tlp;
            new_tet_tail->next=NULL;
  
            /* now set the tail_target_node for the next sorting pass */
            tp=new_tet_tail->tet;
            tail_target_found=0;
            for (k=0; k<4 && !tail_target_found; k++)
            {
              if (tp->node_index[k]!=ni1
                  && tp->node_index[k]!=ni2
                  && tp->node_index[k]!=tail_target_node)
              {
                tail_target_node=tp->node_index[k];
                tail_target_found=1;
              }
            }
            x=cx-tet_mesh->nodes[tail_target_node]->x;
            y=cy-tet_mesh->nodes[tail_target_node]->y;
            z=cz-tet_mesh->nodes[tail_target_node]->z;
            rt=sqrt(x*x + y*y + z*z);
  /*
            if (rt < rc)
  */
            if (rc-rt > 1e-15)
            {
              tp->delaunay=0;
              fprintf(stderr,"Delaunay circumcenter violation at tets: %d %d; r: %g %g; edge: %d %d\n",new_tet_tail->tet->tet_index,tmp_tp->tet_index,rc,rt,ni1,ni2);
            }
          }
  
        }
  
        /* attach sorted tet_list to this edge */
        ep->tet_head=new_tet_head;
      }

    }
  }

  return(0);
}



int construct_material_info(struct tet_mesh *tet_mesh)
{
struct tet_list *tlp,*new_tlp;
struct tet *tp;
struct node **nodes;
struct node *np,*np2;
struct edge_list *elp,*new_elp;
struct edge *ep;
struct edge_material_info *emip;
struct voronoi_facet *vfp;
struct polygon *pop;
u_int n_voronoi_facet_nodes;
u_int ni;
u_int i,j,k,m;
byte not_found,prev_tet_type,boundary_tet,boundary_tet_next;

  nodes=tet_mesh->nodes;
  for (i=0; i<tet_mesh->nnodes; i++)
  {
    np=nodes[i];
    if (np->node_type==INTERFACE)
    {
      for (elp=np->edge_head; elp!=NULL; elp=elp->next)
      {
        ep=elp->edge;
  
        /* copy interior edge into edge_list of matching node_material_info */
        /* interior edges are one material so don't need edge_material_info */
        if (ep->edge_type==INTERIOR)
        {
          /* find partner node sharing this edge */
          ni=ep->node_index[0];
          if (i==ni)
          {
            ni=ep->node_index[1];
          }
          np2=nodes[ni];


          /* determine material of INTERIOR edge */
          /* edge material of an interior edge is always
             the same as its shared tets */
          ep->matnum=ep->tet_head->tet->matnum;

          if (np2->node_type==INTERIOR)
          {
            /* increment interface neighbor count of neighbor node */
            np2->ninterfaceneighbors++;
          }

          /* determine material of INTERIOR edge */
          /* silly broken method: */
#if 0
          if (np2->node_type==INTERIOR)
          {
            /* set edge material to same type as partner node material */
            /* and increment interface neighbor count */
            ep->matnum=np2->material_head->matnum;
            np2->ninterfaceneighbors++;
          }
          else if (np2->node_type==INTERFACE)
          {
            /* find child node with matching material */
            /* set edge material to same type as matching child node */
            for (j=0; j<np->material_count; j++)
            {
              for (m=0; m<np2->material_count; m++)
              {
                if (np->node_material[j].matnum==np2->node_material[m].matnum)
                {
                  ep->matnum=np2->node_material[m].matnum;
                  np2->ninterfaceneighbors++;
                }
              }
            }
          }
#endif /* 0 */
          

          for (j=0; j<np->material_count; j++)
          {
            if (np->node_material[j].matnum==ep->matnum)
            {
              np->node_material[j].nedges++;
              np->node_material[j].ntets=ep->ntets;
              np->node_material[j].nfacets++;
              if ((new_elp=(struct edge_list *)malloc
                 (sizeof(struct edge_list)))==NULL)
              {
                fprintf(stderr,
                "tetgen2FV: out of memory while constructing material info\n");
                return(1);
              }
              new_elp->edge=ep;
              new_elp->next=np->node_material[j].edge_head;
              np->node_material[j].edge_head=new_elp;
            }
          }
        }
        /* copy virtual edge into edge_list of each node_material_info */
        /* virtual edges are special so edge_material_info is not needed */
        else if (ep->edge_type==VIRTUAL)
        {
          for (j=0; j<np->material_count; j++)
          {
/* do not include virtual edges in the count */
/*
            np->node_material[j].nedges++;
*/
            np->node_material[j].ntets=ep->ntets;
            np->node_material[j].nfacets++;
            if ((new_elp=(struct edge_list *)malloc
               (sizeof(struct edge_list)))==NULL)
            {
              fprintf(stderr,
              "tetgen2FV: out of memory while constructing material info\n");
              return(1);
            }
            new_elp->edge=ep;
            new_elp->next=np->node_material[j].edge_head;
            np->node_material[j].edge_head=new_elp;
          }
        }
        /* copy interface edge into edge_list of matching node_material_info */
        /* interface edges are multi-material so we need edge_material_info */
        else if (ep->edge_type==INTERFACE)
        {
          if (ep->edge_material==NULL)
          {
            /* allocate memory for edge_material_info */
            ep->material_count=np->material_count;
            if ((emip=(struct edge_material_info *)malloc
               (ep->material_count*sizeof(struct edge_material_info)))==NULL)
            {
              fprintf(stderr,
              "tetgen2FV: out of memory while constructing material info\n");
              return(1);
            }
            ep->edge_material=emip;
            for (j=0; j<ep->material_count; j++)
            {
              emip[j].matnum=np->node_material[j].matnum;
              emip[j].ntets=0;
              emip[j].tet_head=NULL;
              emip[j].facet=NULL;
            }

            /* build material-specific tet_list and Voronoi facet
               for each edge material */
            for (j=0; j<ep->material_count; j++)
            {
              /* build material-specific tet_list for each edge material */
              /* new tet_list will already be in perimeter order */
              for (tlp=ep->tet_head; tlp!=NULL; tlp=tlp->next)
              {
                tp=tlp->tet;
                if (tp->matnum==ep->edge_material[j].matnum)
                {
                  ep->edge_material[j].ntets++;
                  if ((new_tlp=(struct tet_list *)malloc
                     (sizeof(struct tet_list)))==NULL)
                  {
                    fprintf(stderr,
                    "tetgen2FV: out of memory while constructing material info\n");
                    return(1);
                  }
                  new_tlp->tet=tp;
                  new_tlp->next=ep->edge_material[j].tet_head;
                  ep->edge_material[j].tet_head=new_tlp;
                }
              }

              /* now build material-specific Voronoi facet for this edge */
              n_voronoi_facet_nodes=ep->edge_material[j].ntets+3;
              /* make the one Voronoi facet associated with this edge --
                 the facet is a planar-polygon whose vertices are the
                 tet circumcenter points associated with this edge */
    
              /* allocate memory for the Voronoi facet
                 associated with this edge */
              if ((vfp=(struct voronoi_facet *)malloc
                 (sizeof(struct voronoi_facet)))==NULL)
              {
                fprintf(stderr,
                 "tetgen2FV: out of memory while constructing material info\n");
                return(1);
              }
              vfp->area=0;
              vfp->nnodes=n_voronoi_facet_nodes;
              vfp->node_index=NULL;
              vfp->next=NULL;
              ep->edge_material[j].facet=vfp;
  
              if ((vfp->node_index=(u_int *)malloc
                 (vfp->nnodes*sizeof(u_int)))==NULL)
              {
                fprintf(stderr,
                 "tetgen2FV: out of memory while constructing material info\n");
                return(1);
              }

              if (ep->edge_material[j].ntets==1)
              {
                /* circumcenter points form partial perimeter around edge --
                   close the perimeter by using additional points lying on the
                   boundary faces sharing this edge and the edge mid-point */
                k=0;
                vfp->node_index[k++]=ep->boundary_face[0]->voronoi_index;
                vfp->node_index[k++]=ep->edge_material[j].tet_head->tet->tet_index;
                vfp->node_index[k++]=ep->boundary_face[1]->voronoi_index;
                vfp->node_index[k++]=ep->voronoi_index;
              }
              else if (ep->edge_material[j].ntets==2)
              {
                /* circumcenter points form partial perimeter around edge --
                   close the perimeter by using additional points lying on the
                   boundary faces sharing this edge and the edge mid-point */
                k=0;
                for (tlp=ep->edge_material[j].tet_head; tlp!=NULL; tlp=tlp->next)
                {
                  tp=tlp->tet;
  
                  if (tlp==ep->edge_material[j].tet_head)
                  {
                    /* include tet circumcenter before additional points */
                    vfp->node_index[k++]=tp->tet_index;
  
                    /* determine perimeter order of boundary faces and
                       include additional points lying on the boundary faces
                       sharing this edge and the edge mid_point */
                    /* find boundary face of this edge
                       which is also shared by this tet */
                    not_found=1;
                    for (m=0; m<2 && not_found; m++)
                    {
                      pop=ep->boundary_face[m];
                      if (pop->shared_tet[0]==tp || pop->shared_tet[1]==tp)
                      {
                        not_found=0;
                        if (m==0)
                        {
                          vfp->node_index[k++]=ep->boundary_face[0]->voronoi_index;
                          vfp->node_index[k++]=ep->voronoi_index;
                          vfp->node_index[k++]=ep->boundary_face[1]->voronoi_index;
                        }
                        else
                        {
                          vfp->node_index[k++]=ep->boundary_face[1]->voronoi_index;
                          vfp->node_index[k++]=ep->voronoi_index;
                          vfp->node_index[k++]=ep->boundary_face[0]->voronoi_index;
                        }
                      }
                    }
                  }
                  else
                  {
                    /* include tet circumcenter */
                    vfp->node_index[k++]=tp->tet_index;
                  }
                }
              }
              else if (ep->edge_material[j].ntets > 2)
              {
  
                /* circumcenter points form partial perimeter around edge --
                   close the perimeter by using additional points lying on the
                   boundary faces sharing this edge and the edge mid-point */
                k=0;
                prev_tet_type=INTERIOR;
                for (tlp=ep->edge_material[j].tet_head; tlp!=NULL; tlp=tlp->next)
                {
                  tp=tlp->tet;
  
                  boundary_tet=0;
                  if (tp->tet_type==INTERFACE)
                  {
                    /* boundary status only valid
                       if tet and edge share boundary face */
                    /* find boundary face of this edge
                       which is also shared by this tet, if there is one */
                    not_found=1;
                    for (m=0; m<2 && not_found; m++)
                    {
                      pop=ep->boundary_face[m];
                      if (pop->shared_tet[0]==tp || pop->shared_tet[1]==tp)
                      {
                        not_found=0;
                        boundary_tet=1;
                      }
                    }
                  }
  
                  if (!boundary_tet)
                  {
                    /* include tet circumcenter */
                    vfp->node_index[k++]=tp->tet_index;
                  }
                  if (boundary_tet && prev_tet_type==INTERFACE)
                  {
                    /* include tet circumcenter */
                    vfp->node_index[k++]=tp->tet_index;
                  }
                  if (boundary_tet && tlp==ep->edge_material[j].tet_head)
                  {
                    /* include tet circumcenter before additional points */
                    vfp->node_index[k++]=tp->tet_index;
  
                    boundary_tet_next=0;
                    if (tlp->next->tet->tet_type==INTERFACE)
                    {
                      /* boundary status only valid
                         if next tet and edge share boundary face */
                      /* find boundary face of this edge
                         which is also shared by next tet, if there is one */
                      not_found=1;
                      for (m=0; m<2 && not_found; m++)
                      {
                        pop=ep->boundary_face[m];
                        if (pop->shared_tet[0]==tlp->next->tet
                            || pop->shared_tet[1]==tlp->next->tet)
                        {
                          not_found=0;
                          boundary_tet_next=1;
                        }
                      }
                    }
  
                    if (boundary_tet_next)
                    {
                      /* determine perimeter order of boundary faces and
                         include additional points lying on the boundary faces
                         sharing this edge and the edge mid_point */
                      /* find boundary face of this edge
                         which is also shared by this tet */
                      not_found=1;
                      for (m=0; m<2 && not_found; m++)
                      {
                        pop=ep->boundary_face[m];
                        if (pop->shared_tet[0]==tp || pop->shared_tet[1]==tp)
                        {
                          not_found=0;
                          if (m==0)
                          {
                            vfp->node_index[k++]=ep->boundary_face[0]->voronoi_index;
                            vfp->node_index[k++]=ep->voronoi_index;
                            vfp->node_index[k++]=ep->boundary_face[1]->voronoi_index;
                          }
                          else
                          {
                            vfp->node_index[k++]=ep->boundary_face[1]->voronoi_index;
                            vfp->node_index[k++]=ep->voronoi_index;
                            vfp->node_index[k++]=ep->boundary_face[0]->voronoi_index;
                          }
                        }
                      }
                    }
                  }
                  if (boundary_tet
                      && tlp!=ep->edge_material[j].tet_head
                      && prev_tet_type!=INTERFACE)
                  {
                    /* include tet circumcenter before additional points */
                    vfp->node_index[k++]=tp->tet_index;
  
                    /* determine perimeter order of boundary faces and
                       include additional points lying on the boundary faces
                       sharing this edge and the edge mid_point */
                    /* find boundary face of this edge
                       which is also shared by this tet */
                    not_found=1;
                    for (m=0; m<2 && not_found; m++)
                    {
                      pop=ep->boundary_face[m];
                      if (pop->shared_tet[0]==tp || pop->shared_tet[1]==tp)
                      {
                        not_found=0;
                        if (m==0)
                        {
                          vfp->node_index[k++]=ep->boundary_face[0]->voronoi_index;
                          vfp->node_index[k++]=ep->voronoi_index;
                          vfp->node_index[k++]=ep->boundary_face[1]->voronoi_index;
                        }
                        else
                        {
                          vfp->node_index[k++]=ep->boundary_face[1]->voronoi_index;
                          vfp->node_index[k++]=ep->voronoi_index;
                          vfp->node_index[k++]=ep->boundary_face[0]->voronoi_index;
                        }
                      }
                    }
                  }
                  if (boundary_tet)
                  {
                    prev_tet_type=INTERFACE;
                  }
                  else
                  {
                    prev_tet_type=INTERIOR;
                  }
  
                }
  
              }

              /* calculate area of Voronoi facet (dA) */
	      vfp->area=polyarea(tet_mesh,vfp);
              if (vfp->area==0)
              {
                fprintf(stderr,"warning: zero Voronoi facet area at node: %d %g\n",i,vfp->area);
              }
            }
          }

          /* find matching edge and node materials */
          for (j=0; j<ep->material_count; j++)
          {
            for (k=0; k<np->material_count; k++)
            {
              if (ep->edge_material[j].matnum==np->node_material[k].matnum)
              {
                np->node_material[k].nedges++;
                np->node_material[k].ntets=ep->edge_material[j].ntets;
                np->node_material[k].nfacets++;
                if ((new_elp=(struct edge_list *)malloc
                   (sizeof(struct edge_list)))==NULL)
                {
                  fprintf(stderr,
                  "tetgen2FV: out of memory while constructing material info\n");
                  return(1);
                }
                new_elp->edge=ep;
                new_elp->next=np->node_material[k].edge_head;
                np->node_material[k].edge_head=new_elp;
              }
            }
          }
        }

      }
      /* calculate volume of each material-specific Voronoi cell
         associated with this interface node */
      for (j=0;j<np->material_count;j++)
      {
        np->node_material[j].volume=polyvol_by_material_1(tet_mesh,np,j);
        if (np->node_material[j].volume<=0)
        {
          fprintf(stderr,"warning: negative or zero volume for material specific Voronoi element at node: %d %g\n",np->node_index,np->node_material[j].volume);
        }
      }
    }
  }

  return(0);
}



/*
  Output shape of the Voronoi element associated with a given node
  Output will be in Hughes Hoppe mesh format:

  Vertex i x y z
  ...
  Face i v1 v2 v3
  ...

*/
int output_voronoi_element(struct tet_mesh *tet_mesh, struct node *node)
{
  struct node **nodes;
  struct edge *ep;
  struct edge_list *elp;
  struct voronoi_facet *vfp;
  struct vector3 p1,p2,pint,v1,v2,v_node,v_norm;
  struct vector3 **vnp;
  u_int i;
  u_int n1;
  u_int n2;
  u_int n3;
  u_int v,nnodes;
  double dot;
  byte done;

  nodes=tet_mesh->nodes;
  vnp=tet_mesh->voronoi_nodes;
 
  printf("# Output Voronoi element for node %u %.17g %.17g %.17g\n",node->node_index,node->x,node->y,node->z);

  for (i=0;i<tet_mesh->nvoronoinodes;i++)
  {
    printf("Vertex %u %.17g %.17g %.17g\n",i+1,vnp[i]->x,vnp[i]->y,vnp[i]->z);
  }

  /* first find a point on the interior of the Voronoi cell */
  /* this is done by finding an interior Delaunay edge  */
  /* and placing a point 1/4 of the way along this edge away from node */
  done=0;
  for (elp=node->edge_head; elp!=NULL && !done; elp=elp->next)
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


  i=0;
  for (elp=node->edge_head; elp!=NULL; elp=elp->next)
  {
    for (vfp=elp->edge->facet; vfp!=NULL; vfp=vfp->next)
    {

    /* orient and triangulate each face.*/

      n3 = vfp->node_index[vfp->nnodes-1]; 

    /* determine orientation of Voronoi facet */
      /* construct vector from node to first vertex of facet */
      p1.x = vnp[vfp->node_index[0]]->x; 
      p1.y = vnp[vfp->node_index[0]]->y; 
      p1.z = vnp[vfp->node_index[0]]->z; 
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

          i++;

          n1 = vfp->node_index[v];
          n2 = vfp->node_index[v+1];

          printf("Face %u %u %u %u\n",i,n1+1,n2+1,n3+1);

        }
      }
      else
      {
        /* facet normal points inward -- need to reverse it */
        n3 = vfp->node_index[1]; 
        for (v=0; v<nnodes-2; v++)
        {

          i++;

          if (v==0)
          {
            n1 = vfp->node_index[0];
          }
          else
          {
            n1 = vfp->node_index[nnodes-v];
          }
          n2 = vfp->node_index[nnodes-v-1];

          printf("Face %u %u %u %u\n",i,n1+1,n2+1,n3+1);
      
        }
      }
    }
  }

  return(0);
}



/*
  Output shape of the material-specific Voronoi element
  associated with a given node.
  Output will be in Hughes Hoppe mesh format:

  Vertex i x y z
  ...
  Face i v1 v2 v3
  ...

*/
int output_voronoi_element_by_material(struct tet_mesh *tet_mesh,
                                          struct node *node,
                                          u_int node_matindex)
{
  struct node **nodes;
  struct edge_list *elp;
  struct edge *ep;
  struct voronoi_facet *facet_head,*vfp;
  struct vector3 p1,p2,pint,v1,v2,v_node,v_norm;
  struct vector3 **vnp;
  double dot;
  u_int n1;
  u_int n2;
  u_int n3;
  u_int v,nnodes,node_count;
  u_int edge_matindex;
  u_int i,j;
  byte found;

  nodes=tet_mesh->nodes;
  vnp=tet_mesh->voronoi_nodes;

  printf("# Output material-specific Voronoi element for node %u (%d) %.17g %.17g %.17g\n",node->node_index,node_matindex,node->x,node->y,node->z);

  for (i=0;i<tet_mesh->nvoronoinodes;i++)
  {
    printf("Vertex %u %.17g %.17g %.17g\n",i+1,vnp[i]->x,vnp[i]->y,vnp[i]->z);
  }

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

  /* now output material-specific Voronoi cell */
  i=0;
  for (elp=node->node_material[node_matindex].edge_head; elp!=NULL; elp=elp->next)
  {
    ep=elp->edge;

    if (ep->material_count>1)
    {
      edge_matindex=0;
      for (j=0;j<ep->material_count;j++)
      {
        if (node->node_material[node_matindex].matnum==ep->edge_material[j].matnum)
        {
          edge_matindex=j;
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
          i++;

          n1 = vfp->node_index[v];
          n2 = vfp->node_index[v+1];

          printf("Face %u %u %u %u\n",i,n1+1,n2+1,n3+1);
        }
      }
      else
      {
        /* facet normal points inward -- need to reverse it */
        n3 = vfp->node_index[1]; 
        for (v=0; v<nnodes-2; v++)
        {
          i++;

          if (v==0)
          {
            n1 = vfp->node_index[0];
          }
          else
          {
            n1 = vfp->node_index[nnodes-v];
          }
          n2 = vfp->node_index[nnodes-v-1];
      
          printf("Face %u %u %u %u\n",i,n1+1,n2+1,n3+1);
        }
      }
    }
  }

  return(0); 
}


 
/*
  Output shape of the Voronoi element associated with a given node
  Output will be in dx format.
*/
int output_voronoi_element_dx(struct tet_mesh *tet_mesh,
                              struct node *node,
                              char *outfile_name)
{
  FILE *outfile;
  struct node **nodes;
  struct edge *ep;
  struct edge_list *elp;
  struct voronoi_facet *vfp;
  struct vector3 p1,p2,pint,v1,v2,v_node,v_norm;
  struct vector3 **vnp,*vp;
  double dot;
  double l1,l2,edge_ratio;
  float nx,ny,nz,fi;
  u_int i;
  u_int n1;
  u_int n2;
  u_int n3;
  u_int v,nnodes;
  u_int voronoi_node_count;
  u_int voronoi_triangle_count;
  u_int data_offset;
  byte found;
  char obj_name[64];

  nodes=tet_mesh->nodes;
  vnp=tet_mesh->voronoi_nodes;

  voronoi_node_count=tet_mesh->nvoronoinodes;

  voronoi_triangle_count=0;
  for (elp=node->edge_head; elp!=NULL; elp=elp->next)
  {
    for (vfp=elp->edge->facet; vfp!=NULL; vfp=vfp->next)
    {
      nnodes=vfp->nnodes;
      voronoi_triangle_count+=nnodes-2;
    }
  }

  if ((outfile=fopen(outfile_name,"w"))==NULL)
  {
    fprintf(stderr,"Cannot open DX outfile %s",outfile_name);
    return(1);
  }

  data_offset=0;
 
  fprintf(outfile,"# Output Voronoi element for node %u %.17g %.17g %.17g\n",node->node_index,node->x,node->y,node->z);

  sprintf(obj_name,"voronoi mesh nodes");
  fprintf(outfile,"# %s\n",obj_name);
  fprintf(outfile,"object \"%s\" class array type float rank 1 shape 3 items %d lsb ieee data %d\n",obj_name,voronoi_node_count,data_offset);
  fprintf(outfile,"attribute \"dep\" string \"positions\"\n");
  fprintf(outfile,"#\n");
  data_offset+=voronoi_node_count*3*4;

  sprintf(obj_name,"voronoi mesh facets");
  fprintf(outfile,"# %s\n",obj_name);
  fprintf(outfile,"object \"%s\" class array type int rank 1 shape 3 items %d lsb ieee data %d\n",obj_name,voronoi_triangle_count,data_offset);
  fprintf(outfile,"attribute \"dep\" string \"connections\"\n");
  fprintf(outfile,"attribute \"ref\" string \"positions\"\n");
  fprintf(outfile,"attribute \"element type\" string \"triangles\"\n");
  fprintf(outfile,"#\n");
  data_offset+=voronoi_triangle_count*3*4;

  sprintf(obj_name,"voronoi mesh triangle facet indices");
  fprintf(outfile,"# %s\n",obj_name);
  fprintf(outfile,"object \"%s\" class array type float rank 0 items %d lsb ieee data %d\n",obj_name,voronoi_triangle_count,data_offset);
  fprintf(outfile,"attribute \"dep\" string \"connections\"\n");
  fprintf(outfile,"#\n");
  data_offset+=voronoi_triangle_count*1*4;

/*
  sprintf(obj_name,"voronoi mesh edges");
  fprintf(outfile,"# %s\n",obj_name);
  fprintf(outfile,"object \"%s\" class array type int rank 1 shape 2 items %d lsb ieee data %d\n",obj_name,voronoi_edge_count,data_offset);
  fprintf(outfile,"attribute \"dep\" string \"connections\"\n");
  fprintf(outfile,"attribute \"ref\" string \"positions\"\n");
  fprintf(outfile,"attribute \"element type\" string \"lines\"\n");
  fprintf(outfile,"#\n");
  data_offset+=voronoi_edge_count*2*4;

  sprintf(obj_name,"voronoi mesh edge facet indices");
  fprintf(outfile,"# %s\n",obj_name);
  fprintf(outfile,"object \"%s\" class array type float rank 0 items %d lsb ieee data %d\n",obj_name,voronoi_edge_count,data_offset);
  fprintf(outfile,"attribute \"dep\" string \"connections\"\n");
  fprintf(outfile,"#\n");
  data_offset+=voronoi_edge_count*1*4;
*/

  sprintf(obj_name,"voronoi mesh facets field");
  fprintf(outfile,"# %s\n",obj_name);
  fprintf(outfile,"object \"%s\" class field\n",obj_name);
  fprintf(outfile,"component \"positions\" value \"%s\"\n",
    "voronoi mesh nodes");
  fprintf(outfile,"component \"connections\" value \"%s\"\n",
    "voronoi mesh facets");
  fprintf(outfile,"component \"facet index\" value \"%s\"\n",
    "voronoi mesh triangle facet indices");
  fprintf(outfile,"#\n");

  fprintf(outfile,"object \"default\" class group\n");
  fprintf(outfile,"member \"voronoi_facets\" value \"%s\"\n",
    "voronoi mesh facets field");
  fprintf(outfile,"#\n");
  fprintf(outfile,"end\n");

  /* output voronoi mesh nodes */
  for (i=0; i<voronoi_node_count; i++)
  {
    vp=vnp[i];
    nx=vp->x;
    ny=vp->y;
    nz=vp->z;
    fwrite(&nx,sizeof(float),1,outfile);
    fwrite(&ny,sizeof(float),1,outfile);
    fwrite(&nz,sizeof(float),1,outfile);
  }

  /* first find a point on the interior of the Voronoi cell */
  /* this is done by finding an interior Delaunay edge  */
  /* and placing a point 1/4 of the way along this edge away from node */
/*
  found=0;
  for (elp=node->edge_head; elp!=NULL && !found; elp=elp->next)
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
*/
  pint.x=node->x;
  pint.y=node->y;
  pint.z=node->z;


  /* output voronoi mesh facets as triangles */
  i=0;
  for (elp=node->edge_head; elp!=NULL; elp=elp->next)
  {
    for (vfp=elp->edge->facet; vfp!=NULL; vfp=vfp->next)
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
          if (edge_ratio>1e-6)
          {
            found=1;
          }
        }

    /* determine orientation of Voronoi facet */
      /* construct vector from interior to first vertex of facet */
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

      nnodes=vfp->nnodes;
      if (dot>=0)
      {
        /* facet normal already points outward */
        n3 = vfp->node_index[nnodes-1]; 
        for (v=0; v<nnodes-2; v++)
        {

          i++;

          n1 = vfp->node_index[v];
          n2 = vfp->node_index[v+1];

          fwrite(&n1,sizeof(u_int),1,outfile);
          fwrite(&n2,sizeof(u_int),1,outfile);
          fwrite(&n3,sizeof(u_int),1,outfile);

        }
      }
      else
      {
        /* facet normal points inward -- need to reverse it */
        n3 = vfp->node_index[1]; 
        for (v=0; v<nnodes-2; v++)
        {

          i++;

          if (v==0)
          {
            n1 = vfp->node_index[0];
          }
          else
          {
            n1 = vfp->node_index[nnodes-v];
          }
          n2 = vfp->node_index[nnodes-v-1];

          fwrite(&n1,sizeof(u_int),1,outfile);
          fwrite(&n2,sizeof(u_int),1,outfile);
          fwrite(&n3,sizeof(u_int),1,outfile);
      
        }
      }
    }
  }

  /* output voronoi mesh triangle facet indices */
  i=0;
  for (elp=node->edge_head; elp!=NULL; elp=elp->next)
  {
    for (vfp=elp->edge->facet; vfp!=NULL; vfp=vfp->next)
    {
      fi=i;
      nnodes=vfp->nnodes;
      for (v=0; v<nnodes-2; v++)
      {
        fwrite(&fi,sizeof(float),1,outfile);
      }
    }
    i++;
  }

  return(0);
}



/*
  Check that tet mesh edges are truly normal to Voronoi faces

*/
int check_normals(struct tet_mesh *tet_mesh)
{

  struct node **nodes,*np;
  struct edge_list *elp;
  struct edge *ep;
  struct voronoi_facet *vfp;
  struct vector3 p1,p2,v1,v2,v_edge,v_norm;
  struct vector3 **vnp;
  u_int i;
  u_int n3;
  double dot,check,tol;

  vnp=tet_mesh->voronoi_nodes;
  nodes=tet_mesh->nodes;

  for (i=0; i<tet_mesh->nnodes; i++)
  {
  
    np=nodes[i];
    for (elp=np->edge_head; elp!=NULL; elp=elp->next)
    {
      ep=elp->edge;

      for (vfp=ep->facet; ep->edge_type !=VIRTUAL && vfp!=NULL; vfp=vfp->next)
      {

        /* construct normal vector of facet */
        p1.x = vnp[vfp->node_index[0]]->x; 
        p1.y = vnp[vfp->node_index[0]]->y; 
        p1.z = vnp[vfp->node_index[0]]->z; 
        p2.x = vnp[vfp->node_index[1]]->x; 
        p2.y = vnp[vfp->node_index[1]]->y; 
        p2.z = vnp[vfp->node_index[1]]->z; 
        vectorize(&p1,&p2,&v1);
        n3 = vfp->node_index[vfp->nnodes-1]; 
        p2.x = vnp[n3]->x; 
        p2.y = vnp[n3]->y; 
        p2.z = vnp[n3]->z; 
        vectorize(&p1,&p2,&v2);
        cross_prod(&v1,&v2,&v_norm);
        normalize(&v_norm);

        
        /* construct edge vector of facet */
        p1.x = nodes[ep->node_index[0]]->x; 
        p1.y = nodes[ep->node_index[0]]->y; 
        p1.z = nodes[ep->node_index[0]]->z; 
        p2.x = nodes[ep->node_index[1]]->x; 
        p2.y = nodes[ep->node_index[1]]->y; 
        p2.z = nodes[ep->node_index[1]]->z; 
        vectorize(&p1,&p2,&v_edge);
        normalize(&v_edge);

        dot=fabs(dot_prod(&v_edge,&v_norm));
        check = fabs(1.0-dot);
        tol=1e-16;
        
        if (check>tol)
        {
          printf("Tet edge not normal to Voronoi facet: %u %u %.17g\n",
            ep->node_index[0],ep->node_index[1],dot);
        }

      }
    }
  }

  return(0);
}
	

