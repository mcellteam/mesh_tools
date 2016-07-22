#include <stdio.h> 
#include <string.h> 
#include <sys/types.h>
#include <math.h>
#include "strfunc.h"
#include "tetgen2FV.h"
#include "voronoi.h"
#include "output_edsim.h"

#ifdef DEBUG
#define no_printf printf
#endif

/* output volume mesh data for use in electrodiffusion simulation */
int output_edsim(struct tet_mesh *tet_mesh, char *outfile_name)
{

  FILE *outfile;
  struct node **nodes;
  struct node *np,*np2;
  struct edge_list *elp;
  struct edge *ep;
  u_int node_count,total_node_count;
  u_int node_index,parent_index,child_index1,child_index2;
  u_int ni;
  u_int tot_coupling_count,coupling_count;
  u_int edge_matindex;
  u_int i,j,k;


  fprintf(stderr,"tetgen2FV: output edsim data:\n");

  node_count=tet_mesh->nnodes;
  nodes=tet_mesh->nodes;
  tot_coupling_count=0;

  if ((outfile=fopen(outfile_name,"w"))==NULL)
  {
    fprintf(stderr,"Cannot open edsim outfile %s",outfile_name);
    return(1);
  }

  /* output part 1, node list */

  /* calculate total node count as node_count + 2*number of interface nodes */
  total_node_count=node_count;
  for (i=0; i<node_count; i++)
  {
    np=nodes[i];
    if (np->node_type==INTERFACE) /* output interface node as a parent node */
    {
      total_node_count+=2;
    }
  }

  /* output total_node_count */
  fprintf(outfile,"%u\n",total_node_count);
 
  /* first output parent(interface), interior, and world_boundary nodes */
  for (i=0; i<node_count; i++)
  {
    np=nodes[i];

    node_index=i+1;
    if (np->node_type==INTERFACE) /* output interface node as a parent node */
    {
      parent_index=node_index;
      child_index1=np->node_material[0].child_index+1;
      child_index2=np->node_material[1].child_index+1;
      fprintf(outfile,"%u %.17g %.17g %.17g %.17g %u %u %u %u %d %u\n",
        node_index,np->x,np->y,np->z,np->volume_element->volume,PARENT,
        parent_index,child_index1,child_index2,np->material_head->matnum,
        np->nedges);
      tot_coupling_count+=np->nedges;
/*
      printf(" %u %u\n",node_index,np->nedges);
*/
    }
    else if (np->node_type==BOUNDARY) /* output world_boundary node */
    {
      parent_index=node_index;
      child_index1=node_index;
      child_index2=node_index;
      fprintf(outfile,"%u %.17g %.17g %.17g %.17g %u %u %u %u %d %u\n",
        node_index,np->x,np->y,np->z,np->volume_element->volume,np->node_type,
        parent_index,child_index1,child_index2,np->material_head->matnum,
        np->nedges);
      tot_coupling_count+=np->nedges;
/*
      printf(" %u %u\n",node_index,np->nedges);
*/
    }
    else if (np->node_type==INTERIOR) /* output interior node */
    {
      parent_index=node_index;
      child_index1=node_index;
      child_index2=node_index;
      fprintf(outfile,"%u %.17g %.17g %.17g %.17g %u %u %u %u %d %u\n",
        node_index,np->x,np->y,np->z,np->volume_element->volume,np->node_type,
        parent_index,child_index1,child_index2,np->material_head->matnum,
        np->nedges+np->ninterfaceneighbors);
      tot_coupling_count+=np->nedges+np->ninterfaceneighbors;
/*
      printf(" %u %u\n",node_index,np->nedges+np->ninterfaceneighbors);
*/
    }
  }

  /* next output child(interface) nodes */
  for (i=0; i<node_count; i++)
  {
    np=nodes[i];

    if (np->node_type==INTERFACE) /* output interface node as two child nodes */
    {
      parent_index=i+1;
      node_index=np->node_material[0].child_index+1;
      child_index1=node_index;
      child_index2=node_index;
      fprintf(outfile,"%u %.17g %.17g %.17g %.17g %u %u %u %u %d %u\n",
        node_index,np->x,np->y,np->z,np->node_material[0].volume,CHILD,
        parent_index,child_index1,child_index2,np->node_material[0].matnum,
        np->node_material[0].nedges);
      tot_coupling_count+=np->node_material[0].nedges;
/*
      printf(" %u %u\n",node_index,np->node_material[0].nedges);
*/

      node_index=np->node_material[1].child_index+1;
      child_index1=node_index;
      child_index2=node_index;
      fprintf(outfile,"%u %.17g %.17g %.17g %.17g %u %u %u %u %d %u\n",
        node_index,np->x,np->y,np->z,np->node_material[1].volume,CHILD,
        parent_index,child_index1,child_index2,np->node_material[1].matnum,
        np->node_material[1].nedges);
      tot_coupling_count+=np->node_material[1].nedges;
/*
      printf(" %u %u\n",node_index,np->node_material[1].nedges);
*/
    }
  }

/*
  printf(" part 1 total coupling count = %u\n",tot_coupling_count);
*/


  /* output part 2, node edge list */

  tot_coupling_count=0;
  /* first output parent(interface), interior, and world_boundary edge list */
  for (i=0; i<node_count; i++)
  {
    np=nodes[i];
    coupling_count=0;
    for (elp=np->edge_head; elp!=NULL; elp=elp->next)
    {
      ep=elp->edge;
      if (ep->edge_type!=VIRTUAL)
      {
        /* determine index of neighbor to node i */
        ni=ep->node_index[0];
        if (i==ni)
        {
          ni=ep->node_index[1];
        }
        fprintf(outfile,"%u %u %.17g %.17g\n",
          i+1,ni+1,ep->facet->area,ep->length);
/*
        printf(" %u %u %u\n",i+1,ni+1,ep->edge_type);
*/
        coupling_count++;
        tot_coupling_count++;
        /* also output edge connecting from interior node
           to child interface node of matching material */
        if ((np->node_type==INTERIOR) && (ep->edge_type==INTERIOR))
        {
          if (nodes[ni]->node_type==INTERFACE)
          {
            np2=nodes[ni];
            for (k=0;k<np2->material_count;k++)
            {
              if (np->material_head->matnum==np2->node_material[k].matnum)
              {
                ni=np2->node_material[k].child_index+1;
              }
            }
            fprintf(outfile,"%u %u %.17g %.17g\n",
              i+1,ni,ep->facet->area,ep->length);
/*
            printf(" %u %u %u\n",i+1,ni,ep->edge_type);
*/
            coupling_count++;
            tot_coupling_count++;
          }
        }
      }
    }
/*
    printf(" %u %u\n",i+1,coupling_count);
*/
  }

  /* finally output child(interface) edge list */
  for (i=0; i<node_count; i++)
  {
    np=nodes[i];
    if (np->node_type==INTERFACE)
    {
      for (j=0;j<np->material_count;j++)
      {
        coupling_count=0;
        node_index=np->node_material[j].child_index+1;
        for (elp=np->node_material[j].edge_head; elp!=NULL; elp=elp->next)
        {
          ep=elp->edge;

          if (ep->edge_type==INTERIOR)
          {
            ni=ep->node_index[0];
            if (i==ni)
            {
              ni=ep->node_index[1];
            }
            np2=nodes[ni];

            if (np2->node_type==INTERIOR)
            {
              fprintf(outfile,"%u %u %.17g %.17g\n",
                node_index,ni+1,ep->facet->area,ep->length);
/*
              printf(" %u %u %u\n",node_index,ni+1,ep->edge_type);
*/
            }
            else if (np2->node_type==INTERFACE)
            {
              for (k=0; k<np2->material_count; k++)
              {
                if (np->node_material[j].matnum==np2->node_material[k].matnum)
                {
                  ni=np2->node_material[k].child_index+1;
                  fprintf(outfile,"%u %u %.17g %.17g\n",
                    node_index,ni,ep->facet->area,ep->length);
/*
                  printf(" %u %u %u\n",node_index,ni,ep->edge_type);
*/
                }
              }
            }
            coupling_count++;
            tot_coupling_count++;
          }
          else if (ep->edge_type==INTERFACE)
          {
            ni=ep->node_index[0];
            if (i==ni)
            {
              ni=ep->node_index[1];
            }

            np2=nodes[ni];
            for (k=0;k<np2->material_count;k++)
            {
              if (np->node_material[j].matnum==np2->node_material[k].matnum)
              {
                ni=np2->node_material[k].child_index+1;
              }
            }

            edge_matindex=0;
            for (k=0;k<ep->material_count;k++)
            {
              if (np->node_material[j].matnum==ep->edge_material[k].matnum)
              {
                edge_matindex=k;
              }
            }
            fprintf(outfile,"%u %u %.17g %.17g\n",
              node_index,ni,ep->edge_material[edge_matindex].facet->area,ep->length);
/*
            printf(" %u %u %u\n",node_index,ni,ep->edge_type);
*/
            coupling_count++;
            tot_coupling_count++;
          }
        }
/*
        printf(" %u %u\n",node_index,coupling_count);
*/
      }
    }
  }

/*
  printf(" part 2 total coupling count = %u\n",tot_coupling_count);
*/


  fprintf(stderr,"tetgen2FV: finished output edsim data.\n");

  fclose(outfile);

  return(0);

}

