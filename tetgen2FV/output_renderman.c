#include <stdlib.h> 
#include <stdio.h> 
#include <string.h> 
#include <sys/types.h>
#include <math.h>
#include "strfunc.h"
#include "tetgen2FV.h"
#include "voronoi.h"
#include "output_renderman.h"

#ifdef DEBUG
#define no_printf printf
#endif

#define MY_PI 3.14159265358979323846


int output_renderman_full(struct tet_mesh *tet_mesh, char *outfile_name)
{

  FILE *outfile;
  struct hash_table *hp;
  struct voronoi_facet *vfp;
  struct node **node_array;
  struct vector3 **voronoi_node_array;
  struct node *np,*np1,*np2,*np3;
  struct node_list *node_head;
  struct edge_list *elp;
  struct edge *ep;
  struct edge_material_info *emip;
  struct node_material_info *nmip;
  struct tet_list *tlp,*tet_head;
  struct polygon_list *plp;
  struct polygon *pop;
  struct vector3 *vp1,*vp2;
  struct vector3 glyph_axis,side_axis,rotation_axis;
  double radius,rad,dot,rotation_angle,axis_length,side_length;
  int output_mat;
  u_int node_count;
  u_int tet_count;
  u_int edge_count;
  u_int boundaryface_count;
  u_int voronoi_node_count;
  u_int voronoi_edge_count;
  u_int voronoi_boundary_edge_count;
  u_int voronoi_triangle_count;
  u_int i,j,k;

  fprintf(stderr,"tetgen2FV: output RenderMan mesh:\n");

  node_count=tet_mesh->nnodes;
  node_head=tet_mesh->node_head;
  node_array=tet_mesh->nodes;
  voronoi_node_array=tet_mesh->voronoi_nodes;
  tet_count=tet_mesh->ntets;
  voronoi_node_count=tet_mesh->ntets+tet_mesh->nboundaryfaces+tet_mesh->nboundaryedges+tet_mesh->nboundarynodes;
  boundaryface_count=tet_mesh->nboundaryfaces;
  voronoi_boundary_edge_count=tet_mesh->voronoi_boundary_nedges_tot;
  tet_head=tet_mesh->tet_head;
 
  edge_count=tet_mesh->nedges;
  /* multiply the edge and triangle counts by 2 because 
     we'll duplicate the shared voronoi facets as we output them */
  voronoi_edge_count=(2*tet_mesh->voronoi_facet_nedges_tot)+tet_mesh->voronoi_cap_nedges_tot;
  voronoi_triangle_count=(2*tet_mesh->voronoi_facet_ntriangles_tot)+tet_mesh->voronoi_cap_ntriangles_tot;


  if ((outfile=fopen(outfile_name,"w"))==NULL)
  {
    fprintf(stderr,"Cannot open RenderMan outfile %s",outfile_name);
    return(1);
  }

  radius = 0.0005;
  rad=360.0/(2.0*MY_PI);

#if 0
  /* output tet mesh nodes */
  fprintf(outfile,"# tet mesh nodes\n");
  for (i=0; i<node_count; i++)
  {
    np=node_array[i];
    fprintf(outfile,"AttributeBegin\n");
    fprintf(outfile,"  Translate %.9g %.9g %.9g\n",np->x,np->y,np->z);
    fprintf(outfile,"  Sphere %.9g %.9g %.9g 360\n",radius,-radius,radius);
    fprintf(outfile,"AttributeEnd\n");
  }


  /* output tet mesh tet edges */
  fprintf(outfile,"# tet mesh edges\n");
  for (i=0; i<tet_mesh->hashsize; i++)
  {
    for (hp=tet_mesh->edge_hashtab[i]; hp!=NULL; hp=hp->next)
    {
      ep=(struct edge *)hp->contents;

      if (ep->edge_type!=VIRTUAL)
      {
        /* output all tets */
        vp1=(struct vector3 *)node_array[ep->node_index[0]];
        vp2=(struct vector3 *)node_array[ep->node_index[1]];
        
        glyph_axis.x = 0;
        glyph_axis.y = 0;
        glyph_axis.z = 1;
        vectorize(vp1,vp2,&side_axis);
        side_length = vect_length(&side_axis);
        normalize(&side_axis);
        cross_prod(&glyph_axis,&side_axis,&rotation_axis);
        dot=dot_prod(&glyph_axis,&side_axis)*0.999999999999;
        rotation_angle = acos(dot)*rad;
        axis_length = vect_length(&rotation_axis);
        fprintf(outfile,"AttributeBegin\n");
        fprintf(outfile,"  Translate %.9g %.9g %.9g\n",vp1->x,vp1->y,vp1->z);
        if (rotation_angle!=0 && axis_length>0) {
          fprintf(outfile,"  Rotate %.9g %.9g %.9g %.9g\n",rotation_angle,
            rotation_axis.x,rotation_axis.y,rotation_axis.z);
        }
        fprintf(outfile,"  Cylinder %.9g 0 %.9g 360\n",radius,side_length);
        fprintf(outfile,"AttributeEnd\n");
      }
      

      /* output all tets */
/*
      if (ep->edge_type==BOUNDARY || ep->edge_type==INTERFACE)
      {
      }
*/

    }
  }


  /* output tet mesh boundary faces */
  fprintf(outfile,"# tet mesh boundary faces\n");
  for (plp=tet_mesh->boundary_head; plp!=NULL; plp=plp->next)
  {
    pop=plp->polygon;

    np1=node_array[pop->node_index[0]];
    np2=node_array[pop->node_index[1]];
    np3=node_array[pop->node_index[2]];

    fprintf(outfile,"AttributeBegin\n");
    fprintf(outfile,"  Polygon\n");
    fprintf(outfile,"    \"P\" [ %.9g %.9g %.9g  %.9g %.9g %.9g  %.9g %.9g %.9g ]\n",np1->x,np1->y,np1->z,np2->x,np2->y,np2->z,np3->x,np3->y,np3->z);
    fprintf(outfile,"AttributeEnd\n");
  }
#endif


#if 0
  /* output voronoi mesh nodes */
  fprintf(outfile,"# voronoi mesh nodes\n");
  for (i=0; i<voronoi_node_count; i++)
  {
    vp1=voronoi_node_array[i];
    if(vp1->y>0)
    {
      fprintf(outfile,"AttributeBegin\n");
      fprintf(outfile,"  Translate %.9g %.9g %.9g\n",vp1->x,vp1->y,vp1->z);
      fprintf(outfile,"  Sphere %.9g %.9g %.9g 360\n",radius,-radius,radius);
      fprintf(outfile,"AttributeEnd\n");
    }
  }


  /* output voronoi mesh edges */
  fprintf(outfile,"# voronoi mesh edges\n");
  for (i=0; i<tet_mesh->hashsize; i++)
  {
    for (hp=tet_mesh->edge_hashtab[i]; hp!=NULL; hp=hp->next)
    {
      ep=(struct edge *)hp->contents;
      if(node_array[ep->node_index[0]]->y>0)
      {
        for (vfp=ep->facet; vfp!=NULL; vfp=vfp->next)
        {
          for (j=0; j<vfp->nnodes; j++)
          {
            vp1=voronoi_node_array[vfp->node_index[j]];
            if (j<vfp->nnodes-1)
            {
              vp2=voronoi_node_array[vfp->node_index[j+1]];
            }
            else
            {
              vp2=voronoi_node_array[vfp->node_index[0]];
            }
            glyph_axis.x = 0;
            glyph_axis.y = 0;
            glyph_axis.z = 1;
            vectorize(vp1,vp2,&side_axis);
            side_length = vect_length(&side_axis);
            normalize(&side_axis);
            cross_prod(&glyph_axis,&side_axis,&rotation_axis);
            dot=dot_prod(&glyph_axis,&side_axis)*0.999999999999;
            rotation_angle = acos(dot)*rad;
            axis_length = vect_length(&rotation_axis);
            fprintf(outfile,"AttributeBegin\n");
            fprintf(outfile,"  Translate %.9g %.9g %.9g\n",vp1->x,vp1->y,vp1->z);
            if (rotation_angle!=0 && axis_length>0) {
              fprintf(outfile,"  Rotate %.9g %.9g %.9g %.9g\n",rotation_angle,
                rotation_axis.x,rotation_axis.y,rotation_axis.z);
            }
            fprintf(outfile,"  Cylinder %.9g 0 %.9g 360\n",radius,side_length);
            fprintf(outfile,"AttributeEnd\n");
          }
        }
      }
    }
  }
#endif


#if 0
#endif
  /* output voronoi mesh facets as polygons */
  fprintf(outfile,"# voronoi mesh facets\n");
  output_mat=-10;
  for (i=0; i<tet_mesh->hashsize; i++)
  {
    for (hp=tet_mesh->edge_hashtab[i]; hp!=NULL; hp=hp->next)
    {
      ep=(struct edge *)hp->contents;
      if(node_array[ep->node_index[0]]->y>0)
      {
        if (ep->edge_type==INTERFACE)
        {
          emip=ep->edge_material;
          for (j=0; j<ep->material_count; j++)
          {
            if (emip[j].matnum==output_mat)
            {
              for (vfp=emip[j].facet; vfp!=NULL; vfp=vfp->next)
              {
        
                fprintf(outfile,"AttributeBegin\n");
                fprintf(outfile,"  Polygon\n");
                fprintf(outfile,"    \"P\" [");
                for (k=0; k<vfp->nnodes; k++)
                {
                  vp1=voronoi_node_array[vfp->node_index[k]];
                  fprintf(outfile," %.9g %.9g %.9g ",vp1->x,vp1->y,vp1->z);
                }
                fprintf(outfile,"]\n");
                fprintf(outfile,"AttributeEnd\n");
              }
            }
          }
        }
        else if (ep->edge_type!=VIRTUAL) /* INTERIOR, BOUNDARY edges */
        {
          if (ep->matnum==output_mat)
          {
            for (vfp=ep->facet; vfp!=NULL; vfp=vfp->next)
            {
  
              fprintf(outfile,"AttributeBegin\n");
              fprintf(outfile,"  Polygon\n");
              fprintf(outfile,"    \"P\" [");
              for (j=0; j<vfp->nnodes; j++)
              {
                vp1=voronoi_node_array[vfp->node_index[j]];
                fprintf(outfile," %.9g %.9g %.9g ",vp1->x,vp1->y,vp1->z);
              }
              fprintf(outfile,"]\n");
              fprintf(outfile,"AttributeEnd\n");
            }
          }
        }
        else /* VIRTUAL edges */
        {
          np=node_array[ep->node_index[0]];
          if (np->node_type==INTERFACE)
          {
            nmip=np->node_material;
            for (j=0; j<np->material_count; j++)
            {
              if (nmip[j].matnum==output_mat)
              {
                for (vfp=ep->facet; vfp!=NULL; vfp=vfp->next)
                {
          
                  fprintf(outfile,"AttributeBegin\n");
                  fprintf(outfile,"  Polygon\n");
                  fprintf(outfile,"    \"P\" [");
                  for (k=0; k<vfp->nnodes; k++)
                  {
                    vp1=voronoi_node_array[vfp->node_index[k]];
                    fprintf(outfile," %.9g %.9g %.9g ",vp1->x,vp1->y,vp1->z);
                  }
                  fprintf(outfile,"]\n");
                  fprintf(outfile,"AttributeEnd\n");
                }
              }
            }
          }
          else /* BOUNDARY node */
          {
            if (-30==output_mat)
            {
              for (vfp=ep->facet; vfp!=NULL; vfp=vfp->next)
              {
          
                fprintf(outfile,"AttributeBegin\n");
                fprintf(outfile,"  Polygon\n");
                fprintf(outfile,"    \"P\" [");
                for (k=0; k<vfp->nnodes; k++)
                {
                  vp1=voronoi_node_array[vfp->node_index[k]];
                  fprintf(outfile," %.9g %.9g %.9g ",vp1->x,vp1->y,vp1->z);
                }
                fprintf(outfile,"]\n");
                fprintf(outfile,"AttributeEnd\n");
              }
            }
          }
        }
      }
    }
  }



#if 0

  /* output tet indices */
  fi=0;
  tlp=tet_head;
  while (tlp!=NULL)
  {
    fwrite(&fi,sizeof(float),1,outfile);
    fi++;
    tlp=tlp->next;
  }


  /* output tet material */
  tlp=tet_head;
  while (tlp!=NULL)
  {
    fi=tlp->tet->matnum;
    fwrite(&fi,sizeof(float),1,outfile);
    tlp=tlp->next;
  }


  /* output multi-material node flags */
  for (i=0; i<node_count; i++)
  {
    fi=node_array[i]->material_count;
    fwrite(&fi,sizeof(float),1,outfile);
  }


  /* output tet mesh tet circumcenters, voronoi mesh nodes */
  tlp=tet_head;
  while (tlp!=NULL)
  {
    nx=tlp->tet->cent.x;
    ny=tlp->tet->cent.y;
    nz=tlp->tet->cent.z;
    fwrite(&nx,sizeof(float),1,outfile);
    fwrite(&ny,sizeof(float),1,outfile);
    fwrite(&nz,sizeof(float),1,outfile);
    tlp=tlp->next;
  }


  /* output tet circumcenter indices */
  fi=0;
  tlp=tet_head;
  while (tlp!=NULL)
  {
    fwrite(&fi,sizeof(float),1,outfile);
    fi++;
    tlp=tlp->next;
  }


  /* output tet/circumcenter delaunay status */
  tlp=tet_head;
  while (tlp!=NULL)
  {
    fi=tlp->tet->delaunay;
    fwrite(&fi,sizeof(float),1,outfile);
    tlp=tlp->next;
  }


  /* output tet mesh connections by tet mesh node */
  tlp=tet_head;
  while (tlp!=NULL)
  {
    for (i=0; i<4; i++)
    {
      fwrite(&tlp->tet->node_index[0],sizeof(u_int),1,outfile);
      fwrite(&tlp->tet->node_index[1],sizeof(u_int),1,outfile);
      fwrite(&tlp->tet->node_index[2],sizeof(u_int),1,outfile);
      fwrite(&tlp->tet->node_index[3],sizeof(u_int),1,outfile);
    }
    tlp=tlp->next;
  }


  /* output tet mesh tet to tet mesh node mapping */
  tlp=tet_head;
  while (tlp!=NULL)
  {
    for (i=0; i<4; i++)
    {
      fi=tlp->tet->node_index[i];
      fwrite(&fi,sizeof(float),1,outfile);
    }
    tlp=tlp->next;
  }


  /* output boundary face circumcenters */
  for (plp=tet_mesh->boundary_head; plp!=NULL; plp=plp->next)
  {
    pop=plp->polygon;
    nx=pop->cent.x;
    ny=pop->cent.y;
    nz=pop->cent.z;
    fwrite(&nx,sizeof(float),1,outfile);
    fwrite(&ny,sizeof(float),1,outfile);
    fwrite(&nz,sizeof(float),1,outfile);
  }
  

  /* output voronoi mesh nodes */
  for (i=0; i<voronoi_node_count; i++)
  {
    vp1=tet_mesh->voronoi_nodes[i];
    nx=vp1->x;
    ny=vp1->y;
    nz=vp1->z;
    fwrite(&nx,sizeof(float),1,outfile);
    fwrite(&ny,sizeof(float),1,outfile);
    fwrite(&nz,sizeof(float),1,outfile);
  }


  /* output voronoi mesh facets as triangles */
  for (i=0;i<node_count;i++)
  {
    for (elp=node_array[i]->edge_head; elp!=NULL; elp=elp->next)
    {
      for (vfp=elp->edge->facet; vfp!=NULL; vfp=vfp->next)
      {
        for (j=0; j<vfp->nnodes-2; j++)
        {
          n1=vfp->node_index[0];
          n2=vfp->node_index[j+1];
          n3=vfp->node_index[j+2];
          fwrite(&n1,sizeof(u_int),1,outfile);
          fwrite(&n2,sizeof(u_int),1,outfile);
          fwrite(&n3,sizeof(u_int),1,outfile);
        }
      }
    }
  }


  /* output voronoi mesh triangle to tet mesh node mapping */
  for (i=0;i<node_count;i++)
  {
    for (elp=node_array[i]->edge_head; elp!=NULL; elp=elp->next)
    {
      for (vfp=elp->edge->facet; vfp!=NULL; vfp=vfp->next)
      {
        for (j=0; j<vfp->nnodes-2; j++)
        {
          fi=i;
          fwrite(&fi,sizeof(float),1,outfile);
        }
      }
    }
  }


  /* output voronoi mesh triangle facet indices */
  for (i=0;i<node_count;i++)
  {
    k=0;
    for (elp=node_array[i]->edge_head; elp!=NULL; elp=elp->next)
    {
      for (vfp=elp->edge->facet; vfp!=NULL; vfp=vfp->next)
      {
        for (j=0; j<vfp->nnodes-2; j++)
        {
          fi=k;
          fwrite(&fi,sizeof(float),1,outfile);
        }
      }
      k++;
    }
  }


  /* output voronoi mesh edges as line segments */
  for (i=0;i<node_count;i++)
  {
    for (elp=node_array[i]->edge_head; elp!=NULL; elp=elp->next)
    {
      for (vfp=elp->edge->facet; vfp!=NULL; vfp=vfp->next)
      {
        for (j=0; j<vfp->nnodes; j++)
        {
          if (j<vfp->nnodes-1)
          {
            n1=vfp->node_index[j];
            n2=vfp->node_index[j+1];
          }
          else
          {
            n1=vfp->node_index[j];
            n2=vfp->node_index[0];
          }
          fwrite(&n1,sizeof(u_int),1,outfile);
          fwrite(&n2,sizeof(u_int),1,outfile);
        }
      }
    }
  }


  /* output voronoi mesh edge to tet mesh node mapping */
  for (i=0;i<node_count;i++)
  {
    for (elp=node_array[i]->edge_head; elp!=NULL; elp=elp->next)
    {
      for (vfp=elp->edge->facet; vfp!=NULL; vfp=vfp->next)
      {
        for (j=0; j<vfp->nnodes; j++)
        {
          fi=i;
          fwrite(&fi,sizeof(float),1,outfile);
        }
      }
    }
  }


  /* output voronoi mesh edge facet indices */
  for (i=0;i<node_count;i++)
  {
    k=0;
    for (elp=node_array[i]->edge_head; elp!=NULL; elp=elp->next)
    {
      for (vfp=elp->edge->facet; vfp!=NULL; vfp=vfp->next)
      {
        for (j=0; j<vfp->nnodes; j++)
        {
          fi=k;
          fwrite(&fi,sizeof(float),1,outfile);
        }
      }
      k++;
    }
  }


  /* output voronoi mesh boundary edges */
  for (i=0; i<node_count; i++)
  {
    np=node_array[i];
    if (np->node_type==BOUNDARY || np->node_type==INTERFACE)
    {
      for (plp=np->boundary_head; plp!=NULL; plp=plp->next)
      {
        if (plp->next!=NULL)
        {
          n1=plp->polygon->polygon_index;
          n2=plp->next->polygon->polygon_index;
        }
        else {
          n1=plp->polygon->polygon_index;
          n2=np->boundary_head->polygon->polygon_index;
        }
        fwrite(&n1,sizeof(u_int),1,outfile);
        fwrite(&n2,sizeof(u_int),1,outfile);
      }
    } 
  }


  /* output voronoi mesh boundary edge to tet mesh node mapping */
  for (i=0; i<node_count; i++)
  {
    np=node_array[i];
    if (np->node_type==BOUNDARY || np->node_type==INTERFACE)
    {
      for (j=0; j<np->nboundaryfaces; j++)
      {
        fi=i;
        fwrite(&fi,sizeof(float),1,outfile);
      }
    } 
  }

#endif


  fclose(outfile);

  fprintf(stderr,"tetgen2FV: finished output RenderMan mesh.\n");

  return(0);

}



int output_renderman_tet_mesh(struct tet_mesh *tet_mesh, char *outfile_name)
{

  FILE *outfile;
  struct voronoi_facet *vfp;
  struct node **node_array;
  struct node *np;
  struct node_list *node_head;
  struct edge_list *elp;
  struct tet_list *tlp,*tet_head;
  struct polygon_list *plp;
  struct polygon *pop;
  struct vector3 *vp;
  u_int node_count;
  u_int tet_count;
  u_int edge_count;
  u_int boundaryface_count;
  u_int voronoi_node_count;
  u_int voronoi_edge_count;
  u_int voronoi_boundary_edge_count;
  u_int voronoi_triangle_count;
  int data_offset;
  float nx,ny,nz;
  float fi;
  u_int n1,n2,n3;
  u_int i,j,k;
  char obj_name[64];

  fprintf(stderr,"tetgen2FV: output DX mesh:\n");

  node_count=tet_mesh->nnodes;
  node_head=tet_mesh->node_head;
  node_array=tet_mesh->nodes;
  tet_count=tet_mesh->ntets;
  voronoi_node_count=tet_mesh->ntets+tet_mesh->nboundaryfaces+tet_mesh->nboundaryedges+tet_mesh->nboundarynodes;
  boundaryface_count=tet_mesh->nboundaryfaces;
  voronoi_boundary_edge_count=tet_mesh->voronoi_boundary_nedges_tot;
  tet_head=tet_mesh->tet_head;
 
  edge_count=tet_mesh->nedges;
  /* multiply the edge and triangle counts by 2 because 
     we'll duplicate the shared voronoi facets as we output them */
  voronoi_edge_count=(2*tet_mesh->voronoi_facet_nedges_tot)+tet_mesh->voronoi_cap_nedges_tot;
  voronoi_triangle_count=(2*tet_mesh->voronoi_facet_ntriangles_tot)+tet_mesh->voronoi_cap_ntriangles_tot;

  if ((outfile=fopen(outfile_name,"w"))==NULL)
  {
    fprintf(stderr,"Cannot open DX outfile %s",outfile_name);
    return(1);
  }

  data_offset=0;
  
  sprintf(obj_name,"tet mesh nodes");
  fprintf(outfile,"# %s\n",obj_name);
  fprintf(outfile,"object \"%s\" class array type float rank 1 shape 3 items %d lsb ieee data %d\n",obj_name,node_count,data_offset);
  fprintf(outfile,"attribute \"dep\" string \"positions\"\n");
  fprintf(outfile,"#\n");
  data_offset+=node_count*3*4;

  sprintf(obj_name,"tet mesh tets");
  fprintf(outfile,"# %s\n",obj_name);
  fprintf(outfile,"object \"%s\" class array type int rank 1 shape 4 items %d lsb ieee data %d\n",obj_name,tet_count,data_offset);
  fprintf(outfile,"attribute \"dep\" string \"connections\"\n");
  fprintf(outfile,"attribute \"ref\" string \"positions\"\n");
  fprintf(outfile,"attribute \"element type\" string \"tetrahedra\"\n");
  fprintf(outfile,"#\n");
  data_offset+=tet_count*4*4;

  sprintf(obj_name,"tet indices");
  fprintf(outfile,"# %s\n",obj_name);
  fprintf(outfile,"object \"%s\" class array type float rank 0 items %d lsb ieee data %d\n",obj_name,tet_count,data_offset);
  fprintf(outfile,"attribute \"dep\" string \"connections\"\n");
  fprintf(outfile,"#\n");
  data_offset+=tet_count*1*4;

  sprintf(obj_name,"tet material");
  fprintf(outfile,"# %s\n",obj_name);
  fprintf(outfile,"object \"%s\" class array type float rank 0 items %d lsb ieee data %d\n",obj_name,tet_count,data_offset);
  fprintf(outfile,"attribute \"dep\" string \"connections\"\n");
  fprintf(outfile,"#\n");
  data_offset+=tet_count*1*4;

  sprintf(obj_name,"tet mesh tets field");
  fprintf(outfile,"# %s\n",obj_name);
  fprintf(outfile,"object \"%s\" class field\n",obj_name);
  fprintf(outfile,"component \"positions\" value \"%s\"\n",
    "tet mesh nodes");
  fprintf(outfile,"component \"connections\" value \"%s\"\n",
    "tet mesh tets");
  fprintf(outfile,"component \"tet index\" value \"%s\"\n",
    "tet indices");
  fprintf(outfile,"component \"tet material\" value \"%s\"\n",
    "tet material");
  fprintf(outfile,"#\n");

  fprintf(outfile,"object \"default\" class group\n");
  fprintf(outfile,"member \"delaunay_tets\" value \"%s\"\n",
    "tet mesh tets field");
  fprintf(outfile,"#\n");
  fprintf(outfile,"end\n");


  /* output tet mesh nodes */
  for (i=0; i<node_count; i++)
  {
    np=node_array[i];
    nx=np->x;
    ny=np->y;
    nz=np->z;
    fwrite(&nx,sizeof(float),1,outfile);
    fwrite(&ny,sizeof(float),1,outfile);
    fwrite(&nz,sizeof(float),1,outfile);
  }


  /* output tet mesh tet connections */
  tlp=tet_head;
  while (tlp!=NULL)
  {
    fwrite(&tlp->tet->node_index[0],sizeof(u_int),1,outfile);
    fwrite(&tlp->tet->node_index[1],sizeof(u_int),1,outfile);
    fwrite(&tlp->tet->node_index[2],sizeof(u_int),1,outfile);
    fwrite(&tlp->tet->node_index[3],sizeof(u_int),1,outfile);
    tlp=tlp->next;
  }


  /* output tet indices */
  fi=0;
  tlp=tet_head;
  while (tlp!=NULL)
  {
    fwrite(&fi,sizeof(float),1,outfile);
    fi++;
    tlp=tlp->next;
  }


  /* output tet material */
  tlp=tet_head;
  while (tlp!=NULL)
  {
    fi=tlp->tet->matnum;
    fwrite(&fi,sizeof(float),1,outfile);
    tlp=tlp->next;
  }



  fclose(outfile);

  fprintf(stderr,"tetgen2FV: finished output DX mesh.\n");

  return(0);

}

