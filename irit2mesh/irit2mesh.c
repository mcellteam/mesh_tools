#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "irit2mesh.h"

extern FILE *yyin;
char *infile;
int line_num;
int total_vertex_count,total_face_count;
struct object *world_obj;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{
  struct polyhedron *php;
  struct meta_object *mop;

	if (argc<2) {
      	  printf("Usage: %s in_file_name [out_file_name]\n",argv[0]);
	  exit(1);
	}

	infile=argv[1];

	world_obj=NULL;
	if ((world_obj=(struct object *)malloc(sizeof(struct object)))==NULL) {
	  printf("irit2mesh: cannot store world object\n");
	  exit(1);
        }
        world_obj->next=NULL;
        world_obj->parent=NULL;
        world_obj->first_child=NULL;
        world_obj->last_child=NULL;
        world_obj->name="world metaobject";
        world_obj->object_type=META_OBJ;
        world_obj->contents=NULL;

	if ((yyin=fopen(infile,"r"))==NULL) {
	  printf("irit2mesh: error opening file: %s\n",infile);
	  exit(1);
	} 
	fflush(stdout);
	if (yyparse()) {
	  printf("irit2mesh: error parsing file: %s\n",infile);
	  exit(1);
	} 
	fclose(yyin);


        php=(struct polyhedron *)world_obj->first_child->contents;

        fprintf(stderr,"Sort vertex array...");
        qsort((void *)php->vertex_array,php->n_verts,sizeof(php->vertex_array),compare_verts);
        fprintf(stderr,"Done\n");


        fprintf(stderr,"Merging duplicate vertices...\n");
        merge_duplicate_verts(php,EPS_D);
        fprintf(stderr,"Done\n");
        output_polyhedron(php);

        total_vertex_count=0;
        total_face_count=0;
/*
        output_object(world_obj);
*/

	exit(0);
}



int compare_verts(const void *p1, const void *p2)
{
  struct vertex *v1,*v2;

/*
  fprintf(stderr,"p1 = %u\n",(unsigned int)p1);
  fprintf(stderr,"p2 = %u\n",(unsigned int)p2);
*/

  v1 = * (struct vertex **)p1;
  v2 = * (struct vertex **)p2;

  if (v1->projection<v2->projection)
  {
    return(-1);
  }
  if (v1->projection>v2->projection)
  {
    return(1);
  }
  return(0);

}



/**********************************************************************
distinguishable -- reports whether two doubles are measurably different

Parameters
        a -- first double
        b -- second double
        eps -- fractional difference that we think is different

Returns
        1 if the numbers are different, 0 otherwise
**********************************************************************/

int distinguishable(double a,double b,double eps)
{
  double c;

  c=a-b;

  if (c<0) c=-c;
  if (a<0) a=-a;
  if (a<1) a=1;
  if (b<0) b=-b;

  if (b<a) eps*=a;
  else eps*=b;
  return (c>eps);
}



void merge_duplicate_verts(struct polyhedron *php, double epsilon)
{
  struct vertex *vp1,*vp2;
  int running_index;
  int identical_verts,epsilon_verts;
  int i,j;

  identical_verts=0;
  epsilon_verts=0;
  running_index=0;
  for(i=0; i<php->n_verts; i++)
  {
    vp1=php->vertex_array[i];
    if (vp1->merged_index==-1)
    {
/*
      fprintf(stderr,"%d running_index = %d\n",i,running_index);
*/
      vp1->merged_index=running_index;
      php->merged_verts++;
      for (j=i+1; j<php->n_verts && !distinguishable(vp1->projection,php->vertex_array[j]->projection,epsilon); j++)
      {
        vp2=php->vertex_array[j];
        if ( !distinguishable(vp1->vertex->x,vp2->vertex->x,epsilon)
             && !distinguishable(vp1->vertex->y,vp2->vertex->y,epsilon)
             && !distinguishable(vp1->vertex->z,vp2->vertex->z,epsilon) )
        {
          vp2->merged_index=running_index;
          if (vp1->vertex->x==vp2->vertex->x
              && vp1->vertex->y==vp2->vertex->y
              && vp1->vertex->z==vp2->vertex->z)
          {
            identical_verts++;
          }
          else
          {
            epsilon_verts++;
          }
        }
      }
      running_index++;
    }
  }

  fprintf(stderr,"  identical_verts = %d\n",identical_verts);
  fprintf(stderr,"  epsilon_verts = %d\n",epsilon_verts);
  fprintf(stderr,"  unique_verts = %d\n",php->merged_verts);

  return;
}



void output_nurbs(op)
  struct object *op;
{
  struct nurbs *nrbp;
  struct double_list *dlp;
  struct ctlpt *cpp;
  struct ctlpt_list *clp;
  double val,umin,umax,vmin,vmax;
  int i,upts,vpts,uorder,vorder;

  nrbp=(struct nurbs *)op->contents;

  printf("NuPatch\n");

  upts=nrbp->upts;
  uorder=nrbp->uorder;
  printf("\t%d %d\n",upts,uorder);
  printf("\t[ ");
  i=0;
  dlp=nrbp->uknots;
  while(dlp) {
    val=dlp->val;
    if(i==uorder-1) {
      umin=val;
    }
    if(i==upts) {
      umax=val;
    }
    i++;
    printf("%.9g ",val);
    dlp=dlp->next;
  }
  printf("]\n");
  printf("\t%.9g %.9g\n",umin,umax);

  vpts=nrbp->vpts;
  vorder=nrbp->vorder;
  printf("\t%d %d\n",vpts,vorder);
  printf("\t[ ");
  i=0;
  dlp=nrbp->vknots;
  while(dlp) {
    val=dlp->val;
    if(i==vorder-1) {
      vmin=val;
    }
    if(i==vpts) {
      vmax=val;
    }
    i++;
    printf("%.9g ",val);
    dlp=dlp->next;
  }
  printf("]\n");
  printf("\t%.9g %.9g\n",vmin,vmax);

  clp=nrbp->ctlpts;
  if (clp) {
    if (clp->ctlpt->rational) {
      printf("\t\"Pw\" [\n");
    }
    else {
      printf("\t\"P\" [\n");
    }
  }
  while(clp) {
    cpp=clp->ctlpt;
    printf("\t%.9g %.9g %.9g",cpp->x,cpp->y,cpp->z);
    if(cpp->rational) {
      printf(" %.9g\n",cpp->w);
    }
    else {
      printf("\n");
    }
    clp=clp->next;
  }
  printf("\t]\n");
  
  return;
}



void output_polyhedron(struct polyhedron *php)
{
  struct vertex *vp;
  struct vertex_list *vlp;
  struct vector3 *vecp;
  struct polygon_list *plp;
  struct polygon *pop;
  int running_index;
  int total_face_count;
  int i;

  running_index=0;
  for (i=0; i<php->n_verts; i++)
  {
    vp=php->vertex_array[i];
    if (vp->merged_index>=running_index)
    {
      vecp=vp->vertex;
      printf("Vertex %d %.15g %.15g %.15g\n",running_index+1,vecp->x,vecp->y,vecp->z);
      running_index++;
    }
  }

  total_face_count=0;
  for (plp=php->polygon_list; plp!=NULL; plp=plp->next)
  {
    pop=plp->polygon;

    total_face_count++;
    printf("Face %d",total_face_count);
    for (vlp=pop->vertex_list; vlp!=NULL; vlp=vlp->next)
    {
      printf(" %d",vlp->vertex->merged_index+1);
    }
    printf("\n");
  }

  return;
}

