#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "stl2mesh.h"

#define NUM_DICES 20
#define EPSILON_1 1.0e-12
#define EPSILON_2 1.0e-9

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
  int i,j,k,l;
  char *argstr;

        if (argc==1)
        {
          infile="stdin";
          yyin=stdin;
        }
        else
        {
          argstr=argv[1];
          if (strcmp(argstr,"-h")==0)
          {
      	    fprintf(stderr,"\nRead ascii stl file and convert to mesh format\n");
      	    fprintf(stderr,"  Read from stdin if stl_file_name is absent\n");
      	    fprintf(stderr,"  Output is written to stdout\n\n");
      	    fprintf(stderr,"  Usage: %s [-h] [stl_file_name]\n\n",argv[0]);
	    exit(1);
          }
	  infile=argv[1];
	  if ((yyin=fopen(infile,"r"))==NULL) {
	    fprintf(stderr,"stl2mesh: error opening file: %s\n",infile);
	    exit(1);
	  } 
        }

	world_obj=NULL;
	if ((world_obj=(struct object *)malloc(sizeof(struct object)))==NULL) {
	  fprintf(stderr,"stl2mesh: cannot store world object\n");
	  exit(1);
        }
        world_obj->next=NULL;
        world_obj->parent=NULL;
        world_obj->first_child=NULL;
        world_obj->last_child=NULL;
        world_obj->name="world metaobject";
        world_obj->object_type=META_OBJ;
        world_obj->contents=NULL;

	fflush(stdout);
	if (yyparse()) {
	  fprintf(stderr,"stl2mesh: error parsing file: %s\n",infile);
	  exit(1);
	} 
	fclose(yyin);

        php=(struct polyhedron *)world_obj->first_child->contents;

        fprintf(stderr,"Sort vertex array...");
        qsort((void *)php->vertex_array,php->n_verts,sizeof(php->vertex_array),compare_verts);
        fprintf(stderr,"Done\n");

/*
        for (i=0; i<php->n_verts; i++)
        {
          fprintf(stderr,"%d  proj = %.9g\n",i,php->vertex_array[i]->projection);
        }
*/

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



