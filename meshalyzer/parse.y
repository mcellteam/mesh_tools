%{
#include <stdlib.h> 
#include <stdio.h> 
#include <string.h> 
#include <float.h>
#include <math.h>
#include "vector.h"
#include "meshalyzer.h"

#include "lex.c"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *yyin;
extern char *infile;
extern struct object *world_obj;

int ival;
double rval;
char *cval;
char *a_str,*rem_str;
char c,fmt_str[128],time_str[128];
char err_msg[128];
struct object *op;
struct surface_mesh *smp;
struct polygon_list *plp,*polygon_head,*polygon_tail;
struct polygon *pop;
struct nurbs *nrbp;
struct double_list *dlp_head,*dlp;
struct ctlpt_list *clp_head,*clp;
struct ctlpt *cpp;
struct vertex *vp;
struct vertex_list *sm_vlp,*sm_vertex_head,*sm_vertex_tail;
struct vertex_list *vlp,*vertex_head,*vertex_tail;
struct vector3 *vecp;
struct vector3 rand_vec;
double vol_infinity;
int vertex_index;

char *my_strdup(s)
  char *s;
{
  char *temp;

  if ((temp=(char *)malloc(strlen(s)+1))!=NULL) {
    strcpy(temp,s);
  }
  return(temp);
}


%}



%union {
int tok;
char *str;
double dbl;
struct vector3 *vec;
struct object *obj;
} 

%token <tok> REAL INTEGER VERTEX FACE
%type <dbl> int_arg real_arg num_arg

%right '='
%left '+' '-'
%left '*' '/'
%left UNARYMINUS



%%



mesh_format:
{
    if ((op=(struct object *)malloc(sizeof(struct object)))==NULL) {
      yyerror("Cannot store object");
      return(1);
    }

    if (world_obj->first_child==NULL) {
      world_obj->first_child=op;
    }
    if (world_obj->last_child!=NULL) {
      world_obj->last_child->next=op;
    }
    world_obj->last_child=op;

    op->next=NULL;
    op->parent=world_obj;
    op->first_child=NULL;
    op->last_child=NULL;
    op->name="mesh";
    op->object_type=POLY;
    op->contents=NULL;

    if ((smp=(struct surface_mesh *)malloc(sizeof(struct surface_mesh)))==NULL) {
      yyerror("Cannot store polygon object");
      return(1);
    }
    smp->n_polys=0;
    smp->n_verts=0;
    smp->max_verts=0;
    smp->merged_verts=0;
    smp->polygon_list=NULL;
    smp->vertex_list=NULL;
    smp->vertex_array=NULL;
    vol_infinity=sqrt(DBL_MAX)/4;
    smp->llf.x=vol_infinity;
    smp->llf.y=vol_infinity;
    smp->llf.z=vol_infinity;
    smp->urb.x=-vol_infinity;
    smp->urb.y=-vol_infinity;
    smp->urb.z=-vol_infinity;
    plp=NULL;
    polygon_head=NULL;
    polygon_tail=NULL;
    sm_vlp=NULL;
    sm_vertex_head=NULL;
    sm_vertex_tail=NULL;
    rand_vec.x=0.236416584579274058342;
    rand_vec.y=0.927225593011826276779;
    rand_vec.z=0.389099507126957178116;
    normalize(&rand_vec);

    op->contents=(void *)smp;
}
	vertex_list
{

  smp->vertex_list=sm_vertex_head;
  if ((smp->vertex_array=(struct vertex **)malloc(smp->max_verts*sizeof(struct vertex *)))==NULL) {
    yyerror("Cannot store vertex array");
    return(1);
  }

  for (sm_vlp=smp->vertex_list; sm_vlp!=NULL; sm_vlp=sm_vlp->next)
  {
    smp->vertex_array[sm_vlp->vertex->vertex_index-1]=sm_vlp->vertex;
  }
}
	face_list
{
    smp->polygon_list=polygon_head;

    fprintf(stderr,"\nNumber of vertices: %d  Number of polygons: %d\n",smp->n_verts,smp->n_polys);

};


int_arg: INTEGER {$$=(double)ival;}
;


real_arg: REAL {$$=rval;}
;


num_arg: INTEGER {$$=(double)ival;}
        | REAL {$$=rval;}
;


vertex_list: vertex
	| vertex_list vertex
;


vertex:  VERTEX int_arg num_arg num_arg num_arg
{
  if ((vecp=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
    yyerror("Cannot store vertex");
    return(1);
  }
  vecp->x=$<dbl>3;
  vecp->y=$<dbl>4;
  vecp->z=$<dbl>5;
  if (vecp->x < smp->llf.x) {
    smp->llf.x=vecp->x;
  }
  if (vecp->y < smp->llf.y) {
    smp->llf.y=vecp->y;
  }
  if (vecp->z < smp->llf.z) {
    smp->llf.z=vecp->z;
  }
  if (vecp->x > smp->urb.x) {
    smp->urb.x=vecp->x;
  }
  if (vecp->y > smp->urb.y) {
    smp->urb.y=vecp->y;
  }
  if (vecp->z > smp->urb.z) {
    smp->urb.z=vecp->z;
  }

  if ((vp=(struct vertex *)malloc(sizeof(struct vertex)))==NULL) {
    yyerror("Cannot store vertex list");
    return(1);
  }
  vertex_index=$<dbl>2;
  if (vertex_index>smp->max_verts) {
    smp->max_verts=vertex_index;
  }

  vp->vertex_index=vertex_index;
  vp->vertex_count=smp->n_verts++;
  vp->merged_index=-1;
  vp->vertex=vecp;
  vp->normal=NULL;
  vp->projection=dot_prod(vp->vertex,&rand_vec);

  if ((sm_vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    yyerror("Cannot store vertex list");
    return(1);
  }
  sm_vlp->vertex=vp;
  if (sm_vertex_tail==NULL) {
    sm_vertex_tail=sm_vlp;
  }
  sm_vertex_tail->next=sm_vlp;
  sm_vlp->next=NULL;
  sm_vertex_tail=sm_vlp;
  if (sm_vertex_head==NULL) {
    sm_vertex_head=sm_vlp;
  }
};


face_list: face
	| face_list face
;


face: FACE int_arg int_arg int_arg int_arg
{
  if ((plp=(struct polygon_list *)malloc
             (sizeof(struct polygon_list)))==NULL) {
    yyerror("Cannot store polygon object");
    return(1);
  }
  if (polygon_tail==NULL) {
    polygon_tail=plp;
  }
  polygon_tail->next=plp;
  plp->next=NULL;
  polygon_tail=plp;
  if (polygon_head==NULL) {
    polygon_head=plp;
  }

  if ((pop=(struct polygon *)malloc(sizeof(struct polygon)))==NULL) {
    yyerror("Cannot store polygon object");
    return(1);
  }
  pop->polygon_index=$<dbl>2;
  pop->area=0.0;
  pop->n_verts=3;

  if ((pop->vertex_array=(struct vertex **)malloc(pop->n_verts*sizeof(struct vertex *)))==NULL) {
    yyerror("Cannot store polygon vertices");
    return(1);
  }

  pop->vertex_array[0]=smp->vertex_array[(int)$<dbl>3-1];
  pop->vertex_array[1]=smp->vertex_array[(int)$<dbl>4-1];
  pop->vertex_array[2]=smp->vertex_array[(int)$<dbl>5-1];

  plp->polygon=pop;
  smp->n_polys++;
};



%%



yyerror(s)
char *s;
{
	fprintf(stderr,"Irit2mesh: error on line: %d of file: %s  %s\n",
	        line_num,infile,s);
	fflush(stderr);
	return(1);
}

