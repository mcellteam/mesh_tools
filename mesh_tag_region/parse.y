%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "mesh_tag_region.h"
#include "vector.h"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *yyin;
extern char *infile;
extern char *outfile;
extern int skip_freq;
extern struct object *world_obj;
extern struct object *op;
extern double cos_threshold;
extern double vol_infinity;

int ival;
double rval;
char *cval;
char *var;
char *a_str,*rem_str;
char c,fmt_str[128],time_str[128];
char err_msg[128];
struct meta_object *mop,*curr_mop,*wmop;
struct object_list *olp;
struct polyhedron *php;
struct polygon_list *plp,*polygon_head,*polygon_tail;
struct polygon *pop;
struct nurbs *nrbp;
struct double_list *dlp_head,*dlp;
struct ctlpt_list *clp_head,*clp;
struct ctlpt *cpp;
struct vertex_list *vlp,*vertex_head,*vertex_tail,**vertex_list_array;
struct vector3 *vecp,v1,v2,normal,y_axis;
int vertex_index;
int vertex_count;
int max_vertex;
int skip_count;
int polygon_count,polygon_index;
double x,y,z,dot;
int vert_0,vert_1,vert_2;
int i;


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

%{
  #include "lex.flex.c"
%}

%output="parse.bison.c"

%token <tok> REAL INTEGER VERTEX FACE COMMENT
%type <dbl> int_arg real_arg num_arg 

%right '='
%left '+' '-'
%left '*' '/'
%left UNARYMINUS

%%

mesh_format: 
{
  skip_count=0;
  vertex_count=0;
  max_vertex=0;
  polygon_count=0;
  vlp=NULL;
  vertex_head=NULL;
  vertex_tail=NULL;
  plp=NULL;
  polygon_head=NULL;
  polygon_tail=NULL;
  if ((op=(struct object *)malloc
       (sizeof(struct object)))==NULL) {
    yyerror("Cannot store object");
    return(1);
  }
  if ((php=(struct polyhedron *)malloc
       (sizeof(struct polyhedron)))==NULL) {
    yyerror("Cannot store polyhedron");
    return(1);
  }
  op->name=NULL;
  op->object_type=POLY;
  op->contents=(void *)php;
  op->parent=NULL;
  php->llf.x= vol_infinity;
  php->llf.y= vol_infinity;
  php->llf.z= vol_infinity;
  php->urb.x= -vol_infinity;
  php->urb.y= -vol_infinity;
  php->urb.z= -vol_infinity;
}
	statement_list
{
  php->n_polys=polygon_count;
  php->n_verts=vertex_count;
  php->polygon_list=polygon_head;
  php->unique_vertex=vertex_head;
  php->vertex_list_array=vertex_list_array;
  php->element_data=NULL;
};

statement_list: statement
	| statement_list statement
;

statement: comment
	| vertex_face_block
;

comment: COMMENT
{
/*
  printf("%s\n",cval);
*/
};

vertex_face_block: vertex_list
{ 
  if ((vertex_list_array=(struct vertex_list **)malloc
       (max_vertex*sizeof(struct vertex_list *)))==NULL) {
    yyerror("Cannot store vertex array");
    return(1);
  }
  vlp=vertex_head;
  while (vlp!=NULL) {
    vertex_list_array[vlp->vertex_index-1]=vlp;
    vecp=vlp->vertex;
    if (vecp->x < php->llf.x) {
      php->llf.x=vecp->x;
    }
    if (vecp->y < php->llf.y) {
      php->llf.y=vecp->y;
    }
    if (vecp->z < php->llf.z) {
      php->llf.z=vecp->z;
    }
    if (vecp->x > php->urb.x) {
      php->urb.x=vecp->x;
    }
    if (vecp->y > php->urb.y) {
      php->urb.y=vecp->y;
    }
    if (vecp->z > php->urb.z) {
      php->urb.z=vecp->z;
    }
    vlp=vlp->next;
  }
}
	face_list
{ 
  fprintf(stderr,"\npolygon mesh:  %d vertices & %d polygons\n",
          vertex_count,polygon_count);
};

vertex_list: vertex
	| vertex_list vertex
;

vertex: VERTEX int_arg num_arg num_arg num_arg
{
  if ((vecp=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
    yyerror("Cannot store normal vector");
    return(1);
  }
  vertex_index=$<dbl>2;
  if (vertex_index>max_vertex) {
    max_vertex=vertex_index;
  }
  vecp->x=$<dbl>3;
  vecp->y=$<dbl>4;
  vecp->z=$<dbl>5;
  if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    yyerror("Cannot store vertex list");
    return(1);
  }
  vlp->vertex_count=vertex_count++;
  vlp->vertex_index=vertex_index;
  vlp->fully_outside_index=0;
  vlp->fully_inside_index=0;
  vlp->on_edge_index=0;
  vlp->vertex_status=FULLY_OUTSIDE;
  vlp->fully_outside_member=0;
  vlp->fully_inside_member=0;
  vlp->on_edge_member=0;
  vlp->vertex=vecp;
  vlp->normal=NULL;
  if (vertex_tail==NULL) {
    vertex_tail=vlp;
  }
  vertex_tail->next=vlp;
  vlp->next=NULL;
  vertex_tail=vlp;
  if (vertex_head==NULL) {
    vertex_head=vlp;
  }
/*
  printf("Vertex %d %g %g %g\n",vertex_index,vecp->x,vecp->y,vecp->z);
*/
};

face_list: face
	| face_list face
;

face: FACE int_arg int_arg int_arg int_arg
{
  polygon_index=$<dbl>2;
  vert_0=$<dbl>3;
  vert_1=$<dbl>4;
  vert_2=$<dbl>5;

  polygon_count++;
  if ((pop=(struct polygon *)malloc(sizeof(struct polygon)))==NULL) {
    yyerror("Cannot store polygon");
    return(1);
  }
  if ((plp=(struct polygon_list *)malloc(sizeof(struct polygon_list)))==NULL) {
    yyerror("Cannot store polygon list");
    return(1); 
  }
  pop->n_verts=3;
  pop->polygon_index=polygon_index;
  pop->polygon_status=FULLY_OUTSIDE;
  pop->vertex_list_array[0]=vertex_list_array[vert_0-1];
  pop->vertex_list_array[1]=vertex_list_array[vert_1-1];
  pop->vertex_list_array[2]=vertex_list_array[vert_2-1];
  pop->vertex_array[0]=vertex_list_array[vert_0-1]->vertex;
  pop->vertex_array[1]=vertex_list_array[vert_1-1]->vertex;
  pop->vertex_array[2]=vertex_list_array[vert_2-1]->vertex;
  plp->polygon=pop;  
  if (polygon_tail==NULL) {
    polygon_tail=plp;
  } 
  polygon_tail->next=plp;
  plp->next=NULL; 
  polygon_tail=plp;
  if (polygon_head==NULL) {
    polygon_head=plp;
  }
};

int_arg: INTEGER {$$=(double)ival;}
;

real_arg: REAL {$$=rval;}
;

num_arg: INTEGER {$$=(double)ival;}
	| REAL {$$=rval;}
;


%%


yyerror(s)
char *s;
{
	fprintf(stderr,"mesh_tag_region: error on line: %d of file: %s  %s\n",
	        line_num,infile,s);
	fflush(stderr);
	return(1);
}

