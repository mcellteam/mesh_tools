%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "vector.h"
#include "mesh2stl.h"

#include "lex.c"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *yyin;
extern char *infile;
extern char *curr_file;
extern int skip_freq;
extern struct object *world_obj;

int ival;
double rval;
char *cval;
char *var;
char *a_str,*rem_str;
char c,fmt_str[128],time_str[128];
char err_msg[128];
struct object *op,*curr_op;
struct meta_object *mop,*curr_mop,*wmop;
struct object_list *olp;
struct polyhedron *php;
struct polygon_list *plp,*polygon_head,*polygon_tail;
struct polygon *pop;
struct nurbs *nrbp;
struct double_list *dlp_head,*dlp;
struct ctlpt_list *clp_head,*clp;
struct ctlpt *cpp;
struct vertex_list *vlp,*vertex_head,*vertex_tail,**vertex_array;
struct vector3 *vecp;
struct vector3 p1,p2,v1,v2,normal;
int vertex_index;
int vertex_count;
int max_vertex;
int skip_count;
int polygon_count;
int polygon_value_count;
double x,y,z;
int vert_1,vert_2,vert_3,vert_4;
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

%token <tok> REAL INTEGER VERTEX FACE
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
  polygon_value_count=0;
  vlp=NULL;
  vertex_head=NULL;
  vertex_tail=NULL;
  plp=NULL;
  polygon_head=NULL;
  polygon_tail=NULL;
}
	vertex_list
{ 
  if ((vertex_array=(struct vertex_list **)malloc
       (max_vertex*sizeof(struct vertex_list *)))==NULL) {
    yyerror("Cannot store vertex array");
    return(1);
  }
  vlp=vertex_head;
  while (vlp!=NULL) {
    vertex_array[vlp->vertex_index-1]=vlp;
    vlp=vlp->next;
  }
}
	face_list
{ 
  printf("solid\n");
  plp=polygon_head;
  while (plp!=NULL) {
    vectorize(vertex_array[plp->polygon->vertex_index[0]-1]->vertex,vertex_array[plp->polygon->vertex_index[1]-1]->vertex,&v1);
    vectorize(vertex_array[plp->polygon->vertex_index[0]-1]->vertex,vertex_array[plp->polygon->vertex_index[2]-1]->vertex,&v2);
    cross_prod(&v1,&v2,&normal);
    if (vect_length(&normal)==0)
    {
      normal.x=0;
      normal.y=0;
      normal.z=0;
    }
    else 
    {
      normalize(&normal);
    }
    printf("  facet normal %.15g %.15g %.15g\n",normal.x,normal.y,normal.z);
    printf("    outer loop\n");
    printf("      vertex %.15g %.15g %.15g\n",
      vertex_array[plp->polygon->vertex_index[0]-1]->vertex->x,
      vertex_array[plp->polygon->vertex_index[0]-1]->vertex->y,
      vertex_array[plp->polygon->vertex_index[0]-1]->vertex->z);
    printf("      vertex %.15g %.15g %.15g\n",
      vertex_array[plp->polygon->vertex_index[1]-1]->vertex->x,
      vertex_array[plp->polygon->vertex_index[1]-1]->vertex->y,
      vertex_array[plp->polygon->vertex_index[1]-1]->vertex->z);
    printf("      vertex %.15g %.15g %.15g\n",
      vertex_array[plp->polygon->vertex_index[2]-1]->vertex->x,
      vertex_array[plp->polygon->vertex_index[2]-1]->vertex->y,
      vertex_array[plp->polygon->vertex_index[2]-1]->vertex->z);
    printf("    endloop\n");
    printf("  endfacet\n");
    plp=plp->next;
  }
  printf("endsolid\n");
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
};

face_list: face
	| face_list face
;

face: FACE int_arg int_arg int_arg int_arg
{
  vert_1=$<dbl>3;
  vert_2=$<dbl>4;
  vert_3=$<dbl>5;
  polygon_count++;
  polygon_value_count=polygon_value_count+4;
  if ((pop=(struct polygon *)malloc(sizeof(struct polygon)))==NULL) {
    yyerror("Cannot store polygon");
    return(1);
  }
  if ((plp=(struct polygon_list *)malloc(sizeof(struct polygon_list)))==NULL) {
    yyerror("Cannot store polygon list");
    return(1);
  }
  pop->n_verts=3;
  pop->vertex_index[0]=vertex_array[vert_1-1]->vertex_index;
  pop->vertex_index[1]=vertex_array[vert_2-1]->vertex_index;
  pop->vertex_index[2]=vertex_array[vert_3-1]->vertex_index;
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
}
	| FACE int_arg int_arg int_arg int_arg int_arg
{
  vert_1=$<dbl>3;
  vert_2=$<dbl>4;
  vert_3=$<dbl>5;
  vert_4=$<dbl>6;
  polygon_count++;
  polygon_value_count=polygon_value_count+5;
  if ((pop=(struct polygon *)malloc(sizeof(struct polygon)))==NULL) {
    yyerror("Cannot store polygon");
    return(1);
  }
  if ((plp=(struct polygon_list *)malloc(sizeof(struct polygon_list)))==NULL) {
    yyerror("Cannot store polygon list");
    return(1);
  }
  pop->n_verts=4;
  pop->vertex_index[0]=vertex_array[vert_1-1]->vertex_index;
  pop->vertex_index[1]=vertex_array[vert_2-1]->vertex_index;
  pop->vertex_index[2]=vertex_array[vert_3-1]->vertex_index;
  pop->vertex_index[3]=vertex_array[vert_4-1]->vertex_index;
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
	fprintf(stderr,"mesh2stl: error on line: %d of file: %s  %s\n",
	        line_num,curr_file,s);
	fflush(stderr);
	return(1);
}

