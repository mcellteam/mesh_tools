%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "meshfilter.h"
#include "vector.h"

#include "lex.c"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *yyin;
extern char *infile;
extern int skip_freq;
extern struct object *world_obj;
extern double cos_threshold;

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
struct vector3 *vecp,v1,v2,normal,y_axis;
int vertex_index;
int vertex_count;
int max_vertex;
int skip_count;
int polygon_count,polygon_index;
double x,y,z,dot;
int vert_1,vert_2,vert_3;
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
  y_axis.x=0;
  y_axis.y=1;
  y_axis.z=0;
}
	statement_list
;

statement_list: statement
	| statement_list statement
;

statement: comment
	| vertex_face_block
;

comment: COMMENT
{
  printf("%s\n",cval);
};

vertex_face_block: vertex_list
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
  vertex_count++;
  if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    yyerror("Cannot store vertex list");
    return(1);
  }
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
  printf("Vertex %d %g %g %g\n",vertex_index,vecp->x,vecp->y,vecp->z);
};

face_list: face
	| face_list face
;

face: FACE int_arg int_arg int_arg int_arg
{
  polygon_index=$<dbl>2;
  vert_1=$<dbl>3;
  vert_2=$<dbl>4;
  vert_3=$<dbl>5;
  vectorize(vertex_array[vert_1-1]->vertex,vertex_array[vert_2-1]->vertex,&v1);
  vectorize(vertex_array[vert_1-1]->vertex,vertex_array[vert_3-1]->vertex,&v2);
  cross_prod(&v1,&v2,&normal);
  normalize(&normal);
  dot=fabs(dot_prod(&normal,&y_axis));
  if (dot>cos_threshold) {
    fprintf(stderr,"Y-Axis oriented triangle found.  Face Index = %d  Cosine = %g\n",polygon_index,dot);
  }
  else {
    printf("Face %d %d %d %d\n",polygon_index,vert_1,vert_2,vert_3);
  }
  polygon_count++;
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
	fprintf(stderr,"meshfilter: error on line: %d of file: %s  %s\n",
	        line_num,infile,s);
	fflush(stderr);
	return(1);
}

