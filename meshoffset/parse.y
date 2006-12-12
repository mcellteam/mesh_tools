%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "meshoffset.h"
#include "vector.h"

#include "lex.c"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *yyin;
extern char *infile;
extern int skip_freq;
extern struct object *world_obj;
extern double offset;

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
struct vector3 *vecp,v1,v2,*v_normal,*p_normal,y_axis;
int vertex_index;
int vertex_count;
int max_vertex;
int skip_count;
int polygon_count,polygon_index;
double x,y,z,dot;
double normal_scale;
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
  plp=polygon_head;
  while (plp!=NULL) {
    pop=plp->polygon;
    for (i=0; i<pop->n_verts; i++) {
      vlp=vertex_array[pop->vertex_index[i]-1];
      v_normal=vlp->normal;
      p_normal=pop->normal;
      normal_scale=1.0/vlp->polygon_count;
      v_normal->x=v_normal->x+(normal_scale*p_normal->x);
      v_normal->y=v_normal->y+(normal_scale*p_normal->y);
      v_normal->z=v_normal->z+(normal_scale*p_normal->z);
    }
    plp=plp->next;
  }
  vlp=vertex_head;
  while (vlp!=NULL) {
    vecp=vlp->vertex;
    v_normal=vlp->normal;
    v1.x=vecp->x+(offset*v_normal->x);
    v1.y=vecp->y+(offset*v_normal->y);
    v1.z=vecp->z+(offset*v_normal->z);
    printf("Vertex %d %.15g %.15g %.15g\n",vlp->vertex_index,v1.x,v1.y,v1.z);
    vlp=vlp->next;
  }
  plp=polygon_head;
  while (plp!=NULL) {
    pop=plp->polygon;
    printf("Face %d %d %d %d\n",pop->polygon_index,
      pop->vertex_index[0],pop->vertex_index[1],pop->vertex_index[2]);
    plp=plp->next;
  }
  
  fprintf(stderr,"\npolygon mesh:  %d vertices & %d polygons\n",
          vertex_count,polygon_count);
};

vertex_list: vertex
	| vertex_list vertex
;

vertex: VERTEX int_arg num_arg num_arg num_arg
{
  if ((vecp=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
    yyerror("Cannot store vertex");
    return(1);
  }
  if ((v_normal=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
    yyerror("Cannot store vertex normal");
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
  vlp->vertex_index=vertex_index;
  vlp->vertex_count=vertex_count++;
  vlp->polygon_count=0;
  vlp->vertex=vecp;
  vlp->normal=v_normal;
  v_normal->x=0;
  v_normal->y=0;
  v_normal->z=0;
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
  polygon_index=$<dbl>2;
  vert_0=$<dbl>3;
  vert_1=$<dbl>4;
  vert_2=$<dbl>5;
  if ((pop=(struct polygon *)malloc(sizeof(struct polygon)))==NULL) {
    yyerror("Cannot store polygon");
    return(1);
  }
  if ((plp=(struct polygon_list *)malloc(sizeof(struct polygon_list)))==NULL) {
    yyerror("Cannot store polygon");
    return(1);
  }
  pop->n_verts=3;
  pop->polygon_index=polygon_index;
  pop->polygon_count=polygon_count++;
  pop->vertex_index[0]=vert_0;
  pop->vertex_index[1]=vert_1;
  pop->vertex_index[2]=vert_2;
  vertex_array[vert_0-1]->polygon_count++;
  vertex_array[vert_1-1]->polygon_count++;
  vertex_array[vert_2-1]->polygon_count++;
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

  if ((p_normal=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
    yyerror("Cannot store polygon normal");
    return(1);
  }
  pop->normal=p_normal;
  vectorize(vertex_array[vert_0-1]->vertex,vertex_array[vert_1-1]->vertex,&v1);
  vectorize(vertex_array[vert_0-1]->vertex,vertex_array[vert_2-1]->vertex,&v2);
  cross_prod(&v1,&v2,p_normal);
  normalize(p_normal);
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
	fprintf(stderr,"meshoffset: error on line: %d of file: %s  %s\n",
	        line_num,infile,s);
	fflush(stderr);
	return(1);
}

