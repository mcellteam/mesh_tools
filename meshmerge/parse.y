%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "meshmerge.h"

#include "lex.c"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *yyin;
extern char *infile_1;
extern char *infile_2;
extern char *curr_file;
extern int parse_seq;
extern int skip_freq;
extern struct object *world_obj;
extern struct vector3 translate;

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
struct vertex_list *vlp,*vertex_head,*vertex_tail;
struct vector3 *vecp;
int vertex_index;
int vertex_count;
int max_vertex;
int skip_count;
int polygon_count;
double x,y,z;
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
  vlp=NULL;
  vertex_head=NULL;
  vertex_tail=NULL;
  plp=NULL;
  polygon_head=NULL;
  polygon_tail=NULL;
}
	vertex_list
	face_list
	vertex_list
	face_list
{ 
  vlp=vertex_head;
  while (vlp!=NULL) {
    printf("Vertex %d %.15g %.15g %.15g\n",
      vlp->vertex_index,vlp->vertex->x,vlp->vertex->y,vlp->vertex->z);
    vlp=vlp->next;
  }
  plp=polygon_head;
  polygon_count=0;
  while (plp!=NULL) {
    polygon_count++;
    printf("Face %d %d %d %d\n",polygon_count,plp->polygon->vertex_index[0],
      plp->polygon->vertex_index[1],plp->polygon->vertex_index[2]);
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
    yyerror("Cannot store normal vector");
    return(1);
  }
  vertex_index=$<dbl>2;
  if (parse_seq==0) {
    if (vertex_index>max_vertex) {
      max_vertex=vertex_index;
    }
  }
  if (parse_seq==0) {
    vecp->x=$<dbl>3;
    vecp->y=$<dbl>4;
    vecp->z=$<dbl>5;
  }
  else {
    vecp->x=$<dbl>3+translate.x;
    vecp->y=$<dbl>4+translate.y;
    vecp->z=$<dbl>5+translate.z;
  }
  if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    yyerror("Cannot store vertex list");
    return(1);
  }
  vlp->vertex_count=vertex_count++;
  if (parse_seq==0) {
    vlp->vertex_index=vertex_index;
  }
  else {
    vlp->vertex_index=vertex_index+max_vertex;
  }
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
  if (parse_seq==0) {
    vert_1=$<dbl>3;
    vert_2=$<dbl>4;
    vert_3=$<dbl>5;
  }
  else {
    vert_1=$<dbl>3+max_vertex;
    vert_2=$<dbl>4+max_vertex;
    vert_3=$<dbl>5+max_vertex;
  }
  if ((pop=(struct polygon *)malloc(sizeof(struct polygon)))==NULL) {
    yyerror("Cannot store polygon");
    return(1);
  }
  if ((plp=(struct polygon_list *)malloc(sizeof(struct polygon_list)))==NULL) {
    yyerror("Cannot store polygon list");
    return(1);
  }
  pop->n_verts=3;
  pop->vertex_index[0]=vert_1;
  pop->vertex_index[1]=vert_2;
  pop->vertex_index[2]=vert_3;
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

yywrap()
{
        if (parse_seq==0) {
          if ((yyin=fopen(infile_2,"r"))==NULL) {
            sprintf(err_msg,"%s %s","Cannot open mesh file:",infile_2);
            yyerror(err_msg);
            return(1);
          }
          curr_file=infile_2;
          yyless(0);
          yyclearin;
          yy_switch_to_buffer(yy_create_buffer(yyin,YY_BUF_SIZE));
          parse_seq=1;
          line_num=0;
          return(0);
        }
        else {
          return(1);
        }
}

yyerror(s)
char *s;
{
	fprintf(stderr,"meshmerge: error on line: %d of file: %s  %s\n",
	        line_num,curr_file,s);
	fflush(stderr);
	return(1);
}

