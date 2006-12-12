%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "vtk2rib.h"

#include "lex.c"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *yyin;
extern char *infile;
extern char *object_name;
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
struct vertex_list *vlp,*vertex_head,*vertex_tail;
struct vector3 *vecp;
int vertex_index;
int vertex_count;
int max_vertex;
int skip_count;
int polygon_count;
int output_normals;
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

%token <tok> REAL INTEGER FLOAT POINTS POINT_DATA POLYGONS NORMALS
%token <tok> NORMALS_LOWER
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
  output_normals=0;
}
	vertex_block
	face_block
{
  vlp=vertex_head;
}
        normals_block
{ 
  printf("PointsPolygons\n");
  printf("    [\n");
  plp=polygon_head;
  polygon_count=0; 
  while (plp!=NULL) {
    polygon_count++;
    printf("\t%d\n",3);
    plp=plp->next;
  }
  printf("    ]\n");

  printf("    [\n");
  plp=polygon_head;
  polygon_count=0; 
  while (plp!=NULL) {
    polygon_count++;
    printf("\t%d %d %d\n",plp->polygon->vertex_index[0],
      plp->polygon->vertex_index[1],plp->polygon->vertex_index[2]);
    plp=plp->next;
  }
  printf("    ]\n");

  printf("    \"P\"  [\n");
  vlp=vertex_head;
  while (vlp!=NULL) {
    printf("\t%.9g %.9g %.9g\n",
      vlp->vertex->x,vlp->vertex->y,vlp->vertex->z);
    vlp=vlp->next;
  }
  printf("         ]\n");

  if (output_normals) {
    printf("    \"N\"  [\n");
    vlp=vertex_head;
    while (vlp!=NULL) {
      if (vlp->normal!=NULL) {
        printf("\t%.9g %.9g %.9g\n",
          vlp->normal->x,vlp->normal->y,vlp->normal->z);
      }
      else {
        yyerror("Number of normals is less than the number of vertices");
        return(1);
      }
      vlp=vlp->next;
    }
    printf("         ]\n");
  }

  fprintf(stderr,"\npolygon mesh:  %d vertices & %d polygons\n",
          vertex_count,polygon_count);
};

vertex_block: POINTS int_arg FLOAT vertex_list
{
};

vertex_list: vertex
	| vertex_list vertex
;

vertex: num_arg num_arg num_arg
{
  if ((vecp=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
    yyerror("Cannot store vertex");
    return(1);
  }
  vecp->x=$<dbl>1;
  vecp->y=$<dbl>2;
  vecp->z=$<dbl>3;
  if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    yyerror("Cannot store vertex list");
    return(1);
  }
  vlp->vertex_count=vertex_count++;
  vertex_index=vertex_count;
  if (vertex_index>max_vertex) {
    max_vertex=vertex_index;
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
};

face_block: POLYGONS int_arg int_arg face_list
{
};

face_list: face
	| face_list face
;

face: int_arg int_arg int_arg int_arg
{
  vert_1=$<dbl>2;
  vert_2=$<dbl>3;
  vert_3=$<dbl>4;
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

normals_block: /* empty */
	| POINT_DATA int_arg NORMALS NORMALS_LOWER FLOAT normals_list
{
  output_normals = 1;
};

normals_list: normal
	| normals_list normal
;

normal: num_arg num_arg num_arg
{
  if (vlp==NULL) {
    yyerror("Number of normals exceeds number of vertices");
    return(1);
  }
  if ((vecp=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
    yyerror("Cannot store normal vector");
    return(1);
  }
  vecp->x=$<dbl>1;
  vecp->y=$<dbl>2;
  vecp->z=$<dbl>3;
  vlp->normal=vecp;
  vlp=vlp->next;
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
	fprintf(stderr,"vtk2rib: error on line: %d of file: %s  %s\n",
	        line_num,infile,s);
	fflush(stderr);
	return(1);
}

