%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "strfunc.h"
#include "netgen2mesh.h"

#include "netgenlex.c"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *netgenin;
extern char *infile;
extern char *curr_file;
extern int skip_freq;
extern struct object *world_obj;
extern struct vertex_list **vertex_array;

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
struct polygon *pop;
struct polygon_list *plp,*polygon_head,*polygon_tail;
struct tet *tetp;
struct tet_list *tlp,*tet_head,*tet_tail;
struct nurbs *nrbp;
struct double_list *dlp_head,*dlp;
struct ctlpt_list *clp_head,*clp;
struct ctlpt *cpp;
struct vertex_list *vlp,*vertex_head,*vertex_tail;
struct vector3 *vecp;
int skip_count;
int vertex_count;
int vertex_index;
int polygon_count;
int tet_count;
int output_count;
double x,y,z;
int vert_1,vert_2,vert_3,vert_4;
int i;


%}

%union {
int tok;
char *str;
double dbl;
struct vector3 *vec;
struct object *obj;
} 

%token <tok> DIMENSION EDGE_SEGMENTS MESH3D POINTS SURFACE_ELEMENTS
%token <tok> VOLUME_ELEMENTS NEWLINE REAL INTEGER 
%type <dbl> int_arg real_arg num_arg 

%right '='
%left '+' '-'
%left '*' '/'
%left UNARYMINUS

%%

netgen_format: 
{
  skip_count=0;
  vertex_count=0;
  vertex_index=0;
  polygon_count=0;
  tet_count=0;
  vlp=NULL;
  vertex_head=NULL;
  vertex_tail=NULL;
  plp=NULL;
  polygon_head=NULL;
  polygon_tail=NULL;
  tlp=NULL;
  tet_head=NULL;
  tet_tail=NULL;
}
	netgen_header
        surface_elements_block
        volume_elements_block
        edge_segments_block
        points_block
{ 
  if ((vertex_array=(struct vertex_list **)malloc
       (vertex_count*sizeof(struct vertex_list *)))==NULL) {
    netgenerror("Cannot store vertex array");
    return(1);
  }
  vlp=vertex_head;
  while (vlp!=NULL) {
    vertex_array[vlp->vertex_index-1]=vlp;
    vlp=vlp->next;
  }

  vlp=vertex_head;
  while (vlp!=NULL) {
    printf("Vertex %d  %.15g %.15g %.15g\n",vlp->vertex_index,
      vlp->vertex->x,vlp->vertex->y,vlp->vertex->z);
    vlp=vlp->next;
  }
  plp=polygon_head;
  polygon_count=0;
  while (plp!=NULL) {
    polygon_count++;
    printf("Face %d  %d %d %d\n",polygon_count,
      plp->polygon->vertex_index[0],plp->polygon->vertex_index[1],
      plp->polygon->vertex_index[2]);
    plp=plp->next;
  }

  fprintf(stderr,"\npolygon mesh:  %d vertices & %d polygons\n",
          vertex_count,polygon_count);
};

newline_list: NEWLINE
	| newline_list NEWLINE
;

netgen_header: MESH3D newline_list
	DIMENSION newline_list
	int_arg newline_list
;

surface_elements_block:  SURFACE_ELEMENTS newline_list
	int_arg newline_list
	surface_element_list
;

surface_element_list: surface_element
	| surface_element_list surface_element
;

surface_element:	int_arg int_arg int_arg int_arg int_arg int_arg int_arg int_arg int_arg int_arg int_arg newline_list
{
  vert_1=$<dbl>6;
  vert_2=$<dbl>7;
  vert_3=$<dbl>8;
  polygon_count++;
  if ((pop=(struct polygon *)malloc(sizeof(struct polygon)))==NULL) {
    netgenerror("Cannot store polygon");
    return(1);
  }
  if ((plp=(struct polygon_list *)malloc(sizeof(struct polygon_list)))==NULL) {
    netgenerror("Cannot store polygon list");
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

volume_elements_block:  VOLUME_ELEMENTS newline_list
	int_arg newline_list
	volume_element_list
;

volume_element_list: /* empty */
	| volume_element_list volume_element
;

volume_element:	int_arg int_arg int_arg int_arg int_arg int_arg newline_list
{
  vert_1=$<dbl>3;
  vert_2=$<dbl>4;
  vert_3=$<dbl>5;
  vert_4=$<dbl>6;
  tet_count++;
  if ((tetp=(struct tet *)malloc(sizeof(struct tet)))==NULL) {
    netgenerror("Cannot store tet");
    return(1);
  }
  if ((tlp=(struct tet_list *)malloc(sizeof(struct tet_list)))==NULL) {
    netgenerror("Cannot store tet list");
    return(1);
  }
  tetp->vertex_index[0]=vert_1;
  tetp->vertex_index[1]=vert_2;
  tetp->vertex_index[2]=vert_3;
  tetp->vertex_index[3]=vert_4;
  tlp->tet=tetp;
  if (tet_tail==NULL) {
    tet_tail=tlp;
  }
  tet_tail->next=tlp;
  tlp->next=NULL;
  tet_tail=tlp;
  if (tet_head==NULL) {
    tet_head=tlp;
  }
};

edge_segments_block:  EDGE_SEGMENTS newline_list
	int_arg newline_list
	edge_segment_list
;

edge_segment_list: edge_segment
	| edge_segment_list edge_segment
;

//edge_segment:	num_arg num_arg num_arg num_arg num_arg num_arg newline_list
edge_segment:	num_arg_list newline_list
;

num_arg_list:  num_arg
	| num_arg_list num_arg
;

points_block:  POINTS newline_list
	int_arg newline_list
	point_list
;

point_list: point
	| point_list point
;

point: num_arg num_arg num_arg newline_list
{
  if ((vecp=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
    netgenerror("Cannot store normal vector");
    return(1);
  }
  vecp->x=$<dbl>1;
  vecp->y=$<dbl>2;
  vecp->z=$<dbl>3;
  if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    netgenerror("Cannot store vertex list");
    return(1);
  }
  vlp->vertex_count=vertex_count++;
  vertex_index++;
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

int_arg: INTEGER {$$=(double)ival;}
;

real_arg: REAL {$$=rval;}
;

num_arg: INTEGER {$$=(double)ival;}
	| REAL {$$=rval;}
;

%%


netgenerror(char *s)
{
	fprintf(stderr,"netgen2mesh: error on line: %d of file: %s  %s\n",
	        line_num,curr_file,s);
	fflush(stderr);
	return(1);
}

