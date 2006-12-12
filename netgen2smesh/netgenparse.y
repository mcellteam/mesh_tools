%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "strfunc.h"
#include "netgen2smesh.h"

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
struct material_list *mlp;
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
int multimat_tet_count;
int multimat_max;
int mat_found;
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

%token <tok> CSG_SURFACES CYLINDER PLANE SPHERE END_MESH
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
  multimat_tet_count=0;
  multimat_max=0;
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
        end_mesh_cmd
        csg_surfaces_block
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

  plp=polygon_head;
  while (plp!=NULL) {
    pop=plp->polygon; 
    for (i=0;i<3;i++) {
      vertex_array[pop->vertex_index[i]-1]->srfnum=pop->srfnum;
      if (pop->domout==0) {
        vertex_array[pop->vertex_index[i]-1]->vertex_type=BOUNDARY;
      }
      else {
        vertex_array[pop->vertex_index[i]-1]->vertex_type=PARENT;
      }
    }
    plp=plp->next;
  }

  tlp=tet_head;
  while (tlp!=NULL) {
    for (i=0;i<4;i++) {
      mat_found=0;
      mlp=vertex_array[tlp->tet->vertex_index[i]-1]->material;
      while (mlp!=NULL) {
        if (tlp->tet->matnum==mlp->matnum) {  
          mat_found=1;  
        }
        mlp=mlp->next;
      }
      if (!mat_found) {
        if ((mlp=(struct material_list *)malloc
          (sizeof(struct material_list)))==NULL) {
          netgenerror("Cannot store vertex material");
          return(1);
        }
        mlp->matnum=tlp->tet->matnum;
        mlp->next=vertex_array[tlp->tet->vertex_index[i]-1]->material;
        vertex_array[tlp->tet->vertex_index[i]-1]->material=mlp;
        vertex_array[tlp->tet->vertex_index[i]-1]->material_count++;
        tlp->tet->material_count++;
        if (vertex_array[tlp->tet->vertex_index[i]-1]->material_count > multimat_max) {
          multimat_max=vertex_array[tlp->tet->vertex_index[i]-1]->material_count;
        }
      }
    }
    if (tlp->tet->material_count > 1){
      multimat_tet_count++;
    }
    tlp=tlp->next;
  }

  printf("# TetGen Surface Mesh in .smesh format\n\n");
  printf("# part 1 - node list\n");
  printf("#   n nodes, 3 dimensions, 0 attributes, boundary marks on\n");
  printf("%d 3 0 1\n\n",vertex_count);

  vlp=vertex_head;
  for (i=0; i<vertex_count; i++)
  {
    vlp=vertex_array[i];
    printf("%d %.15g %.15g %.15g %d\n",
      i+1,vlp->vertex->x,vlp->vertex->y,vlp->vertex->z,vlp->srfnum);
/*
    if (vlp->material_count==1)
    {
      printf("%d %.15g %.15g %.15g %d\n",
        i+1,vlp->vertex->x,vlp->vertex->y,vlp->vertex->z,0);
    }
    else
    {
      printf("%d %.15g %.15g %.15g %d\n",
        i+1,vlp->vertex->x,vlp->vertex->y,vlp->vertex->z,1);
    }
*/
  }

  printf("\n");
  printf("# part 2 - facet list\n");
  printf("#   n facets, boundary marks on\n");
  printf("%d 1\n\n",polygon_count);

  plp=polygon_head;
  while (plp!=NULL) {
    pop=plp->polygon; 
    printf("3 %d %d %d %d\n",
      pop->vertex_index[0],
      pop->vertex_index[1],
      pop->vertex_index[2],
      pop->srfnum);
    plp=plp->next;
  }

  printf("\n");
  printf("# part 3 - hole list\n");
  printf("0\n\n");
  printf("# part 4 - region list\n");
  printf("0\n\n");

  fprintf(stderr,"\nvolume mesh:  %d vertices,  %d polygons,  %d tets\n",vertex_count,polygon_count,tet_count);
  fprintf(stderr,"              %d multi-material tets,  %d maximum materials per vertex\n",multimat_tet_count,multimat_max);
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
  polygon_count++;
  if ((pop=(struct polygon *)malloc(sizeof(struct polygon)))==NULL) {
    netgenerror("Cannot store polygon");
    return(1);
  }
  if ((plp=(struct polygon_list *)malloc(sizeof(struct polygon_list)))==NULL) {
    netgenerror("Cannot store polygon list");
    return(1);
  }
  pop->srfnum=$<dbl>1;
  pop->bcnum=$<dbl>2;
  pop->domin=$<dbl>3;
  pop->domout=$<dbl>4;
  pop->n_verts=$<dbl>5;
  if (pop->n_verts != 3) {
    netgenerror("Non-triangular surface element found: %d",polygon_count);
    return(1);
  }
  vert_1=$<dbl>6;
  vert_2=$<dbl>7;
  vert_3=$<dbl>8;
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
  tetp->matnum=$<dbl>1;
  tetp->material_count=0;
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

edge_segment:	num_arg num_arg num_arg num_arg num_arg num_arg newline_list
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
  vlp->vertex_type=INTERIOR;
  vlp->srfnum=0;
  vlp->material_count=0;
  vlp->material=NULL;
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

end_mesh_cmd:  END_MESH newline_list
;

csg_surfaces_block: CSG_SURFACES int_arg newline_list
                    csg_surfaces_list
;

csg_surfaces_list: csg_surface
                  | csg_surfaces_list csg_surface
;

csg_surface: cylinder
             | sphere
             | plane
;

cylinder: CYLINDER int_arg newline_list
          num_arg num_arg num_arg num_arg num_arg num_arg num_arg newline_list
;

sphere: SPHERE int_arg newline_list
          num_arg num_arg num_arg num_arg num_arg num_arg num_arg newline_list
;

plane: PLANE int_arg newline_list
          num_arg num_arg num_arg num_arg num_arg num_arg newline_list
;

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
	fprintf(stderr,"netgen2smesh: error on line: %d of file: %s  %s\n",
	        line_num,curr_file,s);
	fflush(stderr);
	return(1);
}

