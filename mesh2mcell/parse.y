%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "mesh2mcell.h"

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
struct vertex_list *vlp,*vertex_head,*vertex_tail,**vertex_array;
struct vertex_index_list *vilp,*vertex_index_head,*vertex_index_tail;
struct vector3 *vecp;
int n_verts;
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


%{
  #include "lex.flex.c"
%}

%output="parse.bison.c"


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
  printf("%s POLYGON_LIST {\n",object_name);
  printf("  VERTEX_LIST {\n");
  vlp=vertex_head;
  while (vlp!=NULL) {
    printf("    [ %.9g, %.9g, %.9g ]\n",vlp->vertex->x,vlp->vertex->y,vlp->vertex->z);
    vlp=vlp->next;
  }
  printf("  }\n");
  printf("  ELEMENT_CONNECTIONS {\n");
  plp=polygon_head;
  while (plp!=NULL) {
    n_verts=plp->polygon->n_verts;
    vilp=plp->polygon->vertex_index;
    printf("    [ ");
    for (i=0;i<n_verts;i++) {
      if (i>0) {
        printf(", ");
      }
      printf("%d",vilp->vertex_index);
      vilp=vilp->next; 
    } 
    printf(" ]\n");
    plp=plp->next;
  }
  printf("  }\n");
/*
  printf("  FULLY_CLOSED = NO\n");
*/
  printf("}\n");
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

face: FACE int_arg
{
  polygon_count++;
  if ((pop=(struct polygon *)malloc(sizeof(struct polygon)))==NULL) {
    yyerror("Cannot store polygon");
    return(1);
  }
  if ((plp=(struct polygon_list *)malloc(sizeof(struct polygon_list)))==NULL) {
    yyerror("Cannot store polygon list");
    return(1);
  }
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

  pop->n_verts=0;
  vertex_index_head=NULL;
  vertex_index_tail=NULL;

}
  list_vertex_indices
{
  pop->vertex_index=vertex_index_head;
  if (pop->n_verts<3) {
    printf("n_verts = %d\n",n_verts);
    yyerror("Less than 3 vertices for Face");
    return(1);
  }
};

list_vertex_indices: list_vertex_indices int_arg
{
  vertex_index=$<dbl>2;

  if ((vilp=(struct vertex_index_list *)malloc
       (sizeof(struct vertex_index_list)))==NULL) {
    yyerror("Cannot store vertex index list");
    return(1);
  }

  vilp->vertex_index=vertex_array[vertex_index-1]->vertex_count;
  pop->n_verts++;
  
  if (vertex_index_tail==NULL) {
    vertex_index_tail=vilp;
  }
  vertex_index_tail->next=vilp;
  vilp->next=NULL;
  vertex_index_tail=vilp;
  if (vertex_index_head==NULL) {
    vertex_index_head=vilp;
  }
}
	| int_arg
{
  vertex_index=$<dbl>1;

  if ((vilp=(struct vertex_index_list *)malloc
       (sizeof(struct vertex_index_list)))==NULL) {
    yyerror("Cannot store vertex index list");
    return(1);
  }

  vilp->vertex_index=vertex_array[vertex_index-1]->vertex_count;
  pop->n_verts++;
  
  if (vertex_index_tail==NULL) {
    vertex_index_tail=vilp;
  }
  vertex_index_tail->next=vilp;
  vilp->next=NULL;
  vertex_index_tail=vilp;
  if (vertex_index_head==NULL) {
    vertex_index_head=vilp;
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


#undef yywrap

yywrap()
{
        return(1);
}


yyerror(s)
char *s;
{
	fprintf(stderr,"mesh2mcell: error on line: %d of file: %s  %s\n",
	        line_num,infile,s);
	fflush(stderr);
	return(1);
}

