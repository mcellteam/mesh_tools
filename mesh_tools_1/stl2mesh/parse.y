%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "vector.h"
#include "stl2mesh.h"

#include "lex.c"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *yyin;
extern char *infile;
extern struct object *world_obj;
extern double vol_infinity;

int ival;
double rval;
char *cval;
char *a_str,*rem_str;
char c,fmt_str[128],time_str[128];
char err_msg[128];
struct object *op;
struct polyhedron *php;
struct polygon_list *plp,*polygon_head,*polygon_tail;
struct polygon *pop;
struct nurbs *nrbp;
struct double_list *dlp_head,*dlp;
struct ctlpt_list *clp_head,*clp;
struct ctlpt *cpp;
struct vertex *vp;
struct vertex_list *ph_vlp,*ph_vertex_head,*ph_vertex_tail;
struct vertex_list *vlp,*vertex_head,*vertex_tail;
struct vector3 *vecp;
struct vector3 rand_vec;
int vertex_index;
int num_verts;

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

%token <tok> SOLID OBJ_NAME FACET NORMAL OUTER LOOP VERTEX
%token <tok> ENDLOOP ENDFACET ENDSOLID
%token <tok> INTEGER REAL
%token <tok> EOF_TOK
%type <tok> normal
%type <tok> polyhedron facet outer_loop
%type <str> object_name
%type <dbl> num_arg
%type <vec> vert

%right '='
%left '+' '-'
%left '*' '/'
%left UNARYMINUS

%%

stl_format: solid_def
	| stl_format solid_def
;


solid_def:  SOLID object_name
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
    op->name=$<str>2; 
    op->object_type=POLY;
    op->contents=NULL;

    if ((php=(struct polyhedron *)malloc(sizeof(struct polyhedron)))==NULL) {
      yyerror("Cannot store polygon object");
      return(1);
    }
    php->n_polys=0;
    php->polygon_list=NULL;
    php->n_verts=0;
    php->merged_verts=0;
    php->vertex_list=NULL;
    php->vertex_array=NULL;
    php->llf.x=vol_infinity;
    php->llf.y=vol_infinity;
    php->llf.z=vol_infinity;
    php->urb.x=-vol_infinity;
    php->urb.y=-vol_infinity;
    php->urb.z=-vol_infinity;
    plp=NULL;
    polygon_head=NULL;
    polygon_tail=NULL;
    ph_vlp=NULL;
    ph_vertex_head=NULL;
    ph_vertex_tail=NULL;
    vertex_index=0;
    rand_vec.x=0.236416584579274058342;
    rand_vec.y=0.927225593011826276779;
    rand_vec.z=0.389099507126957178116;
    normalize(&rand_vec);

    op->contents=(void *)php;
}
	polyhedron
{
    php->polygon_list=polygon_head;
    php->vertex_list=ph_vertex_head;

    if ((php->vertex_array=(struct vertex **)malloc(php->n_verts*sizeof(struct vertex *)))==NULL) {
      yyerror("Cannot store vertex array");
      return(1);
    }

    for (ph_vlp=php->vertex_list; ph_vlp!=NULL; ph_vlp=ph_vlp->next)
    {
      php->vertex_array[ph_vlp->vertex->vertex_index]=ph_vlp->vertex;
/*
      fprintf(stderr,"%d  projection = %.9g\n",ph_vlp->vertex->vertex_index,ph_vlp->vertex->projection);
*/
    }

    fprintf(stderr,"number of verts: %d  number of polygons: %d\n",php->n_verts,php->n_polys);
}
	end_of_stl_file
;


end_of_stl_file: EOF_TOK
{
  return(0);
}
	| ENDSOLID object_name EOF_TOK
{
  return(0);
};


object_name: /* empty */
{
  $$ = "";
}
	| OBJ_NAME
{
  $$ = cval;
};


polyhedron: facet
        | polyhedron facet
;


facet: FACET normal outer_loop ENDFACET
;


normal: NORMAL num_arg num_arg num_arg
;


num_arg: INTEGER {$$=(double)ival;}
        | REAL {$$=rval;}
;


outer_loop: OUTER LOOP
{
    num_verts=0;
    vlp=NULL;
    vertex_head=NULL;
    vertex_tail=NULL;
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
    pop->n_verts=0;
    plp->polygon=pop;
}
	vert_list
{
    pop->vertex_list=vertex_head;
    php->n_polys++;
}
	ENDLOOP
;


vert_list: vert_spec
	| vert_list vert_spec
;


vert_spec: VERTEX vert
{
    if ((vp=(struct vertex *)malloc(sizeof(struct vertex)))==NULL) {
      yyerror("Cannot store vertex list");
      return(1);
    }
    vp->vertex_index=vertex_index++;
    vp->merged_index=-1;
    vp->vertex=$<vec>2;
    vp->normal=NULL;
    vp->projection=dot_prod(vp->vertex,&rand_vec);

    if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
      yyerror("Cannot store vertex list");
      return(1);
    }
    vlp->vertex=vp;
    if (vertex_tail==NULL) {
      vertex_tail=vlp;
    }
    vertex_tail->next=vlp;
    vlp->next=NULL;
    vertex_tail=vlp;
    if (vertex_head==NULL) {
      vertex_head=vlp;
    }
    pop->n_verts++;

    if ((ph_vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
      yyerror("Cannot store vertex list");
      return(1);
    }
    ph_vlp->vertex=vp;
    if (ph_vertex_tail==NULL) {
      ph_vertex_tail=ph_vlp;
    }
    ph_vertex_tail->next=ph_vlp;
    ph_vlp->next=NULL;
    ph_vertex_tail=ph_vlp;
    if (ph_vertex_head==NULL) {
      ph_vertex_head=ph_vlp;
    }
    php->n_verts++;
}; 


vert: num_arg num_arg num_arg
{
    if ((vecp=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
      yyerror("Cannot store vertex");
      return(1);
    }
    vecp->x=$<dbl>1;
    vecp->y=$<dbl>2;
    vecp->z=$<dbl>3;
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

    $$ = vecp;
};



%%



yyerror(s)
char *s;
{
	fprintf(stderr,"stl2mesh: error on line: %d of file: %s  %s\n",
	        line_num,infile,s);
	fflush(stderr);
	return(1);
}
