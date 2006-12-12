%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "dx2mesh.h"

#include "lex.c"

#ifdef DEBUG
#define no_printf printf
#endif

#define VERTEX_MODE 0
#define FACE_MODE 1

extern FILE *yyin;
extern char *infile_1;
extern char *infile_2;
extern char *curr_file;
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
double n1,n2,n3;
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
int mode;
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

%token <tok> REAL INTEGER ARRAY
%type <dbl> int_arg real_arg num_arg 

%right '='
%left '+' '-'
%left '*' '/'
%left UNARYMINUS

%%

mesh_format: 
{
  mode=VERTEX_MODE;
  skip_count=0;
  vertex_count=0;
  max_vertex=0;
  vlp=NULL;
  vertex_head=NULL;
  vertex_tail=NULL;
  plp=NULL;
  polygon_head=NULL;
  polygon_tail=NULL;
/*
  fprintf(stderr,"Vertex Mode...\n");
*/
}
	array_list
{
};

array_list: array
	| array_list array
;

array: ARRAY
{
  if (mode==VERTEX_MODE) {
    if ((vecp=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
      yyerror("Cannot store normal vector");
      return(1);
    }
    vecp->x=n1;
    vecp->y=n2;
    vecp->z=n3;
    if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
      yyerror("Cannot store vertex list");
      return(1);
    }
    vlp->vertex_count=vertex_count++;
    vertex_index=vertex_count;
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
  }
  else {
    vert_1=n1;
    vert_2=n2;
    vert_3=n3;
    if ((pop=(struct polygon *)malloc(sizeof(struct polygon)))==NULL) {
      yyerror("Cannot store polygon");
      return(1);
    }
    if ((plp=(struct polygon_list *)malloc(sizeof(struct polygon_list)))==NULL) {
      yyerror("Cannot store polygon list");
      return(1);
    }
    pop->n_verts=3;
    pop->vertex_index[0]=vert_1+1;
    pop->vertex_index[1]=vert_2+1;
    pop->vertex_index[2]=vert_3+1;
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
        if (mode==VERTEX_MODE) {
          vlp=vertex_head;
          while (vlp!=NULL) {
            printf("Vertex %d  %.15g %.15g %.15g\n",vlp->vertex_index,
              vlp->vertex->x,vlp->vertex->y,vlp->vertex->z);
            vlp=vlp->next;
          }

          if ((yyin=fopen(infile_2,"r"))==NULL) {
            sprintf(err_msg,"%s %s","Cannot open mesh file:",infile_2);
            yyerror(err_msg);
            return(1);
          }
          curr_file=infile_2;
          yyless(0);
          yyclearin;
          yy_switch_to_buffer(yy_create_buffer(yyin,YY_BUF_SIZE));
          mode=FACE_MODE;
/*
          fprintf(stderr,"Face Mode...\n");
*/
          line_num=0;
          return(0);
        }
        else {
          plp=polygon_head;
          polygon_count=0;
          while (plp!=NULL) {
            polygon_count++;
            printf("Face %d  %d %d %d\n",polygon_count,
              plp->polygon->vertex_index[0],plp->polygon->vertex_index[1],
              plp->polygon->vertex_index[2]);
            plp=plp->next;
          }
          fprintf(stderr,"\npolygon mesh: %d vertices & %d polygons\n",
            vertex_count,polygon_count);
          return(1);
        }
}

yyerror(s)
char *s;
{
	fprintf(stderr,"dx2mesh: error on line: %d of file: %s  %s\n",
	        line_num,curr_file,s);
	fflush(stderr);
	return(1);
}

