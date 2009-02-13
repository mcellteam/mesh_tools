%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "hash.h"
#include "mesh2gts.h"

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
struct edge *ep;
struct edge_list *edge_head,*edge_tail,*elp;
struct hash_table **edge_hashtab, *hp;
int vertex_index;
int vertex_count;
int max_vertex;
int skip_count;
int polygon_count;
int edge_count;
int polygon_value_count;
double x,y,z;
int vert_1,vert_2,vert_3,vert_4;
int v1,v2,vtmp;
int i;
char ckey[9];
char *key;

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
  edge_count=0;
  polygon_count=0;
  polygon_value_count=0;
  vlp=NULL;
  vertex_head=NULL;
  vertex_tail=NULL;
  edge_head=NULL;
  edge_tail=NULL;
  plp=NULL;
  polygon_head=NULL;
  polygon_tail=NULL;
  if ((edge_hashtab=init_hashtab(EDGE_HASHSIZE))==NULL)
  {
    yyerror("Out of memory while initializing edge hash table");
    return(1);
  }
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
  printf("%d %d %d\n",vertex_count,edge_count,polygon_count);
  vlp=vertex_head;
  while (vlp!=NULL) {
    printf("%.15g %.15g %.15g\n",vlp->vertex->x,vlp->vertex->y,vlp->vertex->z);
    vlp=vlp->next;
  }

  for(elp=edge_head;elp!=NULL;elp=elp->next) {
    ep=elp->edge;
    printf("%d %d\n",ep->vertex_index[0],ep->vertex_index[1]);
  }
  
  plp=polygon_head;
  while (plp!=NULL) {
    pop=plp->polygon;
    for (i=0;i<pop->n_verts;i++) {
      printf("%d ",pop->edges[i]->edge_index);
    }
    printf("\n");
    plp=plp->next;
  }
  fprintf(stderr,"\npolygon mesh:  %d vertices, %d edges, %d polygons\n",
          vertex_count,edge_count,polygon_count);
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
  vertex_count++;
  vlp->vertex_count=vertex_count;
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
  pop->vertex_index[0]=vertex_array[vert_1-1]->vertex_count;
  pop->vertex_index[1]=vertex_array[vert_2-1]->vertex_count;
  pop->vertex_index[2]=vertex_array[vert_3-1]->vertex_count;
  pop->polygon_id=pop->vertex_index[0]^pop->vertex_index[1]^pop->vertex_index[2];
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

  /* construct edge list */
  for (i=0;i<pop->n_verts;i++) {
    /* construct edge key */
    if (i<pop->n_verts-1) {
      v1=pop->vertex_index[i];
      v2=pop->vertex_index[i+1];
    }
    else {
      v1=pop->vertex_index[i];
      v2=pop->vertex_index[0];
    }
    if (v1>v2) {
      vtmp=v1;
      v1=v2;
      v2=vtmp;
    }
    sprintf(ckey,"%04x%04x",v1,v2);

    /* look-up edge key in edge hash table */
    hp=retrieve_key(ckey,EDGE_HASHMASK,edge_hashtab);

    /* if this edge key has not already been inserted into edge hash table */
    if (hp==NULL)
    {
      edge_count++;

      /* make permanent copy of key */
      if ((key=my_strdup(ckey))==NULL)
      {
        yyerror("Out of memory while duplicating edge key");
        return(1);
      }

      /* store edge key in edge hash table */
      if ((hp=store_key(key,EDGE_HASHMASK,edge_hashtab))==NULL)
      {
        yyerror("Out of memory while storing key in edge hash table");
        return(1);
      }

      /* allocate new edge and insert into contents of edge hash table */
      if ((ep=(struct edge *)malloc(sizeof(struct edge)))==NULL)
      {
        yyerror("Out of memory while creating edge");
        return(1);
      }
      if ((elp=(struct edge_list *)malloc(sizeof(struct edge_list)))==NULL)
      {
        yyerror("Out of memory while creating edge list");
        return(1);
      }

      hp->contents=(void *)ep;
      ep->edge_id=v1^v2;
      ep->edge_index=edge_count;
      ep->vertex_index[0]=v1;
      ep->vertex_index[1]=v2;
      elp->edge=ep;
      elp->next=NULL;
      if (edge_head==NULL) {
        edge_head=elp;
        edge_tail=elp;
      }
      edge_tail->next=elp;
      edge_tail=elp;
    }

    pop->edges[i]=hp->contents;
 
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
  pop->vertex_index[0]=vertex_array[vert_1-1]->vertex_count;
  pop->vertex_index[1]=vertex_array[vert_2-1]->vertex_count;
  pop->vertex_index[2]=vertex_array[vert_3-1]->vertex_count;
  pop->vertex_index[3]=vertex_array[vert_4-1]->vertex_count;
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
	fprintf(stderr,"mesh2gts: error on line: %d of file: %s  %s\n",
	        line_num,curr_file,s);
	fflush(stderr);
	return(1);
}

