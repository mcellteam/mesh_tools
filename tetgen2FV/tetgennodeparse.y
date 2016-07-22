%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "strfunc.h"
#include "tetgen2FV.h"
#include "tetgennodeparse.h"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *tetgennodein;
extern char *infile_1;
extern char *infile_2;
extern char *outfile_name;
extern char *curr_file;
extern int line_num;
extern struct tet_mesh *tet_mesh;

FILE *outfile;
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
struct nurbs *nrbp;
struct double_list *dlp_head,*dlp;
struct ctlpt_list *clp_head,*clp;
struct ctlpt *cpp;
struct node *np;
struct node **node_array;
struct node_list *nlp,*node_head,*node_tail;
u_int node_count;
int node_index;
double x,y,z;
int vert_1,vert_2,vert_3,vert_4;
int i;


%}

%union {
int tok;
char *str;
double dbl;
struct object *obj;
} 

%{
  #include "tetgennodelex.c"
%}

%token <tok> REAL INTEGER NEWLINE
%type <dbl> int_arg real_arg num_arg 

%right '='
%left '+' '-'
%left '*' '/'
%left UNARYMINUS


%%


tetgennode_format: 
{
  fprintf(stderr,"tetgen2FV: parsing tetgen node file %s...\n",curr_file);
  node_count=0;
  node_index=0;
  nlp=NULL;
  node_head=NULL;
  node_tail=NULL;
}
	node_header
	node_list
{ 
  if ((node_array=(struct node **)malloc
       (node_count*sizeof(struct node *)))==NULL) {
    tetgennodeerror("Out of memory while storing node array");
    return(1);
  }
  nlp=node_head;
  while (nlp!=NULL) {
    node_array[nlp->node->node_index]=nlp->node;
    nlp=nlp->next;
  }

  tet_mesh->nnodes=node_count;
  tet_mesh->node_head=node_head;
  tet_mesh->nodes=node_array;
  fprintf(stderr,"          ... %d nodes found\n",node_count);

};


newline_list: NEWLINE
        | newline_list NEWLINE
;


node_header: int_arg int_arg int_arg int_arg newline_list
;


node_list: node 
	| node_list node
;


node: int_arg num_arg num_arg num_arg newline_list
{

  if ((np=(struct node *)malloc(sizeof(struct node)))==NULL)
  {
    tetgennodeerror("Out of memory while storing node");
    return(1);
  }
  if ((nlp=(struct node_list *)malloc(sizeof(struct node_list)))==NULL)
  {
    tetgennodeerror("Out of memory while storing node list");
    return(1);
  }

  np->x=$<dbl>2;
  np->y=$<dbl>3;
  np->z=$<dbl>4;
  np->srfnum=0;
  np->node_index=node_index;
  np->voronoi_index=0;
  np->node_type=INTERIOR;
  np->nedges=0;
  np->ntets=0;
  np->nfacets=0;
  np->nboundaryfaces=0;
  np->ninterfaceneighbors=0;
  np->material_count=0;
  np->material_head=NULL;
  np->node_material=NULL;
  np->volume_element=NULL;
  np->edge_head=NULL;
  np->boundary_head=NULL;

  nlp->node=np;
  if (node_tail==NULL)
  {
    node_tail=nlp;
  }
  node_tail->next=nlp;
  nlp->next=NULL;
  node_tail=nlp;
  if (node_head==NULL)
  {
    node_head=nlp;
  }

  node_count++;
  node_index++;

}
	| int_arg num_arg num_arg num_arg int_arg newline_list
{

  if ((np=(struct node *)malloc(sizeof(struct node)))==NULL)
  {
    tetgennodeerror("Out of memory while storing node");
    return(1);
  }
  if ((nlp=(struct node_list *)malloc(sizeof(struct node_list)))==NULL)
  {
    tetgennodeerror("Out of memory while storing node list");
    return(1);
  }

  np->x=$<dbl>2;
  np->y=$<dbl>3;
  np->z=$<dbl>4;
  np->srfnum=$<dbl>5;
  np->node_index=node_index;
  if (np->srfnum==0)
  {
    np->node_type=INTERIOR;
  }
  else
  {
    np->node_type=BOUNDARY;
    tet_mesh->nboundarynodes++;
  }
  np->nedges=0;
  np->ntets=0;
  np->nfacets=0;
  np->nboundaryfaces=0;
  np->ninterfaceneighbors=0;
  np->material_count=0;
  np->material_head=NULL;
  np->volume_element=NULL;
  np->edge_head=NULL;
  np->boundary_head=NULL;

  nlp->node=np;
  if (node_tail==NULL)
  {
    node_tail=nlp;
  }
  node_tail->next=nlp;
  nlp->next=NULL;
  node_tail=nlp;
  if (node_head==NULL)
  {
    node_head=nlp;
  }

  node_count++;
  node_index++;

};


int_arg: INTEGER {$$=(double)ival;}
;


real_arg: REAL {$$=rval;}
;


num_arg: INTEGER {$$=(double)ival;}
	| REAL {$$=rval;}
;


%%


#undef tetgennodewrap

int tetgennodewrap()
{
        return(1);
} 


int tetgennodeerror(char *s)
{
	fprintf(stderr,"tetgen2FV: error on line: %d of file: %s  %s\n",
	        line_num,curr_file,s);
	fflush(stderr);
	return(1);
}

