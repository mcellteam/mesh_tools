%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "strfunc.h"
#include "tetgen2FV.h"
#include "tetgenfaceparse.h"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *tetgenfacein;
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
struct polygon_list *plp,*polygon_head,*polygon_tail;
struct polygon *pop;
struct nurbs *nrbp;
struct double_list *dlp_head,*dlp;
struct ctlpt_list *clp_head,*clp;
struct ctlpt *cpp;
int polygon_count;
double x,y,z;
int node_1,node_2,node_3,node_4;
u_int srfnum;
int i;


%}

%union {
int tok;
char *str;
double dbl;
struct object *obj;
} 

%{
  #include "tetgenfacelex.c"
%}


%token <tok> REAL INTEGER NEWLINE
%type <dbl> int_arg real_arg num_arg 

%right '='
%left '+' '-'
%left '*' '/'
%left UNARYMINUS


%%


tetgenface_format: 
{
  fprintf(stderr,"tetgen2FV: parsing tetgen face file %s...\n",curr_file);
  polygon_count=0;
  plp=NULL;
  polygon_head=NULL;
  polygon_tail=NULL;
}
        face_header
	face_list
{

  tet_mesh->nboundaryfaces=polygon_count;
  tet_mesh->boundary_head=polygon_head;
  
  fprintf(stderr,"          ... %d boundary faces found\n",polygon_count);
}; 


newline_list: NEWLINE
        | newline_list NEWLINE
;

face_header: int_arg int_arg newline_list
;


face_list: face
	| face_list face
;


face: int_arg int_arg int_arg int_arg newline_list
{
  node_1=$<dbl>2;
  node_2=$<dbl>3;
  node_3=$<dbl>4;
  if ((pop=(struct polygon *)malloc(sizeof(struct polygon)))==NULL) {
    tetgenfaceerror("Cannot store polygon");
    return(1);
  }
  if ((plp=(struct polygon_list *)malloc(sizeof(struct polygon_list)))==NULL) {
    tetgenfaceerror("Cannot store polygon list");
    return(1);
  }
  pop->polygon_index=polygon_count;
  polygon_count++;
  pop->voronoi_index=0;
  pop->nnodes=3;
  pop->ntets=0;
  pop->node_index[0]=node_1-1;
  pop->node_index[1]=node_2-1;
  pop->node_index[2]=node_3-1;
  pop->shared_tet[0]=NULL;
  pop->shared_tet[1]=NULL;
  pop->srfnum=0;
  pop->bcnum=0;
  pop->domin=0;
  pop->domout=0;
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
	| int_arg int_arg int_arg int_arg int_arg newline_list
{
  node_1=$<dbl>2;
  node_2=$<dbl>3;
  node_3=$<dbl>4;
  srfnum=$<dbl>5;
  if ((pop=(struct polygon *)malloc(sizeof(struct polygon)))==NULL) {
    tetgenfaceerror("Cannot store polygon");
    return(1);
  }
  if ((plp=(struct polygon_list *)malloc(sizeof(struct polygon_list)))==NULL) {
    tetgenfaceerror("Cannot store polygon list");
    return(1);
  }
  pop->polygon_index=polygon_count;
  polygon_count++;
  pop->voronoi_index=0;
  pop->nnodes=3;
  pop->ntets=0;
  pop->node_index[0]=node_1-1;
  pop->node_index[1]=node_2-1;
  pop->node_index[2]=node_3-1;
  pop->shared_tet[0]=NULL;
  pop->shared_tet[1]=NULL;
  pop->srfnum=srfnum;
  pop->bcnum=0;
  pop->domin=0;
  pop->domout=0;
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

/*
	| NEWLINE
{
};
*/

int_arg: INTEGER {$$=(double)ival;}
;


real_arg: REAL {$$=rval;}
;


num_arg: INTEGER {$$=(double)ival;}
	| REAL {$$=rval;}
;


%%


#undef tetgenfacewrap

int tetgenfacewrap()
{
        return(1);
} 


int tetgenfaceerror(char *s)
{
	fprintf(stderr,"tetgen2FV: error on line: %d of file: %s  %s\n",
	        line_num,curr_file,s);
	fflush(stderr);
	return(1);
}

