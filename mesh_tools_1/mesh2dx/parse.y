%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "mesh2dx.h"

#include "lex.c"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *yyin;
extern FILE *outfile;
extern char *outfile_name;
extern char *infile_name;
extern char *curr_file_name;
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
int vertex_index;
int vertex_count;
int max_vertex;
int skip_count;
int polygon_count;
int polygon_value_count;
double x,y,z;
float xval,yval,zval;
int vert_1,vert_2,vert_3,vert_4;
int i;
int word;
byte *word_p;
char my_byte_order[8];


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
  polygon_count=0;
  polygon_value_count=0;
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

  word_p=(unsigned char *)&word;
  word=0x04030201;

  if (word_p[0]==1) {
    sprintf(my_byte_order,"lsb");
  }
  else {
    sprintf(my_byte_order,"msb");
  }

  fprintf(outfile,"object \"1\" array type float shape 3 items %d %s binary\n",vertex_count,my_byte_order);
  fprintf(outfile,"  data 0\n");
  fprintf(outfile,"  attribute \"dep\" string \"positions\"\n\n");
  fprintf(outfile,"object \"2\" array type int shape 3 items %d %s binary\n",polygon_count,my_byte_order);
  fprintf(outfile,"  data %d\n",3*sizeof(float)*vertex_count);
  fprintf(outfile,"  attribute \"ref\" string \"positions\"\n");
  fprintf(outfile,"  attribute \"element type\" string \"triangles\"\n\n");
  fprintf(outfile,"object \"3\" field\n");
  fprintf(outfile,"  component \"positions\" \"1\"\n");
  fprintf(outfile,"  component \"connections\" \"2\"\n\n");
/*
  fprintf(outfile,"object \"null_positions\" array\n\n");
  fprintf(outfile,"object \"null_connections\" array\n\n");
  fprintf(outfile,"object \"null_object\" field\n");
  fprintf(outfile,"  component \"positions\" \"null_positions\"\n");
  fprintf(outfile,"  component \"connections\" \"null_connections\"\n\n");
  fprintf(outfile,"object \"4\" group\n");
  fprintf(outfile,"  member \"null_object (default)\" \"null_object\"\n");
  fprintf(outfile,"  member \"0\" \"3\"\n\n");
*/
  fprintf(outfile,"end\n");
  vlp=vertex_head;
  while (vlp!=NULL) {
    xval=vlp->vertex->x;
    yval=vlp->vertex->y;
    zval=vlp->vertex->z;
    fwrite(&xval,sizeof(float),1,outfile);
    fwrite(&yval,sizeof(float),1,outfile);
    fwrite(&zval,sizeof(float),1,outfile);
    vlp=vlp->next;
  }
  plp=polygon_head;
  while (plp!=NULL) {
    if (plp->polygon->n_verts==3) {
      fwrite(plp->polygon->vertex_index,sizeof(int),3,outfile);
    }
    else {
      printf("Warning: skipping non-triangular polygon\n\n");
    }
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

#undef yywrap

yywrap()
{
        return(1);
} 


yyerror(s)
char *s;
{
	fprintf(stderr,"mesh2dx: error on line: %d of file: %s  %s\n",
	        line_num,curr_file_name,s);
	fflush(stderr);
	return(1);
}

