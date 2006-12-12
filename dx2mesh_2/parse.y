%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "dx2mesh_2.h"

#include "lex.c"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *yyin;
extern char *infile;
extern int skip_freq;
extern int positions_item_count,positions_data_format,positions_data_offset;
extern int connections_item_count,connections_data_format,connections_data_offset;
extern struct object *world_obj;

int ival;
double rval;
char *cval,*strval;
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
double x,y,z;
int vert_1,vert_2,vert_3;
int dx_file_data_offset;
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

%token <tok> END_STMT INTEGER REAL OBJECT POSITIONS_ARRAY CONNECTIONS_ARRAY
%token <tok> ITEMS LSB_ORDER MSB_ORDER BINARY_DATA DATA_STMT 
%token <tok> POSITIONS_DEP_ATTRIB POSITIONS_REF_ATTRIB CONNECTIONS_DEP_ATTRIB
%token <tok> TRIANGLES_ATTRIB POSITIONS_COMPONENT CONNECTIONS_COMPONENT
%token <tok> FIELD VALUE STR_VALUE
%type <dbl> int_arg real_arg num_arg item_count data_offset
%type <tok> byte_order data_format_def

%right '='
%left '+' '-'
%left '*' '/'
%left UNARYMINUS

%%

dx_format: 
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
  position_object
  connection_object
  default_field_object
  END_STMT
{ 
  line_num++;
  char_num++;
  dx_file_data_offset=char_num;
  positions_data_offset+=dx_file_data_offset;
  connections_data_offset+=dx_file_data_offset;
/*
  printf("\"end\" found on line number = %d\n",line_num);
  printf("File offset = %d\n",dx_file_data_offset);
  printf("Positions data offset = %d\n",positions_data_offset);
  printf("Connections data offset = %d\n",connections_data_offset);
*/
  return(0);
};

position_object: OBJECT int_or_string positions_array_def item_count data_format_def data_offset POSITIONS_DEP_ATTRIB
{
  positions_item_count=$<dbl>4;
  positions_data_format=$<tok>5;
  positions_data_offset=$<dbl>6;
};

positions_array_def: POSITIONS_ARRAY
{
};

item_count: ITEMS int_arg {$$=$<dbl>2;};

data_format_def: byte_order BINARY_DATA
{
  $$=$<tok>1;
};

byte_order: MSB_ORDER
{
  $$=0;
}
	| LSB_ORDER
{
  $$=1;
};

data_offset: DATA_STMT int_arg {$$=$<dbl>2;};

connection_object: OBJECT int_or_string connections_array_def item_count data_format_def data_offset connections_attribs
{
  connections_item_count=$<dbl>4;
  connections_data_format=$<tok>5;
  connections_data_offset=$<dbl>6;
};

connections_array_def: CONNECTIONS_ARRAY
{
};

connections_attribs: CONNECTIONS_DEP_ATTRIB POSITIONS_REF_ATTRIB TRIANGLES_ATTRIB
{
}
	| POSITIONS_REF_ATTRIB TRIANGLES_ATTRIB
{
}
	| TRIANGLES_ATTRIB CONNECTIONS_DEP_ATTRIB POSITIONS_REF_ATTRIB
{
}
	| TRIANGLES_ATTRIB POSITIONS_REF_ATTRIB
{
};


default_field_object: OBJECT int_or_string FIELD positions_field_component
connections_field_component
{
};

positions_field_component: POSITIONS_COMPONENT object_ref
{
};

connections_field_component: CONNECTIONS_COMPONENT object_ref
{
};

object_ref: VALUE int_arg
{
};
	| STR_VALUE
{
};

int_arg: INTEGER {$$=(double)ival;}
;

real_arg: REAL {$$=rval;}
;

num_arg: INTEGER {$$=(double)ival;}
	| REAL {$$=rval;}
;

int_or_string: int_arg
{
}
	| STR_VALUE
{
};

%%


yyerror(s)
char *s;
{
	fprintf(stderr,"dx2mesh_2: error on line: %d of file: %s  %s\n",
	        line_num,infile,s);
	fflush(stderr);
	return(1);
}

