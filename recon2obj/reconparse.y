%{
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "strfunc.h"
#include "reconstruct_transform.h"
#include "recon2obj.h"

#ifdef DEBUG
#define no_printf printf
#endif

extern int line_num;
extern FILE *reconin;
extern struct object *objp;
extern struct section *section_head, *section_tail;
extern double section_thickness;
extern int start_slice_number, end_slice_number, curr_slice_number;
extern int vesicles_opt;
extern char *object_name;
extern char *curr_file;
extern struct object *world_obj;
extern struct vertex_list **vertex_array;

int ival;
double rval;
char *strval;
char *tempstr;
char *cval;
char *var;
char *a_str,*rem_str;
char c,fmt_str[128],time_str[128];
char err_msg[128];
struct object_list *olp;
struct section *sp;
struct contour *cp, *contour_head;
struct double_list *dlp_head,*dlp;
struct vertex_list *vlp,*vertex_head,*vertex_tail;
struct vector3 *vecp;
int vertex_count;
double x,y,z,x_sum,y_sum;
double xcoef[6], ycoef[6];
int transform_dim;
int vert_1,vert_2,vert_3,vert_4;
int i;
int found_object;


%}

%union {
int tok;
char *str;
double dbl;
struct vector3 *vec;
struct object *obj;
} 


%{
  #include "reconlex.flex.c"
%}


%name-prefix="recon"
%output="reconparse.bison.c"


%token <tok> SECTION_BEGIN SECTION_END TAG_END
%token <tok> TRANSFORM_BEGIN TRANSFORM_END IMAGE_BEGIN
%token <tok> DIM XCOEF YCOEF
%token <tok> CONTOUR_BEGIN
%token <tok> ATTRIBUTE_NAME NAME POINTS HANDLES EOF_TOK
%token <tok> REAL INTEGER STR_VALUE
%type <dbl> int_arg real_arg num_arg 

%right '='
%left '+' '-'
%left '*' '/'
%left UNARYMINUS

%%

recon_format:
{

  transform_dim = 0;
  xcoef[0] = 0;
  xcoef[1] = 1;
  xcoef[2] = 0;
  xcoef[3] = 0;
  xcoef[4] = 0;
  xcoef[5] = 0;

  ycoef[0] = 0;
  ycoef[1] = 0;
  ycoef[2] = 1;
  ycoef[3] = 0;
  ycoef[4] = 0;
  ycoef[5] = 0;

}
  section EOF_TOK
{
  
  if (curr_slice_number == end_slice_number && objp != NULL)
  {
    if (objp->object_type==SECTIONS)
    {
      section_head = (struct section *)objp->contents;
      for (sp=section_head; sp!=NULL; sp=sp->next)
      {
        contour_head = sp->contour;
        for (cp=contour_head; cp!=NULL; cp=cp->next)
        {
          vertex_head = cp->vertex_list;
          x_sum = 0;
          y_sum = 0;
          for (vlp=vertex_head; vlp!=NULL; vlp=vlp->next)
          {
            x = vlp->vertex->x;
            y = vlp->vertex->y;
            z = vlp->vertex->z;
            if (vesicles_opt == 1)
            {
              x_sum += x;
              y_sum += y;
            }
            else
            {
              printf("v %.17g %.17g %.17g\n",x,y,z);
            }
          }

          vertex_count = cp->vertex_count;
          if (vesicles_opt == 1)
          {
            x = x_sum/vertex_count;
            y = y_sum/vertex_count;
            printf("v %.17g %.17g %.17g\n",x,y,z);
          }
          else
          {
            printf("l");
            for (i=0; i<vertex_count; i++)
            {
              printf(" %d",i-(vertex_count));
            }
            printf(" %d", -(vertex_count));
            printf("\n");
          }

//          fprintf(stderr,"\nContour: %s  vertices: %d\n",objp->name, vertex_count);
        }
      }
    }
  }

  return(0);
};


section: section_tag transform_list SECTION_END
;

section_tag: SECTION_BEGIN attribute_list TAG_END
{
  if (objp!=NULL)
  {
    /* Append a new section to end of section list*/
    if ((sp=(struct section *)malloc(sizeof(struct section)))==NULL)
    {
      reconerror("Out of memory while creating section");
      return(1);
    }
    section_tail->next=sp;
    section_tail=sp;
    section_tail->contour=NULL;
    section_tail->next=NULL;
  }
}
;

attribute_list: /* empty */
	| attribute_list attribute
;

attribute: ATTRIBUTE_NAME '=' STR_VALUE
;

transform_list: /* empty */
	| transform_list transform
;

transform: TRANSFORM_BEGIN transform_spec TAG_END
	image
	contour_list
  TRANSFORM_END
;

transform_spec: DIM int_arg '"'
                XCOEF num_arg num_arg num_arg num_arg num_arg num_arg '"'
                YCOEF num_arg num_arg num_arg num_arg num_arg num_arg '"'
{
  transform_dim = $<dbl>2;
  xcoef[0] = $<dbl>5;
  xcoef[1] = $<dbl>6;
  xcoef[2] = $<dbl>7;
  xcoef[3] = $<dbl>8;
  xcoef[4] = $<dbl>9;
  xcoef[5] = $<dbl>10;

  ycoef[0] = $<dbl>13;
  ycoef[1] = $<dbl>14;
  ycoef[2] = $<dbl>15;
  ycoef[3] = $<dbl>16;
  ycoef[4] = $<dbl>17;
  ycoef[5] = $<dbl>18;
}
;

image: /* empty */
	| IMAGE_BEGIN attribute_list TAG_END
;

contour_list: /* empty */
  | contour
	| contour_list contour
;

contour: CONTOUR_BEGIN name_spec attribute_list handles_spec points_spec TAG_END
;

name_spec: NAME '=' STR_VALUE
{
  found_object = 0;
  if ((tempstr=strip_quotes(strval))==NULL)
  {
    reconerror("Cannot store tempstr");
    return(1);
  }
  if (strcmp(tempstr,object_name)==0)
  {
//    fprintf(stderr,"Found object %s in file %s\n",object_name,curr_file);
    found_object = 1;
    if (objp==NULL)
    {
      if ((objp=(struct object *)malloc(sizeof(struct object)))==NULL)
      {
        reconerror("Out of memory while creating object");
        return(1);
      }
      objp->name=tempstr;
      objp->object_type=SECTIONS;
      objp->parent=NULL;
      /* Add first section of new object */
      if ((sp=(struct section *)malloc(sizeof(struct section)))==NULL)
      {
        reconerror("Out of memory while creating section");
        return(1);
      }
      section_head=sp;
      section_tail=sp;
      section_tail->contour=NULL;
      section_tail->next=NULL;
      objp->contents=section_head;
    }
    /* Prepend new contour to contour list of this section */
    if ((cp=(struct contour *)malloc(sizeof(struct contour)))==NULL)
    {
      reconerror("Out of memory while creating contour");
      return(1);
    }
    cp->next=section_tail->contour;
    section_tail->contour=cp;
    cp->vertex_list = NULL;
    cp->vertex_count = 0;
    vertex_head = NULL;
    vertex_tail = NULL;
    vertex_count = 0;
  }
}
;

handles_spec:
   /* empty */
   |  HANDLES discard_list '"'
;

discard_list: discard_item4
	| discard_list discard_item4
;

discard_item4: num_arg num_arg num_arg num_arg ','
;

points_spec:
  /* empty */
  | POINTS vertex2_list '"'
{
  if (found_object)
  {
    cp->vertex_count = vertex_count;
    cp->vertex_list = vertex_head;
  }
}
;

vertex2_list: vertex2
	| vertex2_list vertex2
;

vertex2: num_arg num_arg ','
{
  if (found_object)
  {
    if ((vecp=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
      reconerror("Cannot store normal vector");
      return(1);
    }

    x=$<dbl>1;
    y=$<dbl>2;

    XYinverse(transform_dim, xcoef, ycoef, &x, &y);
//    x = -xcoef[0] + x*xcoef[1] + y*xcoef[2] + x*y*xcoef[3] + x*x*xcoef[4] * y*y*xcoef[5];
//    y = -ycoef[0] + x*ycoef[1] + y*ycoef[2] + x*y*ycoef[3] + x*x*ycoef[4] * y*y*ycoef[5];

    vecp->x=x;
    vecp->y=y;
    vecp->z=section_thickness*(curr_slice_number+0.5);

    if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
      reconerror("Cannot store vertex list");
      return(1);
    }
    vertex_count++;
    vlp->vertex=vecp;
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
};


int_arg: INTEGER {$$=(double)ival;}
;


real_arg: REAL {$$=rval;}
;


num_arg: int_arg {$$=$<dbl>1;}
	| real_arg {$$=$<dbl>1;}
;


%%


#undef reconwrap

reconwrap()
{
        return(1);
} 


reconerror(char *s)
{
	fprintf(stderr,"recon2obj: error on line: %d of file: %s  %s\n",
	        line_num,curr_file,s);
	fflush(stderr);
	return(1);
}

