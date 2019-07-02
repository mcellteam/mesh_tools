%{
#include <stdlib.h> 
#include <stdio.h> 
#include <string.h> 
#include <limits.h>
#include <math.h>
#include "imod2reconstruct.h"
#include "parse.h"

#ifdef DEBUG
#define no_printf printf
#endif

extern FILE *yyin;
extern char *infile;
extern char *recon_ser_prefix;

int ival;
double rval;
char *cval;
char *var;
char *a_str,*rem_str;
char c,fmt_str[128],time_str[128];
char err_msg[128];

FILE *recon_ser_file;
char recon_ser_file_name[256];
char *obj_name;
struct section *section_array;
struct contour *contp;
struct vertex_list *vlp,*vertex_head,*vertex_tail,**vertex_array;
struct vertex_index_list *vilp,*vertex_index_head,*vertex_index_tail;
struct vector3 *vecp;
int n_verts;
int vertex_index;
int vertex_count;
unsigned int num_sections, min_section, max_section, contour_z;
double x,y,z;
double pixsize;
double section_thickness;
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


%token <tok> IMOD
%token <tok> MAX
%token <tok> OFFSETS
%token <tok> ANGLES
%token <tok> REFCURSCALE
%token <tok> REFCURTRANS
%token <tok> REFCURROT
%token <tok> REFOLDTRANS
%token <tok> SCALE
%token <tok> MOUSEMODE
%token <tok> DRAWMODE
%token <tok> BW_LEVEL
%token <tok> RESOLUTION
%token <tok> THRESHOLD
%token <tok> PIXSIZE
%token <tok> UNITS
%token <tok> FLIPPED
%token <tok> SLICERANGLE
%token <tok> CURRENTVIEW
%token <tok> VIEW
%token <tok> VIEWFOVY
%token <tok> VIEWCNEAR
%token <tok> VIEWCFAR
%token <tok> VIEWFLAGS
%token <tok> VIEWTRANS
%token <tok> VIEWROT
%token <tok> VIEWLIGHT
%token <tok> VIEWLABEL
%token <tok> DEPTHCUE
%token <tok> GLOBALCLIPS
%token <tok> OBJECT
%token <tok> NAME
%token <tok> COLOR
%token <tok> FILLCOLOR
%token <tok> LINEWIDTH
%token <tok> SURFSIZE
%token <tok> POINTSIZE
%token <tok> AXIS
%token <tok> WIDTH2D
%token <tok> SYMBOL
%token <tok> SYMSIZE
%token <tok> SYMFLAGS
%token <tok> AMBIENT
%token <tok> DIFFUSE
%token <tok> SPECULAR
%token <tok> SHININESS
%token <tok> OBQUALITY
%token <tok> VALBLACK
%token <tok> VALWHITE
%token <tok> MATFLAGS2
%token <tok> OBJCLIPS
%token <tok> OPEN
%token <tok> CLOSED
%token <tok> SCATTERED
%token <tok> FILL
%token <tok> DRAWMESH
%token <tok> NOLINES
%token <tok> NODRAW
%token <tok> BOTHSIDES
%token <tok> INSIDEOUT
%token <tok> USEFILL
%token <tok> PNTUSEFILL
%token <tok> PNTONSEC
%token <tok> ANTIALIAS
%token <tok> HASTIMES
%token <tok> USEVALUE
%token <tok> VALCOLOR
%token <tok> CONTOUR
%token <tok> CONTFLAGS
%token <tok> CONTTIME
%token <tok> MESH
%token <tok> MESHFLAGS
%token <tok> MESHSURF
%token <tok> MESHTIME
%token <tok> NL

%token <str> STR_VALUE
%token <str> VAR

%token <tok> REAL INTEGER
%type <dbl> int_arg num_arg 
/*
%type <dbl> real_arg
*/

%right '='
%left '+' '-'
%left '*' '/'
%left UNARYMINUS


%%


imod_format: 
{
  num_sections = 10000;
  min_section = UINT_MAX;
  max_section = 0;
  section_thickness = 1;
  pixsize = 1;
  if ((section_array=(struct section *)malloc(num_sections*sizeof(struct section)))==NULL)
  {
    yyerror("Cannot store section array");
    return(1);
  }
  for (i=0;i<num_sections;i++) 
  {
    section_array[i].z = i;
    section_array[i].contour_head = NULL;
    section_array[i].contour_tail = NULL;
    section_array[i].next = NULL;
  }
  obj_name = NULL;
  printf("\n");
}
  list_model_stmts
{
  sprintf(recon_ser_file_name,"%s.ser", recon_ser_prefix);
  if ((recon_ser_file=fopen(recon_ser_file_name,"w"))==NULL)
  {
    yyerror("Cannot open Reconstruct Series ser file");
    return(1);
  }
  section_thickness *= pixsize;
  fprintf(recon_ser_file,
           "<?xml version=\"1.0\"?>\n"
           "<!DOCTYPE Series SYSTEM \"series.dtd\">\n"
           "<Series index=\"152\" viewport=\"11.741 7.1974 0.00210015\"\n"
           "  units=\"microns\"\n"
           "  defaultThickness=\"%.6g\"\n"
           "  offset3D=\"0 0 0\"\n"
           "  type3Dobject=\"1\"\n"
           "  first3Dsection=\"0\"\n"
           "  last3Dsection=\"2147483647\"\n"
           "  >\n"
           "</Series>\n",
           section_thickness
         );
  fclose(recon_ser_file);

  printf("\nmin_section: %d   max_section: %d\n",min_section,max_section);
  printf("with capping sections: min_section: %d   max_section: %d\n\n",min_section-1,max_section+1);
  for (i=min_section-1;i<=max_section+1;i++)
  {
//    printf("SECTION: %d\n",i);
    sprintf(recon_ser_file_name, "%s.%d", recon_ser_prefix, i);
    if ((recon_ser_file=fopen(recon_ser_file_name,"w"))==NULL)
    {
      yyerror("Cannot open Reconstruct Series section file");
      return(1);
    }
    fprintf(recon_ser_file,
             "<?xml version=\"1.0\"?>\n"
             "<!DOCTYPE Section SYSTEM \"section.dtd\">\n"
             "\n"
             "<Section index=\"%d\" thickness=\"%.6g\" alignLocked=\"true\">\n"
             "<Transform dim=\"0\"\n"
             " xcoef=\" 0 1 0 0 0 0\"\n"
             " ycoef=\" 0 0 1 0 0 0\">\n",
             i,
             section_thickness
          );

    for (contp=section_array[i].contour_head;contp!=NULL;contp=contp->next)
    {
//      printf("  CONTOUR: %s\n",contp->name);
      fprintf(recon_ser_file,
               "<Contour name=\"%s\" hidden=\"false\" closed=\"true\" simplified=\"true\" border=\"1 0 0\" fill=\"1 0 0\" mode=\"11\"\n"
               "points=\"",
               contp->name
            );
      for (vlp=contp->vertex_list;vlp!=NULL;vlp=vlp->next)
      {
//        printf("    VERTEX: %.4g %.4g %.4g\n",vlp->vertex->x,vlp->vertex->y,vlp->vertex->z);
        fprintf(recon_ser_file,"%.6g %.6g,\n",pixsize*vlp->vertex->x,pixsize*vlp->vertex->y);
      }
      fprintf(recon_ser_file,"\"/>\n");
    }
    fprintf(recon_ser_file,
            "</Transform>\n"
            "\n"
            "</Section>\n"
          );
    fclose(recon_ser_file);
  }
}
;


list_model_stmts:
  model_stmt
  | list_model_stmts model_stmt
;


model_stmt:
  /* empty */ NL
  | IMOD int_arg NL
  | model_directive
  | view_directive
  | object_directive
  | contour_directive
  | mesh_directive
;


model_directive:
  OFFSETS num_arg num_arg num_arg NL
  | MAX num_arg num_arg num_arg NL
  | SCALE num_arg num_arg num_arg NL
    {
      section_thickness = $<dbl>4;
    }
  | MOUSEMODE int_arg NL
  | ANGLES num_arg num_arg num_arg NL
  | REFCURSCALE num_arg num_arg num_arg NL
  | REFCURTRANS num_arg num_arg num_arg NL
  | REFCURROT num_arg num_arg num_arg NL
  | REFOLDTRANS num_arg num_arg num_arg NL
  | DRAWMODE int_arg NL
  | BW_LEVEL int_arg ',' int_arg NL
  | RESOLUTION int_arg NL
  | THRESHOLD int_arg NL
  | PIXSIZE num_arg NL
    {
      pixsize = $<dbl>2;
    }
  | UNITS VAR NL
  | SLICERANGLE num_arg num_arg num_arg num_arg num_arg num_arg num_arg VAR NL
  | FLIPPED int_arg NL
  | CURRENTVIEW int_arg NL
;


view_directive:
  VIEW int_arg NL
  | VIEWFOVY num_arg NL
  | VIEWCNEAR num_arg NL
  | VIEWCFAR num_arg NL
  | VIEWFLAGS int_arg NL
  | VIEWTRANS num_arg num_arg num_arg NL
  | VIEWROT num_arg num_arg num_arg NL
  | VIEWLIGHT num_arg num_arg NL
  | DEPTHCUE num_arg num_arg NL
  | VIEWLABEL VAR NL
  | global_clip_stmt
; 


global_clip_stmt: GLOBALCLIPS int_arg int_arg num_arg int_arg NL list_clip_planes
;


list_clip_planes:
  clip_plane
  | list_clip_planes clip_plane
;


clip_plane: num_arg num_arg num_arg num_arg num_arg num_arg NL
;


object_directive:
  OBJECT int_arg int_arg int_arg NL
    {
      obj_name = NULL;
      printf("Object ID: %g  num_contours: %g \n", $2, $3);
    }
  | COLOR  num_arg num_arg num_arg num_arg NL
  | FILLCOLOR num_arg num_arg num_arg NL
  | NAME VAR NL
    {
      obj_name = cval;
      printf("  Object name: %s\n", cval);
    }
  | LINEWIDTH int_arg NL
  | SURFSIZE int_arg NL
  | POINTSIZE int_arg NL
  | AXIS int_arg NL
  | WIDTH2D int_arg NL
  | SYMBOL int_arg NL
  | SYMSIZE int_arg NL
  | SYMFLAGS int_arg NL
  | AMBIENT num_arg NL
  | DIFFUSE num_arg NL
  | SPECULAR num_arg NL
  | SHININESS num_arg NL
  | OBQUALITY num_arg NL
  | VALBLACK num_arg NL
  | VALWHITE num_arg NL
  | MATFLAGS2 int_arg NL
  | obj_clip_stmt
;

  
obj_clip_stmt: OBJCLIPS int_arg int_arg num_arg int_arg NL list_clip_planes
;


contour_directive:
 CONTOUR int_arg int_arg int_arg NL
   {
//     printf("CONTOUR %g %g %g\n", $<dbl>2, $<dbl>3, $<dbl>4);
     if ((contp=(struct contour *)malloc(sizeof(struct contour)))==NULL) {
       yyerror("Cannot store contour");
       return(1);
     }
     contp->name = obj_name;
     contp->vertex_count = 0;
     contp->vertex_list = NULL;
     contp->next = NULL;
     contour_z = 0;
     vertex_count=0;
     vertex_head=NULL;
     vertex_tail=NULL;
   }
   vertex_list
   {
     contp->vertex_count = vertex_count;
     contp->vertex_list = vertex_head;
     if (contour_z < min_section) { min_section = contour_z; }
     if (contour_z > max_section) { max_section = contour_z; }
     if (section_array[contour_z].contour_head == NULL) {
       section_array[contour_z].contour_head = contp;
     }
     if (section_array[contour_z].contour_tail == NULL) {
       section_array[contour_z].contour_tail = contp;
     }
     else {
       section_array[contour_z].contour_tail->next = contp;
     }
     section_array[contour_z].contour_tail = contp;
   }
 | CONTOUR int_arg int_arg int_arg int_arg NL
   {
//     printf("CONTOUR %g %g %g %g\n", $<dbl>2, $<dbl>3, $<dbl>4, $<dbl>5);
     if ((contp=(struct contour *)malloc(sizeof(struct contour)))==NULL) {
       yyerror("Cannot store contour");
       return(1);
     }
     vertex_count=0;
     vertex_head=NULL;
     vertex_tail=NULL;
   }
   vertex_list
   {
     contp->vertex_count = vertex_count;
     contp->vertex_list = vertex_head;
     if (section_array[contour_z].contour_head == NULL) {
       section_array[contour_z].contour_head = contp;
     }
     if (section_array[contour_z].contour_tail == NULL) {
       section_array[contour_z].contour_tail = contp;
     }
     else {
       section_array[contour_z].contour_tail->next = contp;
     }
     section_array[contour_z].contour_tail = contp;
   }
;


vertex_list: vertex
	| vertex_list vertex
;


vertex: num_arg num_arg num_arg NL
{
  if (obj_name == NULL)
  {
    yyerror("Cannot create contours for unamed object");
    return(1);
  }
  if ((vecp=(struct vector3 *)malloc(sizeof(struct vector3)))==NULL) {
    yyerror("Cannot store normal vector");
    return(1);
  }
  vecp->x=$<dbl>1;
  vecp->y=$<dbl>2;
  vecp->z=$<dbl>3+1;
  contour_z = $<dbl>3+1;
//  printf("VERTEX  %.4g %.4g %.4g\n",vecp->x, vecp->y, vecp->z);
  if ((vlp=(struct vertex_list *)malloc(sizeof(struct vertex_list)))==NULL) {
    yyerror("Cannot store vertex list");
    return(1);
  }
  vertex_count++;
  vertex_index=vertex_count;
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
}
;


mesh_directive:
  MESH int_arg int_arg int_arg NL mesh_vertex_normal_list mesh_index_list
;


mesh_vertex_normal_list:
  mesh_vertex_normal
  | mesh_vertex_normal_list mesh_vertex_normal
;


mesh_vertex_normal: num_arg num_arg num_arg num_arg num_arg num_arg NL
;


mesh_index_list: 
  mesh_index
  | mesh_index_list mesh_index
;


mesh_index: int_arg NL
;


int_arg: INTEGER {$$=(double)ival;}
;


/*
real_arg: REAL {$$=rval;}
;
*/


num_arg: INTEGER {$$=(double)ival;}
	| REAL {$$=rval;}
;


%%


int yyerror(s)
char *s;
{
	fprintf(stderr,"imod2reconstruct: error on line: %d of file: %s  %s\n",
	        line_num,infile,s);
	fflush(stderr);
	return(1);
}

