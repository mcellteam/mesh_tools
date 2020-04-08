%{

#include "strfunc.h"

extern int line_num;
extern char *cval;
extern int ival;
extern double rval;
extern char err_msg[128];

%}

%option prefix="netgen"
%option outfile="netgenlex.flex.c"
%option noyywrap

%p 10000
%a 10000 
%o 10000
%e 2000
%n 1000

R ([-]*[0-9]+\.)|([-]*[0-9]*\.[0-9]+)|([-]*[0-9]+[eE][+-]?[0-9]+)|([-]*[0-9]+\.[0-9]*[eE][+-]?[0-9]+)|([-]*[0-9]*\.[0-9]+[eE][+-]?[0-9]+)
I ([-]*[0-9]+)


%%

[ \t]+			;
^"#"[^\n]*$		;
"csgsurfaces"		{return(CSGSURFACES);}
"cylinder"		{return(CYLINDER);}
"dimension"		{return(DIMENSION);}
"edgesegmentsgi2"	{return(EDGE_SEGMENTS);}
"endmesh"	        {return(ENDMESH);}
"face_colours"	        {return(FACE_COLOURS);}
"geomtype"	        {return(GEOMTYPE);}
"mesh3d"		{return(MESH3D);}
"plane"			{return(PLANE);}
"points"		{return(POINTS);}
"surfaceelementsgi"	{return(SURFACE_ELEMENTSGI);}
"surfaceelements"	{return(SURFACE_ELEMENTS);}
"undef"			{return(UNDEF);}
"volumeelements"	{return(VOLUME_ELEMENTS);}
{R}			{sscanf((const char *)netgentext,"%lf",&rval); return(REAL);}
{I}			{sscanf((const char *)netgentext,"%d",&ival); return(INTEGER);}
.			;
\n			{line_num++;return(NEWLINE);}
