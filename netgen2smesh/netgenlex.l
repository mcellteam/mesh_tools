%{

#include "strfunc.h"

extern int line_num;
extern char *cval;
extern int ival;
extern double rval;
extern char err_msg[128];

%}

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
"csgsurfaces"		{return(CSG_SURFACES);}
"cylinder"		{return(CYLINDER);}
"dimension"		{return(DIMENSION);}
"edgesegmentsgi2"	{return(EDGE_SEGMENTS);}
"endmesh"	        {return(END_MESH);}
"mesh3d"		{return(MESH3D);}
"plane"		        {return(PLANE);}
"points"		{return(POINTS);}
"sphere"		{return(SPHERE);}
"surfaceelementsgi"	{return(SURFACE_ELEMENTS);}
"volumeelements"	{return(VOLUME_ELEMENTS);}
{R}			{sscanf((const char *)netgentext,"%lf",&rval); return(REAL);}
{I}			{sscanf((const char *)netgentext,"%d",&ival); return(INTEGER);}
.			;
\n			{line_num++;return(NEWLINE);}
