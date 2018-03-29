%{

#include "strfunc.h"

extern int line_num;
extern char *strval;
extern int ival;
extern double rval;
extern char err_msg[128];

%}

%option prefix="recon"
%option outfile="reconlex.flex.c"
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
^"<?"[^\n]*$            ;
^"<!"[^\n]*$            ;
"<Contour"              {return(CONTOUR_BEGIN);}
[\"]"/>"                {return(CONTOUR_END);}
"name"                  {return(NAME);}
"points="[\"]           {return(POINTS);}
"dim="[\"]              {return(DIM);}
"xcoef="[\"]            {return(XCOEF);}
"ycoef="[\"]            {return(YCOEF);}
"<Section"              {return(SECTION_BEGIN);}
"</Section>"            {return(SECTION_END);}
"<Image"                {return(IMAGE_BEGIN);}
"<Transform"            {return(TRANSFORM_BEGIN);}
"</Transform>"          {return(TRANSFORM_END);}
">" |
"/>"                    {return(TAG_END);}
{R}			{sscanf((const char *)recontext,"%lf",&rval); return(REAL);}
{I}			{sscanf((const char *)recontext,"%d",&ival); return(INTEGER);}
\"(([^\"^\n]*)|((\\\")*))*\" {
                               strval=my_strdup(recontext);
                               if(strval == NULL){
                                  fprintf(stderr,"recon2obj: Out of memory while parsing string: %s\n",recontext);
                                  exit(1);
                               }
                               return(STR_VALUE);
                             }
[A-Za-z][A-Za-z0-9\_\.]* {return(ATTRIBUTE_NAME);}
[\,\=\"]                   {return(recontext[0]);}
.			;
\n			{line_num++;}
<<EOF>>                 {return(EOF_TOK);}
