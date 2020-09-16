%{
  #include <stdio.h>
  #include <string.h>
  #include <stdlib.h>
  #include <unistd.h>
  #include "strfunc.h"
  #include "argparse.h"

  extern char *arg_cval;
  extern int arg_ival;
  extern double arg_rval;
%}


%option prefix="arg"
%option outfile="arglex.flex.c"
%option noyywrap

%p 10000
%a 10000 
%o 10000
%e 2000
%n 1000

R ([0-9]+\.)|([0-9]*\.[0-9]+)|([0-9]+[eE][+-]?[0-9]+)|([0-9]+\.[0-9]*[eE][+-]?[0-9]+)|([0-9]*\.[0-9]+[eE][+-]?[0-9]+)
I ([0-9]+)


%%

[ \t]+                  ;
"-help"			{return(HELP_OPT);}
"-object"		{return(OBJ_OPT);}
"-dim"		        {return(DIM_OPT);}
"-pixel_size"           {return(PIXEL_SIZE_OPT);}
"-slices"               {return(SLICES_OPT);}
{I}			{sscanf((const char *)argtext,"%d",&arg_ival); return(INTEGER);}
{R}			{sscanf((const char *)argtext,"%lf",&arg_rval); return(REAL);}
[[:alnum:][:punct:]]+   {arg_cval=my_strdup(argtext); 
			 if(arg_cval == NULL){
                           fprintf(stderr,"MCell: Out of memory while parsing command line arguments.\n");
                           exit(1);
                         }
			 return(TEXT_ARG);}
.			;
<<EOF>>			{return(EOF_TOK);}
