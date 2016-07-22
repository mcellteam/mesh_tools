%{

#include "strfunc.h"

extern int line_num;
extern char *cval;
extern int ival;
extern double rval;
extern char err_msg[128];

%}

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
{R}			{sscanf((const char *)tetgenfacetext,"%lf",&rval); return(REAL);}
{I}			{sscanf((const char *)tetgenfacetext,"%d",&ival); return(INTEGER);}
.			;
\n			{line_num++;return(NEWLINE);}
