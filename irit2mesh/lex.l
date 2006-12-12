%{

extern int line_num;
extern char *cval;
extern int ival;
extern double rval;
extern char err_msg[128];
char *my_strdup();

%}

%option noyywrap

%p 10000
%a 10000 
%o 10000
%e 2000
%n 1000

%%

[ \t]+			;
"BSPLINE"		{return(BSPLINE);}
"color"			{return(COLOR);}
"E3"			{return(E3);}
"KV"			{return(KV);}
"INTERNAL"		{return(INTERNAL);}
"NORMAL"		{return(NORMAL);}
"OBJECT"		{return(OBJECT);}
"P3"			{return(P3);}
"PLANE"			{return(PLANE);}
"POLYGON"		{return(POLYGON);}
"SURFACE"		{return(SURFACE);}
\"[^\"]*\"              {cval=my_strdup(yytext); return(A_STRING);}
([-]*[0-9]*\.[0-9]+) |
([-]*[0-9]+[eE][+-]?[0-9]+) |
([-]*[0-9]*\.[0-9]+[eE][+-]?[0-9]+) {sscanf((const char *)yytext,"%lf",&rval); return(REAL);}
([-]*[0-9]+)            {sscanf((const char *)yytext,"%d",&ival); return(INTEGER);}
[A-Za-z][A-Za-z0-9\_]*  {cval=my_strdup(yytext); return(VAR);}
[\,\(\)\/\-\+\=\[\]\{\}\|\>\*\#\:]      {return(yytext[0]);}
.			;
\n			{line_num++;}
