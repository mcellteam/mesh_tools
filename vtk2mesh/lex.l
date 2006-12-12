%{

extern int line_num;
extern char *cval;
extern int ival;
extern double rval;
extern double n1,n2,n3;
extern char err_msg[128];
extern char buf_str[256];
char junk[256];
char *my_strdup();
double dummy;

%}

%option noyywrap

%p 10000
%a 10000 
%o 10000
%e 2000
%n 1000

R ([-]*[0-9]+\.)|([-]*[0-9]*\.[0-9]+)|([-]*[0-9]+[eE][+-]?[0-9]+)|([-]*[0-9]+\.[0-9]*[eE][+-]?[0-9]+)|([-]*[0-9]*\.[0-9]+[eE][+-]?[0-9]+)
I ([-]*[0-9]+)
N ({R}|{I})

%%

^"#".*		;
"POINTS"	{return(POINTS);}
"POLYGONS"	{return(POLYGONS);}
"float"		{return(FLOAT);}
{R}		{sscanf((const char *)yytext,"%lf",&rval); return(REAL);}
{I}		{sscanf((const char *)yytext,"%d",&ival); return(INTEGER);}
[ \t]+		;
.		;
\n		{line_num++;}
