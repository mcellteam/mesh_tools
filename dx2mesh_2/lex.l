%{

extern int line_num,char_num;
extern char *cval,*strval;
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

^"#".*		{char_num+=yyleng;}
"object"	{char_num+=yyleng; return(OBJECT);}
"array type float shape 3" |
"class array type float rank 1 shape 3" {char_num+=yyleng; return(POSITIONS_ARRAY);}
"array type int shape 3" |
"class array type int rank 1 shape 3" {char_num+=yyleng; return(CONNECTIONS_ARRAY);}
"items"		{char_num+=yyleng; return(ITEMS);}
"lsb"		{char_num+=yyleng; return(LSB_ORDER);}
"msb"		{char_num+=yyleng; return(MSB_ORDER);}
"ieee" |
"binary"	{char_num+=yyleng; return(BINARY_DATA);}
"data"		{char_num+=yyleng; return(DATA_STMT);}
"attribute \"dep\" string \"positions\"" {char_num+=yyleng; return(POSITIONS_DEP_ATTRIB);}
"attribute \"ref\" string \"positions\"" {char_num+=yyleng; return(POSITIONS_REF_ATTRIB);}
"attribute \"dep\" string \"connections\"" {char_num+=yyleng; return(CONNECTIONS_DEP_ATTRIB);}
"attribute \"element type\" string \"triangles\"" {char_num+=yyleng; return(TRIANGLES_ATTRIB);}
"component \"positions\"" {char_num+=yyleng; return(POSITIONS_COMPONENT);}
"component \"connections\"" {char_num+=yyleng; return(CONNECTIONS_COMPONENT);}
"class field" |
"field"		{char_num+=yyleng; return(FIELD);}
"value"		{char_num+=yyleng; return(VALUE);}
"end"		{char_num+=yyleng; return(END_STMT);}
\"(([^\"]*)|((\\\")*))*\"  {char_num+=yyleng; strval=my_strdup(yytext); return(STR_VALUE);}
{R}		{char_num+=yyleng; sscanf((const char *)yytext,"%lf",&rval); return(REAL);}
{I}		{char_num+=yyleng; sscanf((const char *)yytext,"%d",&ival); return(INTEGER);}
[ \t]+		{char_num+=yyleng;}
.		{char_num+=yyleng;}
\n		{char_num+=yyleng; line_num++;}
