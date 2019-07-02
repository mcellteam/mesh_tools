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

R ([-]*[0-9]+\.)|([-]*[0-9]*\.[0-9]+)|([-]*[0-9]+[eE][+-]?[0-9]+)|([-]*[0-9]+\.[0-9]*[eE][+-]?[0-9]+)|([-]*[0-9]*\.[0-9]+[eE][+-]?[0-9]+)
I ([-]*[0-9]+)


%%

[ \t]+                  ;  /* eat white space */
^"#"[^\n]*$             ;  /* a comment */
"imod"                  {return(IMOD);}
"max"                   {return(MAX);}
"offsets"               {return(OFFSETS);}
"angles"                {return(ANGLES);}
"refcurscale"           {return(REFCURSCALE);}
"refcurtrans"           {return(REFCURTRANS);}
"refcurrot"             {return(REFCURROT);}
"refoldtrans"           {return(REFOLDTRANS);}
"scale"                 {return(SCALE);}
"mousemode"             {return(MOUSEMODE);}
"drawmode"              {return(DRAWMODE);}
"b&w_level"             {return(BW_LEVEL);}
"resolution"            {return(RESOLUTION);}
"threshold"             {return(THRESHOLD);}
"pixsize"               {return(PIXSIZE);}
"units"                 {return(UNITS);}
"flipped"               {return(FLIPPED);}
"slicerAngle"           {return(SLICERANGLE);}
"currentview"           {return(CURRENTVIEW);}
"view"                  {return(VIEW);}
"viewfovy"              {return(VIEWFOVY);}
"viewcnear"             {return(VIEWCNEAR);}
"viewcfar"              {return(VIEWCFAR);}
"viewflags"             {return(VIEWFLAGS);}
"viewtrans"             {return(VIEWTRANS);}
"viewrot"               {return(VIEWROT);}
"viewlight"             {return(VIEWLIGHT);}
"viewlabel"             {return(VIEWLABEL);}
"depthcue"              {return(DEPTHCUE);}
"globalclips"           {return(GLOBALCLIPS);}
"object"                {return(OBJECT);}
"name"                  {return(NAME);}
"color"                 {return(COLOR);}
"Fillcolor"             {return(FILLCOLOR);}
"linewidth"             {return(LINEWIDTH);}
"surfsize"              {return(SURFSIZE);}
"pointsize"             {return(POINTSIZE);}
"axis"                  {return(AXIS);}
"width2D"               {return(WIDTH2D);}
"symbol"                {return(SYMBOL);}
"symsize"               {return(SYMSIZE);}
"symflags"              {return(SYMFLAGS);}
"ambient"               {return(AMBIENT);}
"diffuse"               {return(DIFFUSE);}
"specular"              {return(SPECULAR);}
"shininess"             {return(SHININESS);}
"obquality"             {return(OBQUALITY);}
"valblack"              {return(VALBLACK);}
"valwhite"              {return(VALWHITE);}
"matflags2"             {return(MATFLAGS2);}
"objclips"              {return(OBJCLIPS);}
"open"                  {return(OPEN);}
"closed"                {return(CLOSED);}
"scattered"             {return(SCATTERED);}
"fill"                  {return(FILL);}
"drawmesh"              {return(DRAWMESH);}
"nolines"               {return(NOLINES);}
"nodraw"                {return(NODRAW);}
"bothsides"             {return(BOTHSIDES);}
"insideout"             {return(INSIDEOUT);}
"usefill"               {return(USEFILL);}
"pntusefill"            {return(PNTUSEFILL);}
"pntonsec"              {return(PNTONSEC);}
"antialias"             {return(ANTIALIAS);}
"hastimes"              {return(HASTIMES);}
"usevalue"              {return(USEVALUE);}
"valcolor"              {return(VALCOLOR);}
"contour"               {return(CONTOUR);}
"contflags"             {return(CONTFLAGS);}
"conttime"              {return(CONTTIME);}
"mesh"                  {return(MESH);}
"Meshflags"             {return(MESHFLAGS);}
"Meshsurf"              {return(MESHSURF);}
"Meshtime"              {return(MESHTIME);}
{R}                     {sscanf((const char *)yytext,"%lf",&rval); return(REAL);}
{I}                     {sscanf((const char *)yytext,"%d",&ival); return(INTEGER);}
\"[^\"]*\"              {cval=my_strdup(yytext); return(STR_VALUE);}
[A-Za-z][A-Za-z0-9\_]*  {cval=my_strdup(yytext); return(VAR);}
[\,]                    {return(yytext[0]);}   /* a comma */
.                       ;  /* ignore unknown characters */
\n                      {line_num++; return(NL);}  /* new lines */

