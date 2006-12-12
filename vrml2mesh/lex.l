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

[ \t]+			;
^"#"[^\n]*$		;
"Separator"		{return(SEPARATOR);}
"PerspectiveCamera"	{return(PERSPECTIVE_CAMERA);}
"orientation"		{return(ORIENTATION);}
"position"		{return(POSITION);}
"heightAngle"		{return(HEIGHT_ANGLE);}
"ShapeHints"		{return(SHAPE_HINTS);}
"vertexOrdering"	{return(VERTEX_ORDERING);}
"COUNTERCLOCKWISE"	{return(COUNTERCLOCKWISE);}
"shapeType"		{return(SHAPE_TYPE);}
"UNKNOWN_SHAPE_TYPE"	{return(UNKNOWN_SHAPE_TYPE);}
"faceType"		{return(FACE_TYPE);}
"CONVEX"		{return(CONVEX);}
"creaseAngle"		{return(CREASE_ANGLE);}
"NormalBinding"		{return(NORMAL_BINDING);}
"value"			{return(VALUE);}
"PER_VERTEX"		{return(PER_VERTEX);}
"Material"		{return(MATERIAL);}
"ambientColor"		{return(AMBIENT_COLOR);}
"diffuseColor"		{return(DIFFUSE_COLOR);}
"emissiveColor"		{return(EMISSIVE_COLOR);}
"specularColor"		{return(SPECULAR_COLOR);}
"shininess"		{return(SHININESS);}
"transparency"		{return(TRANSPARENCY);}
"Coordinate3"		{return(COORDINATE3);}
"point"			{return(POINT);}
"Normal"		{return(NORMAL);}
"vector"		{return(VECTOR);}
"TransformSeparator"	{return(TRANSFORM_SEPARATOR);}
"IndexedFaceSet"	{return(INDEXED_FACE_SET);}
"coordIndex"		{return(COORD_INDEX);}
"normalIndex"		{return(NORMAL_INDEX);}
{R}			{sscanf((const char *)yytext,"%lf",&rval); return(REAL);}
{I}			{sscanf((const char *)yytext,"%d",&ival); return(INTEGER);}
[\,\{\}\[\]]			{return(yytext[0]);}
.			;
\n			{line_num++;}

