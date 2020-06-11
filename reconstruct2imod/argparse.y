%{
  #include <stdio.h> 
  #include <string.h> 
  #include <stdlib.h>
  #include <unistd.h>
  #include <stdarg.h>
  #include "strfunc.h"
  #include "reconstruct2imod.h"
  #include "argparse.h"

  extern struct name_list *file_name_list;
  extern int xmax, ymax, zmax;
  extern int start_slice_number, end_slice_number;
  extern double pixel_size;
  extern char *object_name;

  char *arg_cval;
  int arg_ival;
  double arg_rval;
  struct name_list *nlp;

  char *arg_err_msg;
%}


%union {
int tok;
char *str;
double dbl;
int intg;
} 


%{
  #include "arglex.flex.c"
%}


%name-prefix="arg"
%output="argparse.bison.c"

%token <tok> REAL INTEGER HELP_OPT OBJ_OPT DIM_OPT PIXEL_SIZE_OPT SLICES_OPT TEXT_ARG
%token <tok> EOF_TOK
%type <intg> int_arg
%type <dbl> real_arg num_arg 

%right '='
%left '+' '-'
%left '*' '/'
%left UNARYMINUS

%%

option_format: 
	option_list
;


option_list: option
	| option_list option
;


option: HELP_OPT
{
  return(1);
}
	| object_opt
	| dim_opt
	| pixel_size_opt
	| slices_opt
	| recon_infile_spec
	| EOF_TOK
{
  if (object_name == NULL) {
    sprintf(arg_err_msg,"No object name specified");
    argerror(arg_err_msg);
    return(1);
  }
  if (file_name_list == NULL) {
    sprintf(arg_err_msg,"No RECONSTRUCT base file name specified");
    argerror(arg_err_msg);
    return(1);
  }
  return(0);
};


object_opt: OBJ_OPT TEXT_ARG
{
  object_name=my_strdup(arg_cval);

  if (object_name == NULL) {
    sprintf(arg_err_msg,"Out of memory while parsing command line arguments: %s",arg_cval);
    argerror(arg_err_msg);
    return(1);
  }
  free((void *)arg_cval);
};


dim_opt: DIM_OPT int_arg int_arg int_arg
{
  xmax = $<intg>2;
  ymax = $<intg>3;
  zmax = $<intg>4;
};


pixel_size_opt: PIXEL_SIZE_OPT num_arg
{
  pixel_size = $<dbl>2;
};


slices_opt: SLICES_OPT int_arg int_arg
{
  start_slice_number = $<intg>2;
  end_slice_number = $<intg>3;
};


recon_infile_spec: TEXT_ARG
{
  if ((nlp=(struct name_list *)malloc(sizeof(struct name_list)))==NULL) {
    sprintf(arg_err_msg,"Out of memory storing base file name: %s",arg_cval);
    argerror(arg_err_msg);
    return(1);
  }
  nlp->next=file_name_list;
  file_name_list=nlp;

  nlp->name=my_strdup(arg_cval);

  if (nlp->name == NULL) {
    sprintf(arg_err_msg,"Out of memory while parsing command line arguments: %s",arg_cval);
    argerror(arg_err_msg);
    return(1);
  }
  free((void *)arg_cval);
};


int_arg: INTEGER {$$=arg_ival;}
;


real_arg: REAL {$$=arg_rval;}
;


num_arg: INTEGER {$$=(double)arg_ival;}
	| REAL {$$=arg_rval;}
;


%%


void argerror(char *s)
{
  FILE *log_file;

  log_file=stderr;

  fprintf(log_file,"\nrecon2obj: command-line argument syntax error: %s\n",s);
  fflush(log_file);
  return;
}


int argparse_init(int argc, char *argv[])
{
  FILE *err_file;
  struct yy_buffer_state *arg_input_buffer;
  char *tempstr,*arg_string;
  int i;


  err_file=stderr;
  file_name_list=NULL;
  object_name = NULL;
  xmax = 0;
  ymax = 0;
  zmax = 0;
  pixel_size = 0.05;
  start_slice_number = 0;
  end_slice_number = 0;

  if ((arg_err_msg=(char *)malloc(1024*sizeof(char)))==NULL) {
    fprintf(err_file,"recon2obj: Out of memory storing arg_err_msg\n");
    return(1);
  }

  arg_string=my_strdup("");
  for (i=1;i<argc;i++) {
    if (i==1) {
      free(arg_string);
      if ((arg_string=my_strdup(argv[i]))==NULL) {
        fprintf(err_file,"recon2obj: Out of memory storing arg_string\n");
        return(1);
      }
    }
    else {
      tempstr=arg_string; 
      if ((arg_string=my_strcat(arg_string," "))==NULL) {
        fprintf(err_file,"recon2obj: Out of memory storing arg_string\n");
        return(1);
      }
      free(tempstr);

      tempstr=arg_string; 
      if ((arg_string=my_strcat(arg_string,argv[i]))==NULL) {
        fprintf(err_file,"recon2obj: Out of memory storing arg_string\n");
        return(1);
      }
      free(tempstr);
    }
  }

  arg_input_buffer=arg_scan_string(arg_string);

  if (argparse()) {
     return(1);
  }

  arg_delete_buffer(arg_input_buffer);

  return(0);
}
