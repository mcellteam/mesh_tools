#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "strfunc.h"
#include "recon2obj.h"

extern FILE *reconin;
struct name_list *file_name_list;
struct object *objp;
struct section *section_head, *section_tail;
int line_num=1;
int start_slice_number, end_slice_number, curr_slice_number;
int vesicles_opt;
struct vector3 translate;
struct vertex_list **vertex_array;
char *object_name;
char *curr_file;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{
  FILE *log_file;
  struct name_list *nlp;
  char *base_file_name;
  char *argstr;
  char intstr[16];

  log_file=stderr;

  if (argparse_init(argc,argv)) {
    fprintf(log_file,"\n");
    fprintf(log_file,"Usage: %s [options] -object object_name recon_base_file_name start_slice_number end_slice_number\n\n",argv[0]);
    fprintf(log_file,"    options:\n");
    fprintf(log_file,"       [-help]                   print this help message\n");
    fprintf(log_file,"       [-contours | -vesicles]   treat contours as plain contours\n");
    fprintf(log_file,"                                 (this is the default) or convert the contours\n");
    fprintf(log_file,"                                 to points representing vesicles\n");
    fprintf(log_file,"\nRead stack of RECONSTRUCT files and convert traces to OBJ format.\n\n ");

    exit(1);
  }

  objp=NULL;
  section_head=NULL;
  section_tail=NULL;
  base_file_name = file_name_list->name;
  for (curr_slice_number=start_slice_number; curr_slice_number<end_slice_number+1; curr_slice_number++)
  {
    sprintf(intstr,"%d",curr_slice_number);
    curr_file = my_strcat(base_file_name, my_strcat(".",intstr));

    if ((reconin=fopen(curr_file,"r"))==NULL) {
      fprintf(log_file,"recon2obj: error opening RECONSTRUCT file: %s\n",curr_file);
      exit(1);
    } 
    fflush(stdout);
    if (reconparse()) {
      fprintf(log_file,"recon2obj: error parsing RECONSTRUCT file %s\n",curr_file);
      exit(1);
    } 
    fclose(reconin);
  }
  
  exit(0);
}

