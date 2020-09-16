#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "strfunc.h"
#include "reconstruct2imod.h"

extern FILE *reconin;
struct name_list *file_name_list;
struct object *objp;
struct section *section_head, *section_tail;
double pixel_size;
int line_num;
int xmax, ymax, zmax, start_slice_number, end_slice_number, curr_slice_number;
int contour_count;
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
    fprintf(log_file,"Usage: %s [options] -object object_name -dim xmax ymax zmax -pixel_size z_value -slices start_slice_number end_slice_number  recon_series_base_file_name\n\n",argv[0]);
    fprintf(log_file,"    options:\n");
    fprintf(log_file,"       [-help]                   print this help message\n");
    fprintf(log_file,"\nRead stack of RECONSTRUCT files and convert traces of named object to IMOD ascii format.\n");
    fprintf(log_file,"  Output is written to stdout\n\n");
    fprintf(log_file,"\n  Example:\n\n ");
    fprintf(log_file,"    reconstruct2imod -object OuterMembrane -dim 1464 1493 1001 -pixel_size 0.00164 -slices 579 699 mito9 > mito9_outermembrane.amod\n\n ");

    exit(1);
  }

  objp=NULL;
  section_head=NULL;
  section_tail=NULL;
  base_file_name = file_name_list->name;
  contour_count = 0;
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

