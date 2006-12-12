#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mesh2dx.h"

extern FILE *yyin;
FILE *outfile;
char *infile_name;
char *outfile_name;
char *curr_file_name;
int line_num;
int skip_freq;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{
  char *argstr;

        
        if (argc < 3) {
      	  fprintf(stderr,"\nRead mesh input and convert to DX output.\n");
      	  fprintf(stderr,"  Read from stdin if mesh_input_file is -.\n");
      	  fprintf(stderr,"  Output is written to dx_output_file\n\n");
      	  fprintf(stderr,"  Usage: %s [-h] mesh_input_file dx_output_file\n\n",argv[0]);
	  fflush(stdout);
	  exit(1);
        }
	argstr=argv[1];
        if (strcmp(argstr,"-h")==0) {
      	  fprintf(stderr,"\nRead mesh input and convert to DX output.\n");
      	  fprintf(stderr,"  Read from stdin if mesh_input_file is -.\n");
      	  fprintf(stderr,"  Output is written to dx_output_file\n\n");
      	  fprintf(stderr,"  Usage: %s [-h] mesh_input_file dx_output_file\n\n",argv[0]);
	  fflush(stdout);
	  exit(1);
        }
	if (argc>=3) {
	  infile_name=argv[1];
          if (strcmp(infile_name,"-")==0) {
            curr_file_name="stdin";
            yyin=stdin;
          }
	  else {
            if ((yyin=fopen(infile_name,"r"))==NULL) {
	      fprintf(stderr,"mesh2dx: error opening file: %s\n",infile_name);
	      fflush(stdout);
	      exit(1);
            }
            curr_file_name=infile_name;
	  } 
	  outfile_name=argv[2];
	  if ((outfile=fopen(outfile_name,"w"))==NULL) {
	    fprintf(stderr,"mesh2dx: error opening file: %s\n",outfile_name);
	    fflush(stdout);
	    exit(1);
	  } 
	}

	if (yyparse()) {
	  fprintf(stderr,"mesh2dx: error parsing file: %s\n",curr_file_name);
	  exit(1);
	} 
	fclose(yyin);
	fclose(outfile);

	exit(0);
}
