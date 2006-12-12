#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mesh2vtk.h"

extern FILE *yyin;
char *infile;
char *curr_file;
int line_num;
int skip_freq;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{
  char *argstr;

        
	if (argc>=2) {
	  argstr=argv[1];
          if (strcmp(argstr,"-h")==0) {
      	    fprintf(stderr,"\nRead mesh input and convert to VTK output.\n");
      	    fprintf(stderr,"  Read from stdin if mesh_file_name is absent.\n");
      	    fprintf(stderr,"  Output is written to stdout.\n\n");
      	    fprintf(stderr,"  Usage: %s [-h] [mesh_file_name]\n\n",argv[0]);
	    fflush(stdout);
	    exit(1);
          }
          else {
	    infile=argv[1];
	    if ((yyin=fopen(infile,"r"))==NULL) {
	      fprintf(stderr,"mesh2vtk: error opening file: %s\n",infile);
	      fflush(stdout);
	      exit(1);
	    } 
            curr_file=infile;
          }
	}
        else {
          yyin=stdin;
          curr_file="stdin";
        }

	if (yyparse()) {
	  fprintf(stderr,"mesh2vtk: error parsing file: %s\n",curr_file);
	  exit(1);
	} 
	fclose(yyin);

	exit(0);
}
