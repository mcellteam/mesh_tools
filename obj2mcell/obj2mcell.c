#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "obj2mcell.h"

extern FILE *yyin;
char *infile;
char *curr_file;
char *object_name;
int line_num;
int skip_freq;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{
  char *argstr;

        if (argc<2) {
          fprintf(stderr,"\nRead mesh input and convert to MCell MDL format.\n");
          fprintf(stderr,"  Read from stdin if mesh_file_name is absent.\n");
          fprintf(stderr,"  Output is written to stdout.\n\n");
      	  fprintf(stderr,"  Usage: %s [-h] object_name [mesh_file_name]\n\n",argv[0]);
          fflush(stdout);
          exit(1);
        }
        if (argc>=2) {
          argstr=argv[1];
          if (strcmp(argstr,"-h")==0) {
            fprintf(stderr,"\nRead mesh input and convert to MCell MDL output.\n");
            fprintf(stderr,"  Read from stdin if mesh_file_name is absent.\n");
            fprintf(stderr,"  Output is written to stdout.\n\n");
      	    fprintf(stderr,"  Usage: %s [-h] object_name [mesh_file_name]\n\n",argv[0]);
            fflush(stdout);
            exit(1);
          }
        }
        if (argc==3) {
	  object_name=argv[1];
	  infile=argv[2];
          if ((yyin=fopen(infile,"r"))==NULL) {
            fprintf(stderr,"mesh2mcell: error opening file: %s\n",infile);
            fflush(stdout);
            exit(1);
          }
          curr_file=infile;
        }
        else {
	  object_name=argv[1];
          yyin=stdin;
          curr_file="stdin";
        }

	if (yyparse()) {
	  fprintf(stderr,"mesh2mcell: error parsing file: %s\n",infile);
	  exit(1);
	} 
	fclose(yyin);

	exit(0);
}
