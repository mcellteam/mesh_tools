#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "dx2mesh.h"

extern FILE *yyin;
char *infile_1;
char *infile_2;
char *curr_file;
int line_num=0;
int skip_freq;
struct vector3 translate;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{

	if (argc<3) {
      	  fprintf(stderr,"\nRead DX vertex file and connection file and convert to mesh format.\n ");
      	  fprintf(stderr,"  Output is written to stdout.\n\n ");
      	  fprintf(stderr,"  Usage: %s vertex_file_name connection_file_name\n\n",argv[0]);
	  exit(1);
	}

	infile_1=argv[1];
	infile_2=argv[2];

	if ((yyin=fopen(infile_1,"r"))==NULL) {
	  fprintf(stderr,"dx2mesh: error opening file: %s\n",infile_1);
	  exit(1);
	} 
        curr_file=infile_1;
	fflush(stdout);
	if (yyparse()) {
	  fprintf(stderr,"dx2mesh: error parsing mesh file %s\n",curr_file);
	  exit(1);
	} 
	fclose(yyin);

	exit(0);
}
