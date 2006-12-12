#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "meshoffset.h"
#include "vector.h"

extern FILE *yyin;
char *infile;
int line_num;
int skip_freq;
double offset;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{

	if (argc<3) {
      	  fprintf(stderr,"Usage: %s offset in_file_name\n",argv[0]);
	  exit(1);
	}

        sscanf(argv[1],"%lf",&offset);
	infile=argv[2];

        fprintf(stderr,"offset set to: %g\n",offset);

	if ((yyin=fopen(infile,"r"))==NULL) {
	  fprintf(stderr,"meshoffset: error opening file: %s\n",infile);
	  exit(1);
	} 
	fflush(stdout);
	if (yyparse()) {
	  fprintf(stderr,"meshoffset: error parsing file: %s\n",infile);
	  exit(1);
	} 
	fclose(yyin);

	exit(0);
}
