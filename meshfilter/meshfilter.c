#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "meshfilter.h"
#include "vector.h"

extern FILE *yyin;
char *infile;
int line_num;
int skip_freq;
double cos_threshold;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{

	if (argc<3) {
      	  fprintf(stderr,"Usage: %s cutoff_angle in_file_name\n",argv[0]);
	  exit(1);
	}

        sscanf(argv[1],"%lf",&cos_threshold);
	infile=argv[2];

        cos_threshold=cos(fabs(cos_threshold)*MY_PI/180.0);
        fprintf(stderr,"Cosine threshold set to: %g\n",cos_threshold);

	if ((yyin=fopen(infile,"r"))==NULL) {
	  fprintf(stderr,"meshfilter: error opening file: %s\n",infile);
	  exit(1);
	} 
	fflush(stdout);
	if (yyparse()) {
	  fprintf(stderr,"meshfilter: error parsing file: %s\n",infile);
	  exit(1);
	} 
	fclose(yyin);

	exit(0);
}
