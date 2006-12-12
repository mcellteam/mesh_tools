#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "meshscale.h"

extern FILE *yyin;
char *infile_1;
char *infile_2;
char *curr_file;
int line_num=0;
int parse_seq=0;
int skip_freq;
struct vector3 scale;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{

	if (argc<5) {
      	  fprintf(stderr,"Usage: %s scale_x scale_y scale_z in_file_name\n",argv[0]);
	  exit(1);
	}

        sscanf(argv[1],"%lf",&scale.x);
        sscanf(argv[2],"%lf",&scale.y);
        sscanf(argv[3],"%lf",&scale.z);
	infile_1=argv[4];

	if ((yyin=fopen(infile_1,"r"))==NULL) {
	  fprintf(stderr,"meshscale: error opening file: %s\n",infile_1);
	  exit(1);
	} 
        curr_file=infile_1;
	fflush(stdout);
	if (yyparse()) {
	  fprintf(stderr,"meshscale: error parsing mesh file %s\n",curr_file);
	  exit(1);
	} 
	fclose(yyin);

	exit(0);
}
