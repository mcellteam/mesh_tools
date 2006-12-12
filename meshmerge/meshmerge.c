#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "meshmerge.h"

extern FILE *yyin;
char *infile_1;
char *infile_2;
char *curr_file;
int line_num=0;
int parse_seq=0;
int skip_freq;
struct vector3 translate;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{

	if (argc<6) {
      	  fprintf(stderr,"Usage: %s translate_x translate_y translate_z in_file_name_1 in_file_name_2\n",argv[0]);
	  exit(1);
	}

        sscanf(argv[1],"%lf",&translate.x);
        sscanf(argv[2],"%lf",&translate.y);
        sscanf(argv[3],"%lf",&translate.z);
	infile_1=argv[4];
	infile_2=argv[5];

	if ((yyin=fopen(infile_1,"r"))==NULL) {
	  fprintf(stderr,"meshmerge: error opening file: %s\n",infile_1);
	  exit(1);
	} 
        curr_file=infile_1;
	fflush(stdout);
	if (yyparse()) {
	  fprintf(stderr,"meshmerge: error parsing mesh file %s\n",curr_file);
	  exit(1);
	} 
	fclose(yyin);

	exit(0);
}
