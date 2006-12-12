#include <stdio.h>
#include <string.h>
#include "vtk2rib.h"

extern FILE *yyin;
char *infile;
char *object_name;
char *curr_file;
int line_num=0;
int skip_freq;
struct vector3 translate;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{

	if (argc<2) {
      	  fprintf(stderr,"Usage: %s in_file_name\n",argv[0]);
	  exit(1);
	}

	infile=argv[1];

	if ((yyin=fopen(infile,"r"))==NULL) {
	  fprintf(stderr,"vtk2rib: error opening file: %s\n",infile);
	  exit(1);
	} 
        curr_file=infile;
	fflush(stdout);
	if (yyparse()) {
	  fprintf(stderr,"vtk2rib: error parsing mesh file %s\n",curr_file);
	  exit(1);
	} 
	fclose(yyin);

	exit(0);
}
