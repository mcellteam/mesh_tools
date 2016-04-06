#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "obj2mesh.h"

extern FILE *yyin;
char *infile;
int line_num;
int skip_freq;

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
	  fprintf(stderr,"obj2mesh: error opening file: %s\n",infile);
	  exit(1);
	} 
	fflush(stdout);
	if (yyparse()) {
	  fprintf(stderr,"obj2mesh: error parsing file: %s\n",infile);
	  exit(1);
	} 
	fclose(yyin);

	exit(0);
}
