#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "obj2mesh.h"

extern FILE *objin;
char *infile;
int line_num;
int skip_freq;

/* Begin main */
int main(int argc, char *argv[])
{

	if (argc<2) {
      	  fprintf(stderr,"Usage: %s in_file_name\n",argv[0]);
	  exit(1);
	}

	infile=argv[1];

	if ((objin=fopen(infile,"r"))==NULL) {
	  fprintf(stderr,"obj2mesh: error opening file: %s\n",infile);
	  exit(1);
	} 
	fflush(stdout);
	if (objparse()) {
	  fprintf(stderr,"obj2mesh: error parsing file: %s\n",infile);
	  exit(1);
	} 
	fclose(objin);

	exit(0);
}
