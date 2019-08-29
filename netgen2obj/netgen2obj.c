#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "netgen2obj.h"

extern FILE *netgenin;
char *infile_1;
char *infile_2;
char *curr_file;
int line_num=0;
int skip_freq;
struct vector3 translate;
struct vertex_list **vertex_array;

/* Begin main */
main(argc,argv)
  int argc; 
  char *argv[];  
{
  char *argstr;

	if (argc<2) {
      	  fprintf(stderr,"\nRead NetGen volume mesh file and convert to obj mesh format.\n ");
      	  fprintf(stderr,"  Output is written to stdout.\n\n ");
      	  fprintf(stderr,"  Usage: %s vol_file_name \n\n",argv[0]);
	  exit(1);
	}
        else {
          argstr=argv[1];
          if (strcmp(argstr,"-h")==0) {
      	    fprintf(stderr,"\nRead NetGen volume mesh file convert to obj mesh format.\n ");
      	    fprintf(stderr,"  Output is written to stdout.\n\n ");
      	    fprintf(stderr,"  Usage: %s vol_file_name \n\n",argv[0]);
	    exit(1);
          }
        }

	infile_1=argv[1];

	if ((netgenin=fopen(infile_1,"r"))==NULL) {
	  fprintf(stderr,"netgen2obj: error opening volume file: %s\n",infile_1);
	  exit(1);
	} 
        curr_file=infile_1;
	fflush(stdout);
	if (netgenparse()) {
	  fprintf(stderr,"netgen2obj: error parsing volume file %s\n",curr_file);
	  exit(1);
	} 
	fclose(netgenin);

	exit(0);
}
