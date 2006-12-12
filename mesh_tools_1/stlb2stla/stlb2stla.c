#include "stdio.h"
#include "stdlib.h"
#include "math.h"

/*
	Convert a binary stl file to ascii stl format
*/
#define TRUE  1
#define FALSE 0

#define SWAP FALSE

typedef struct {
	float x,y,z;
} XYZ;

int ReadFloat(FILE *,float *,int);
int ReadInt(FILE *,int *,int);

int main(int argc,char **argv)
{
	FILE *fptr;
	int i,j;
	int nf = 0;
	float x,y,z;
	XYZ n,p[3];
	char s[256];
        char *argstr;

	if (argc < 2) {
		fprintf(stderr,"\nRead binary stl file and convert to ascii stl format\n");
		fprintf(stderr,"  Output is written to stdout\n\n");
		fprintf(stderr,"  Usage: %s [-h] stl_file_name\n\n",argv[0]);
		exit(-1);
	}
        if (argc>=2)
        {
          argstr=argv[1];
          if (strcmp(argstr,"-h")==0)
          {
		fprintf(stderr,"\nRead binary stl file and convert to ascii stl format\n");
		fprintf(stderr,"  Output is written to stdout\n\n");
		fprintf(stderr,"  Usage: %s [-h] stl_file_name\n\n",argv[0]);
		exit(-1);
          }
        }

	if ((fptr = fopen(argv[1],"r")) == NULL) {
		fprintf(stderr,"Unable to open binary stl file\n");
		exit(-1);
	}

	/* Read the header */
	if (fread(s,sizeof(char),80,fptr) != 80) {
		fprintf(stderr,"Error reading header\n");
		exit(-1);
	}
	s[80] = '\0';
	fprintf(stderr,"Header: %s\n",s);
	ReadInt(fptr,&nf,SWAP);
	fprintf(stderr,"Expecting %d facets\n",nf);

        printf("solid\n");

	for (i=0;i<nf;i++) {

		
		/* Read the normal */
		ReadFloat(fptr,&n.x,SWAP);
      		ReadFloat(fptr,&n.y,SWAP);
      		ReadFloat(fptr,&n.z,SWAP);

          printf("  facet normal %.17g %.17g %.17g\n",n.x,n.y,n.z);
          printf("  outer loop\n");

		/* Read the vertices */
		for (j=0;j<3;j++) {
      	ReadFloat(fptr,&p[j].x,SWAP);
      	ReadFloat(fptr,&p[j].y,SWAP);
      	ReadFloat(fptr,&p[j].z,SWAP);
		}

		/* Read the padding */
		fgetc(fptr);
		fgetc(fptr);

		for (j=0;j<3;j++) {
		  printf("    vertex %.17g %.17g %.17g\n",p[j].x,p[j].y,p[j].z);
                }

          printf("  endloop\n");
          printf("  endfacet\n");

	}

        printf("endsolid\n");

	fprintf(stderr,"Wrote %d facets\n",nf);
	fclose(fptr);
}


/*
   Read a possibly byte swapped integer
*/
int ReadInt(FILE *fptr,int *n,int swap)
{
   unsigned char *cptr,tmp;

   if (fread(n,4,1,fptr) != 1)
      return(FALSE);
   if (swap) {
      cptr = (unsigned char *)n;
      tmp = cptr[0];
      cptr[0] = cptr[3];
      cptr[3] = tmp;
      tmp = cptr[1];
      cptr[1] = cptr[2];
      cptr[2] = tmp;
   }
   return(TRUE);
}


int ReadFloat(FILE *fptr,float *n,int swap)
{
   unsigned char *cptr,tmp;

   if (fread(n,4,1,fptr) != 1)
      return(FALSE);
   if (swap) {
      cptr = (unsigned char *)n;
      tmp = cptr[0];
      cptr[0] = cptr[3];
      cptr[3] =tmp;
      tmp = cptr[1];
      cptr[1] = cptr[2];
      cptr[2] = tmp;
   }
   return(TRUE);
}

