//this file adds the header of the RAWIV data to the existing RAW file.
//converts a RAW to a RAWIV.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <conio.h>

short* vals;
int size;
char ifname[500], ofname[500];
int x, y, z;

size_t getShort(short *shts, size_t n, FILE *fp)
{
	unsigned char *pb = new unsigned char[n*2];
	unsigned char *ps = (unsigned char *)shts;

	size_t nbytes = fread(pb, sizeof(unsigned char), n*2, fp);
	//swap the byte order
	if(nbytes == n*2) {
		for(size_t i = 0; i < n; i++) {
			ps[2*i] = pb[2*i+1];
			ps[2*i+1] = pb[2*i];
		}
	}
	delete pb;
	return nbytes;
}

int putShort(short *sts, int n, FILE *fp)
{
#ifndef WIN32
	return fwrite(sts, sizeof(short), n, fp);
#else

	unsigned char *pb = new unsigned char[2*n]; 
	unsigned char *pf = (unsigned char *)sts;
	int	nbytes;

	//swap the byte order	
	for(int i = 0; i < n; i++) 
	{
		pb[2*i] = pf[2*i+1];
		pb[2*i+1] = pf[2*i];
	}

	nbytes=fwrite(pb, 1, 2*n, fp);
	delete pb;
	return nbytes;
#endif
}

int putInt(int *Ints, int n, FILE *fp)
{
#ifndef WIN32
	return fwrite(Ints, sizeof(int), n, fp);
#else

	unsigned char *pb = new unsigned char[4*n]; 
	unsigned char *pf = (unsigned char *)Ints;
	int	nbytes;

	//swap the byte order	
	for(int i = 0; i < n; i++) 
	{
		pb[4*i] = pf[4*i+3];
		pb[4*i+1] = pf[4*i+2];
		pb[4*i+2] = pf[4*i+1];
		pb[4*i+3] = pf[4*i];
	}

	nbytes=fwrite(pb, 1, 4*n, fp);
	delete pb;
	return nbytes;
#endif
}

int putFloat(float *flts, int n, FILE *fp)
{
#ifndef WIN32
	return fwrite(flts, sizeof(float), n, fp);
#else

	unsigned char *pb = new unsigned char[n*4];
	unsigned char *pf = (unsigned char *)flts;
    int nbytes;
	
	//swap the byte order	
	for(int i = 0; i < n; i++) 
	{ 
		pb[4*i] = pf[4*i+3];
		pb[4*i+1] = pf[4*i+2];
		pb[4*i+2] = pf[4*i+1];
		pb[4*i+3] = pf[4*i];
	}
	
	nbytes = fwrite(pb, 1, 4*n, fp);
	delete pb;
	return nbytes;
#endif
}

void read_file()
{
	FILE* fp;
	
	if ((fp = fopen(ifname, "rb")) == NULL) 
	{
		fprintf(stderr, "ERROR: fopen(%s)\n", ifname);
		exit(0);
    }
	
	size = x*y*z;
	vals = (short*)(malloc(sizeof(short) * size));
	
#ifndef WIN32
	fread(vals, sizeof(short), size, fp);
#else
	getShort(vals, size, fp);
#endif
				
}

void write_file()
{
    FILE* fp;
	char buff[4000];
	int temp;
	float tfloat;
	
	sprintf(buff, "%s.rawiv",ofname );
	
	if ((fp=fopen(buff,"wb")) == NULL)
	{
		printf("Cannot open the Output file for RAW output\n");
		exit(0);
	}
	
	printf("writing head info \n");
	
	//The origin 0,0,0
	tfloat = 0;
	putFloat(&tfloat,1,fp);  putFloat(&tfloat,1,fp);  putFloat(&tfloat,1,fp); 
	
	//The max (size+1),(size+1),(size+1)
	tfloat = x;		putFloat(&tfloat,1,fp);  
	tfloat = y;		putFloat(&tfloat,1,fp);  
	tfloat = z;		putFloat(&tfloat,1,fp); 
	
	//The #of vertices
	temp=size;
	putInt(&temp,1,fp);
	
	//The #of cells
	temp=(x-1)*(y-1)*(z-1);
	putInt(&temp,1,fp);
	
	//The dim of the volume
	temp = x;	putInt(&temp,1,fp); 
	temp = y;	putInt(&temp,1,fp);
	temp = z;	putInt(&temp,1,fp);
	
	//The Origin of the volume
	tfloat = 0;
	putFloat(&tfloat,1,fp);  putFloat(&tfloat,1,fp);  putFloat(&tfloat,1,fp); // origin
	
	//The span of the volume
	tfloat = 1;
	putFloat(&tfloat,1,fp); putFloat(&tfloat,1,fp); putFloat(&tfloat,1,fp); //span
	
	printf("writing data \n");
	
	putShort(&(vals[0]),size,fp);				
	
	fclose(fp);
	
}

int main(int argc, char *argv[])
{
	if (argc != 3) 
	{
		printf("Enter the file name 2 b read\n");
		exit(0);
	}

	strcpy(ifname, argv[1]);
	strcpy(ofname, argv[2]);

	printf("input filename is : %s\n", ifname);
	printf("output filename is : %s\n", ofname);

	x = 256;	 y = 256;	z = 110;
	read_file();
	write_file();

	printf("file written successfully\n");
	return 1;
}