/*
%%%%%%%%%%%%%%%%%%%%%%%%%
The SHASTRA software is not in the Public Domain. It is distributed 
on a person to person basis, solely for educational use and 
permission is NOT granted for its transfer to anyone or 
for its use in any commercial product. 
There is NO warranty on the available software and neither Purdue 
University nor the Applied Algebra and Geometry group directed by C. Bajaj
accept responsibility for the consequences of its use.
%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*
*  slice.c  --  ct slice interface
*
*  programmer: Dan Schikore
*        date: 6/25/93
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <memory.h>
//#include <utils/utils.h> //KLC

#include "ct/ct.h"
#include "ct/sunras.h"
#include "ct/slc.h"
#include "ct/slice.h"

#define DEBUGx

#define MAX(x, y) ((x)>(y)?(x):(y))
#define MIN(x, y) ((x)<(y)?(x):(y))

static int SunRasError(char *);
static void flipl(unsigned char *);
static int rle_read(unsigned char *, int, int, FILE *, int);

extern int STATIC_CHQ;

void my_clear_slice()
{
	rle_read (NULL, STATIC_CHQ, STATIC_CHQ, NULL, STATIC_CHQ);	
}


/*------------------------------------------------------------------------
* CTSliceCompMinMaxD - compute min and max densities for a slice
*------------------------------------------------------------------------
*/
void CTSliceCompMinMaxD(CTSlice slice)
{
	double val;
	int i, j;
	unsigned long nnonzero;
	
	CTSliceMaxD(slice) = -1e100;
	CTSliceMinD(slice) = 1e100;
	nnonzero=0;
	
	if (CTSliceShorts(slice)) 
	{
		for (i=CTSliceMinX(slice); i<=CTSliceMaxX(slice); i++)
		{
			for (j=CTSliceMinY(slice); j<=CTSliceMaxY(slice); j++) 
			{
				if ((val=CTSliceShortData(slice, i, j)) > CTVolNoise(slice->vol))
					continue;
				if (val > CTSliceMaxD(slice))
					CTSliceMaxD(slice) = val;
				if (val < CTSliceMinD(slice))
					CTSliceMinD(slice) = val;
				if (val != 0)
					nnonzero++;
			}
		}
	}
	else if (CTSliceChars(slice)) 
	{
		for (i=CTSliceMinX(slice); i<=CTSliceMaxX(slice); i++)
		{
			for (j=CTSliceMinY(slice); j<=CTSliceMaxY(slice); j++) 
			{
				if ((val=CTSliceCharData(slice, i, j)) > CTVolNoise(slice->vol))
					continue;
				if (val > CTSliceMaxD(slice))
					CTSliceMaxD(slice) = val;
				if (val < CTSliceMinD(slice))
					CTSliceMinD(slice) = val;
				if (val != 0)
					nnonzero++;
			}
		}
	}
	else if (CTSliceFloats(slice)) 
	{
		for (i=CTSliceMinX(slice); i<=CTSliceMaxX(slice); i++)
		{
			for (j=CTSliceMinY(slice); j<=CTSliceMaxY(slice); j++) 
			{
				if ((val=CTSliceFloatData(slice, i, j)) > CTVolNoise(slice->vol))
					continue;
				if (val > CTSliceMaxD(slice))
					CTSliceMaxD(slice) = val;
				if (val < CTSliceMinD(slice))
					CTSliceMinD(slice) = val;
				if (val != 0)
					nnonzero++;
			}
		}
	}
	CTSliceRatio(slice) = nnonzero/
		(double)(CTSliceWidth(slice)*CTSliceHeight(slice));
}

/*---------------------------------------------------------------------
* fskip - skip a specified number of bytes in a file
*---------------------------------------------------------------------
*/
static void fskip(FILE *fp, unsigned long nbytes)
{
	char buf[2048];
	
	while (nbytes > 2048) 
	{
		fread(buf, 1, 2048, fp);
		nbytes-=2048;
	}
	
	if (nbytes > 0)
		fread(buf, 1, nbytes, fp);
}

/*-------------------------------------------------------------------
* CTReadCTF - read a slice from a ctf file
*-------------------------------------------------------------------
*/
static int CTReadCTF(FILE *fp, CTSlice slice, int x1, int y1, int x2, int y2)
{
	unsigned short delta[512], row[512];
	int i, del, left, right;
	long off;
	
	/*
	CTSliceXSize(slice) = 512;
	CTSliceYSize(slice) = 512;
	*/
	
	if (x1==-1) 
	{
		slice->x1=slice->y1=x1=y1=0;
		slice->x2=slice->y2=x2=y2=511;
		slice->data = (unsigned short *)malloc(sizeof(unsigned short)*(x2-x1+1)*
			(y2-y1+1));
	
		memset(slice->data, '\0', sizeof(short)*(x2-x1+1)*(y2-y1+1));
		CTSliceWidth(slice) = CTSliceHeight(slice) = 512;
	}
	
	/* skip past the 2048 byte header */
	fskip(fp, 2048);
	
	/* read the number of data values for each line */
	fread(delta, 2, 512, fp);
	
	/* compute offset for rows up to y1 */
	for (off=0, i=0; i<x1; i++)
		off+=delta[i]*2*2;
	
#ifdef DEBUG
	printf("skipping %ld values\n", off);
#endif
	
	/* skip to the x1'th row of data */
	if (off != 0)
		fskip(fp, off);

	for (; x1<=x2; x1++) 
	{
		/* skip to y1'th column of data in this row */
		del=delta[x1];
		
		/* compute how much data is present */
		left = (512)/2 - del;
		right = (512)/2 + del- 1;
		
		memset(row, '\0', 512*sizeof(short));
		fread(&row[left], 2, 2*del, fp);
		memcpy(&CTSliceShortData(slice, x1, y1), &row[y1],
            sizeof(short)*(y2-y1+1));
	}
	
	CTSliceCompMinMaxD(slice);
	return(0);
}

/*-------------------------------------------------------------------
* CTReadRAS - read a slice from a sun raster file
*-------------------------------------------------------------------
*/
static int CTReadRAS(FILE *fp, CTSlice slice, int x1, int y1, int x2, int y2)
{
	int               linesize,lsize,csize,isize,flipit,i,j,w,h,d,rv;
	unsigned char     *line, r[256], g[256], b[256];
	struct rasterfile sunheader;
	
	rv = 0;
	
	/* read in the Sun Rasterfile picture */
	flipit = 0;
	fread(&sunheader, sizeof(struct rasterfile), 1, fp);
	
	if (sunheader.ras_magic != RAS_MAGIC) 
	{
		flipl( (unsigned char *) &sunheader.ras_magic);
		if (sunheader.ras_magic == RAS_MAGIC)
			flipit = 1;
		else
			flipl( (unsigned char *) &sunheader.ras_magic);
	}
	
	if (sunheader.ras_magic != RAS_MAGIC)
		return( SunRasError("not a Sun rasterfile") );
	
	if (flipit) {
		flipl((unsigned char *) &sunheader.ras_width);
		flipl((unsigned char *) &sunheader.ras_height);
		flipl((unsigned char *) &sunheader.ras_depth);
		flipl((unsigned char *) &sunheader.ras_length);
		flipl((unsigned char *) &sunheader.ras_type);
		flipl((unsigned char *) &sunheader.ras_maptype);
		flipl((unsigned char *) &sunheader.ras_maplength);
	}
	
	/* make sure that the input picture can be dealt with */
	if (sunheader.ras_depth != 8) 
	{
		fprintf (stderr, "Sun rasterfile image has depth %d\n",
			sunheader.ras_depth);
		fprintf (stderr, "Depths supported are 8\n");
		return 1;
	}
	
	if (sunheader.ras_type != RT_OLD && sunheader.ras_type != RT_STANDARD &&
		sunheader.ras_type != RT_BYTE_ENCODED &&
		sunheader.ras_type != RT_FORMAT_RGB) 
	{
		fprintf (stderr, "Sun rasterfile of unsupported type %d\n",
			sunheader.ras_type);
		return 1;
	}
	
	if (sunheader.ras_maptype != RMT_RAW && sunheader.ras_maptype != RMT_NONE &&
		sunheader.ras_maptype != RMT_EQUAL_RGB) 
	{
		fprintf (stderr, "Sun rasterfile colormap of unsupported type %d\n",
			sunheader.ras_maptype);
		return 1;
	}
	
	w = sunheader.ras_width;
	h = sunheader.ras_height;
	d = sunheader.ras_depth;
	isize = sunheader.ras_length ?
		sunheader.ras_length :
	(w * h * d) / 8;
	csize = (sunheader.ras_maptype == RMT_NONE) ? 0 : sunheader.ras_maplength;
	lsize = w * h * ( d == 1 ? d : d/8);
	linesize = w * d;
	/* if ((linesize % 48) && d == 24) linesize += (48 - (linesize % 48)); */
	if (linesize % 16) linesize += (16 - (linesize % 16));
	linesize /= 8;
	
#ifdef DEBUG
	fprintf(stderr,"%s: LoadSunRas() - loading a %dx%d pic, %d planes\n",
		cmd, w, h, d);
	fprintf(stderr,
		"type %d, maptype %d, isize %d, csize %d, lsize %d, linesize %d\n",
		sunheader.ras_type, sunheader.ras_maptype,
		isize, csize, lsize, linesize);
#endif
	
	
	/* read in the colormap, if any */
	if (sunheader.ras_maptype == RMT_EQUAL_RGB && csize) 
	{
		fread (r, sizeof(unsigned char), sunheader.ras_maplength/3, fp);
		fread (g, sizeof(unsigned char), sunheader.ras_maplength/3, fp);
		fread (b, sizeof(unsigned char), sunheader.ras_maplength/3, fp);
	}
	else if (sunheader.ras_maptype == RMT_RAW && csize) 
	{
		/* we don't know how to handle raw colormap, ignore */
		fskip(fp, csize);
	}
	else
	{
		/* no colormap, make one up */
		if (sunheader.ras_depth == 1) 
		{
			r[0] = g[0] = b[0] = 0;
			r[1] = g[1] = b[1] = 255;
		}
		else if (sunheader.ras_depth == 8) 
		{
			for (i = 0; i < 256; i++)
				r[i] = g[i]  = b[i] = i;
		}
	}
	
	/*
	CTSliceXSize(slice) = w;
	CTSliceYSize(slice) = h;
	*/
	
	if (x1==-1) 
	{
		slice->x1=slice->y1=x1=y1=0;
		slice->x2=x2=w-1;
		slice->y2=y2=h-1;
		slice->cdata = (unsigned char *)malloc(sizeof(unsigned char)*(x2-x1+1)*
			(y2-y1+1));
		memset(slice->cdata, '\0', sizeof(char)*(x2-x1+1)*(y2-y1+1));
		CTSliceWidth(slice) = w;
		CTSliceHeight(slice) = h;
	}
	
	/* allocate memory for picture and read it in */
	/* note we may slightly overallocate here (if image is padded) */
	line = (unsigned char *) malloc (linesize);
	if (line == NULL) 
	{
		fprintf(stderr, "Can't allocate memory for image\n");
		exit(1);
	}
	
	for (i = 0; i < h; i++) 
	{
		if (sunheader.ras_type == RT_BYTE_ENCODED) 
		{
			if (rle_read (line, 1, linesize, fp, (i==0)) != linesize) break;
		}
		else 
		{
			if (fread (line, 1, linesize, fp) != (unsigned int)linesize)
				return (SunRasError ("file read error"));
		}
		
		if (slice->y1 <= i && i <= slice->y2) 
		{
			for (j=slice->x1; j<=slice->x2; j++)
				CTSliceCharData(slice, i, j) = line[j];
		}
		/*
		switch (d) {
		case 1:  SunRas1to8 (image + w * i, line, w);       break;
		case 8:  memcpy (image + w * i, line, w);           break;
		case 24: memcpy (image + w * i * 3, line, w * 3); break;
		}
		*/
	}
	
#ifdef DEBUG
	fprintf(stderr,"Sun ras: image loaded!\n");
#endif
	
	/*
	if (d == 24) {
	if (sunheader.ras_type != RT_FORMAT_RGB) fixBGR(image,w,h);
	rv = Conv24to8 (image, w, h, nc);
	free (image);
	return (rv);
	}
	else {
	pic = image;
	pWIDE = w;
	pHIGH = h;
	return (0);
	}
	*/
	CTSliceCompMinMaxD(slice);
	return(0);
}

/*****************************/
static int rle_read (unsigned char* ptr, int size, int nitems, FILE* fp, int init) //KLC
{
	static int count, ch;
	int readbytes, c, read;
	
	if (init)
	{
		count = ch = 0; 
	}
	
	readbytes = size * nitems;
	
	for (read = 0; read < readbytes; read++) 
	{
		if (count) 
		{
			*ptr++ = (unsigned char)ch;
			count--;
		}
		
		else 
		{
			c = getc(fp);
			if (c == EOF) break;
			
			if (c == RAS_RLE)
			{ 
				/* 0x80 */
				count = getc(fp);
				if (count == EOF) break;
				
				if (count < 0) count &= 0xff;
				if (count == 0) *ptr++ = c;
				else
				{
					if ((ch = getc(fp)) == EOF) break;
					*ptr++ = ch;
				}
			}
			else *ptr++ = c;
		}
	}
	return (read/size);
}

/*****************************/
static int SunRasError(char *st)
{
	fprintf(stderr, "LoadSunRas() - %s\n", st);
	return -1;
}

/*****************************/
static void flipl(unsigned char *p)
{
	unsigned char t;
	
	t = p[0];  p[0]=p[3];  p[3] = t;
	t = p[1];  p[1]=p[2];  p[2] = t;
}

/*-------------------------------------------------------------------
* CTReadPGM - read a slice from a portable greymap file
*-------------------------------------------------------------------
*/
static int CTReadPGM(FILE *fp, CTSlice slice, int x1, int y1, int x2, int y2)
{
	int c, i, j, width, height, max;
	unsigned char *line;
	char buf[1024];
	
	/* read for the magic number, skipping comments */
	while ((c=fgetc(fp)) == '#')
		fgets(buf, 1024, fp);

	ungetc(c, fp);
	if (fscanf(fp, "P%d\n", &c) != 1)
		return(-1);
	
	if (c != 2 && c != 5) 
	{
		/* file is not PGM, fail */
		return(-1);
	}
	
	if (c==2) 
	{
		
		while ((c=fgetc(fp)) == '#')
			fgets(buf, 1024, fp);
		ungetc(c, fp);
		
		if (fscanf(fp, "%d%d%d", &width, &height, &max) != 3)
			return(-1);
		
		if (x1==-1) 
		{
			slice->x1=slice->y1=x1=y1=0;
			slice->x2=x2=width-1;
			slice->y2=y2=height-1;
			slice->data = (unsigned short *)malloc(sizeof(unsigned short)*
				(x2-x1+1)*(y2-y1+1));
			memset(slice->data, 0, sizeof(short)*(x2-x1+1)*(y2-y1+1));
			CTSliceWidth(slice) = width;
			CTSliceHeight(slice) = height;
		}
		
		/* eat extra character */
		getc(fp);
		for (i=height-1; i>=0; i--)
		{
			for (j=0; j<width; j++) 
			{
				int val;
				
				fscanf(fp, "%d", &val);
				CTSliceShortData(slice, j, i) = val;
			}
		}

			CTSliceCompMinMaxD(slice);
			return(0);
	}
	else if (c==5) 
	{
		while ((c=fgetc(fp)) == '#')
			fgets(buf, 1024, fp);
		ungetc(c, fp);
		if (fscanf(fp, "%d%d%d", &width, &height, &max) != 3) 
		{
			printf("failed to read width/height/max\n");
			return(-1);
		}
		
		if (x1==-1) 
		{
			slice->x1=slice->y1=x1=y1=0;
			slice->x2=x2=width-1;
			slice->y2=y2=height-1;
			slice->cdata = (unsigned char *)malloc(sizeof(unsigned char)*
				(x2-x1+1)*(y2-y1+1));
			memset(slice->cdata, '\0', sizeof(char)*(x2-x1+1)*(y2-y1+1));
			CTSliceWidth(slice) = width;
			CTSliceHeight(slice) = height;
		}
		/* eat extra character */
		getc(fp);
		
		line = (unsigned char *)malloc(sizeof(char)*width);
		for (i=height-1; i>=0; i--) 
		{
			fread(line, 1, width, fp);
			if (y1 <= i && i <= y2)
			{
				for (j=0; j<width; j++)
				{
					if (x1 <= j && j <= x2)
						CTSliceCharData(slice, j, i) = line[j];
				}
			}
		}
		free(line);
		CTSliceCompMinMaxD(slice);
		return(0);
	}
	return(1);
}

/*-------------------------------------------------------------------
* CTReadVol - read a slice from a volume file
*-------------------------------------------------------------------
*/
static int CTReadVol(FILE *fp, CTSlice slice, int num, int x1, int y1,
                     int x2, int y2)
{
	
	/*
	CTSliceXSize(slice) = 256;
	CTSliceYSize(slice) = 256;
	*/
	CTSliceWidth(slice) = 256;
	CTSliceHeight(slice) = 256;
	
	if (x1==-1) 
	{
		slice->x1=slice->y1=x1=y1=0;
		slice->x2=x2=CTSliceXSize(slice)-1;
		slice->y2=y2=CTSliceYSize(slice)-1;
		slice->data=(unsigned short *)malloc(sizeof(short)*(x2-x1+1)*(y2-y1+1));
		memset(slice->data, '\0', sizeof(short)*(x2-x1+1)*(y2-y1+1));
	}
	
	fskip(fp, 256*256*2*(num-1));
	fread(slice->data, 2, 256*256, fp);
	CTSliceCompMinMaxD(slice);
	return(0);
}

/*-------------------------------------------------------------------
* CTReadSLC - read a slice from an SLC file
*-------------------------------------------------------------------
*/
static int CTReadSLC(FILE *fp, CTSlice slice, int num, int x1, int y1,
                     int x2, int y2)
{
	int i, remaining, current_value;
	unsigned int size;
//	unsigned int xsize, ysize, zsize, bits, magic;
//	float xunit, yunit, zunit;
	unsigned char *data_ptr;
	C_CompressionType compression;
//	C_UnitType unit;
//	C_DataOrigin data;
//	C_DataModification datamod;
	
	if (x1==-1) 
	{
		slice->x1=slice->y1=x1=y1=0;
		slice->x2=x2=CTSliceXSize(slice)-1;
		slice->y2=y2=CTSliceYSize(slice)-1;
		slice->cdata=(unsigned char *)malloc(sizeof(char)*(x2-x1+1)*(y2-y1+1));
		memset(slice->cdata, '\0', sizeof(char)*(x2-x1+1)*(y2-y1+1));
	}
	
	CTSliceWidth(slice) = x2-x1+1;
	CTSliceHeight(slice) = y2-y1+1;
	
	/* skip header */
	fskip(fp, sizeof(int)*5 + sizeof(float)*3 + sizeof(C_UnitType) +
		sizeof(C_DataOrigin) + sizeof(C_DataModification));
	
	/* get compression type */
	fread(&compression, sizeof(C_CompressionType), 1, fp);
	
	/* skip (num-1) slices */
	if (compression == C_NO_COMPRESSION) 
	{
		fskip(fp, CTVolXSize(slice->vol)*CTVolYSize(slice->vol)*(num-1));
	}
	else 
	{
		for (i=0; i<num-1; i++) 
		{
			fread(&size, sizeof(int), 1, fp);
			fskip(fp, size);
		}
	}
	
	/* read slice */
	switch (compression) 
	{
		case C_NO_COMPRESSION:
			fread(slice->cdata, sizeof(char),
				CTVolXSize(slice->vol)*CTVolYSize(slice->vol), fp);
			break;

		case C_RUN_LENGTH_ENCODE:
			fread(&size, sizeof(int), 1, fp);
			data_ptr = slice->cdata;
			while (1) 
			{
				current_value = getc(fp);
				if ( !(remaining = (current_value & 0x7f)) )
					break;
				if ( current_value & 0x80 )
					while (remaining--)
						*(data_ptr++) = getc(fp);
					else 
					{
						current_value = getc(fp);
						while (remaining--)
							*(data_ptr++) = current_value;
					}
			}
			break;

		default:
			fprintf(stderr, "Error in compression type of SLC file\n");
			return(-1);
	}
	
	CTSliceCompMinMaxD(slice);
	return(0);
}

/*-------------------------------------------------------------------
* CTReadH3D - read a slice from an H3D file
*-------------------------------------------------------------------
*/
static int CTReadH3D(FILE *fp, CTSlice slice, int num, int x1, int y1,
                     int x2, int y2)
{
	if (x1==-1) 
	{
		slice->x1=slice->y1=x1=y1=0;
		slice->x2=x2=CTSliceXSize(slice)-1;
		slice->y2=y2=CTSliceYSize(slice)-1;
		slice->fdata=(float *)malloc(sizeof(float)*(x2-x1+1)*(y2-y1+1));
		memset(slice->fdata, '\0', sizeof(float)*(x2-x1+1)*(y2-y1+1));
	}
	
	CTSliceWidth(slice) = x2-x1+1;
	CTSliceHeight(slice) = y2-y1+1;
	
	fskip(fp, sizeof(int)*3 + sizeof(double)*3 +
		CTVolXSize(slice->vol)*CTVolYSize(slice->vol)*(num-1)*sizeof(float));
	fread(slice->fdata, sizeof(float), CTVolXSize(slice->vol)*
		CTVolYSize(slice->vol), fp);
	
	CTSliceCompMinMaxD(slice);
	return(0);
}

/*-------------------------------------------------------------------
* CTOpenSliceNum - open a slice for reading
*-------------------------------------------------------------------
*/
FILE *CTOpenSliceNum(char *prefix, int num, char *suffix)
{
	int zeros;
	FILE *fp;
	char fname[200];
	
	for (zeros=0; zeros<=3; zeros++) 
	{
		if (zeros == 0)
			sprintf(fname, "%s%d%s", prefix, num, suffix);
		else if (zeros == 1)
			sprintf(fname, "%s0%d%s", prefix, num, suffix);
		else if (zeros == 2)
			sprintf(fname, "%s00%d%s", prefix, num, suffix);
		else if (zeros == 3)
			sprintf(fname, "%s000%d%s", prefix, num, suffix);
		//KLC
		//fzopen is used for compressed files ".z" or ".gz" files declared in utils\fileutil.h
		//if ((fp=fzopen(fname, "r")) != NULL)
		//   return(fp);
		if ((fp=fopen(fname, "r")) != NULL)
			return(fp);
	}
	return(NULL);
}

/*-------------------------------------------------------------------
* CTOpenFile - open a slice for reading
*-------------------------------------------------------------------
*/
FILE *CTOpenFile(CTVolume vol, int num)
{
	FILE *fp;
	char fname[200];
	
	fp = NULL;
	
	if (vol->voltype == CTSLC) 
	{
		sprintf(fname, "%s%s", vol->prefix, vol->suffix);
		//KLC
		//fzopen is used for compressed files ".z" or ".gz" files declared in utils\fileutil.h
		//if ((fp=fzopen(fname, "r")) != NULL)
		//   return(fp);
		if ((fp=fopen(fname, "r")) != NULL)
			return(fp);
	}
	else if (vol->voltype == CTVOL) 
	{
		sprintf(fname, "%s", vol->prefix);
		//KLC
		//fzopen is used for compressed files ".z" or ".gz" files declared in utils\fileutil.h
		//if ((fp=fzopen(fname, "r")) != NULL)
		//   return(fp);
		if ((fp=fopen(fname, "r")) != NULL)
			return(fp);
	}
	else if (vol->voltype == CTH3D) 
	{
		sprintf(fname, "%s%s", vol->prefix, vol->suffix);
		//KLC
		//fzopen is used for compressed files ".z" or ".gz" files declared in utils\fileutil.h
		//if ((fp=fzopen(fname, "r")) != NULL)
		//   return(fp);
		if ((fp=fopen(fname, "r")) != NULL)
			return(fp);
	}
	else
		fp = CTOpenSliceNum(vol->prefix, num, vol->suffix);
	
	return(fp);
}

/*-------------------------------------------------------------------
*  CTSliceRead  --  read a ct slice (or part of one)
*-------------------------------------------------------------------
*/
CTSlice CTSliceRead(CTVolume vol, int num, int x1, int y1, int x2,
					int y2)
{
	FILE *fp;
	CTSlice slice;
	int err;
//	int zeros;
//	char path[256], zfile[256];
	
	slice = (CTSlice)malloc(sizeof(CTSliceS));
	memset(slice, '\0', sizeof(CTSliceS));
	
	slice->num = num;
	slice->vol = vol;
	
	if ((fp = CTOpenFile(vol, num)) == NULL) 
	{
		free(slice);
		return(NULL);
	}
	
	if (x1 != -1) 
	{
		if (CTSliceShorts(slice)) 
		{
			slice->data = (unsigned short *)malloc(sizeof(unsigned short)*
				(x2-x1+1)*(y2-y1+1));
			memset(slice->data, '\0', sizeof(short)*(x2-x1+1)*(y2-y1+1));
		}
		else 
		{
			slice->cdata = (unsigned char *)malloc(sizeof(unsigned char)*
				(x2-x1+1)*(y2-y1+1));
			memset(slice->cdata, '\0', sizeof(char)*(x2-x1+1)*(y2-y1+1));
		}
		
		slice->x1 = x1;
		slice->y1 = y1;
		slice->x2 = x2;
		slice->y2 = y2;
		CTSliceWidth(slice) = (x2-x1+1);
		CTSliceHeight(slice) = (y2-y1+1);
	}
	slice->maxd = 0;
	
	switch (vol->voltype) 
	{
		case CTCTF:
			err = CTReadCTF(fp, slice, x1, y1, x2, y2);
			break;

		case CTRAS:
			err = CTReadRAS(fp, slice, x1, y1, x2, y2);
			break;

		case CTPGM:
			err = CTReadPGM(fp, slice, x1, y1, x2, y2);
			break;

		case CTVOL:
			err = CTReadVol(fp, slice, num, x1, y1, x2, y2);
			break;

		case CTSLC:
			err = CTReadSLC(fp, slice, num, x1, y1, x2, y2);
			break;

		case CTH3D:
			err = CTReadH3D(fp, slice, num, x1, y1, x2, y2);
			break;

		default:
			err = 1;
			break;
	}
	
	//fzclose(fp); //KLC
	fclose(fp);
	
	if (err)
	{
		if (CTSliceShorts(slice) && slice->data != NULL)
			free(slice->data);
		else if (CTSliceChars(slice) && slice->cdata != NULL)
			free(slice->cdata);
		else if (CTSliceFloats(slice) && slice->fdata != NULL)
			free(slice->fdata);
		free(slice);
		return(NULL);
	}
	else
		return(slice);
}

/*-------------------------------------------------------------------
* CTWriteCTF - write a CTF file
*-------------------------------------------------------------------
*/
static int CTWriteCTF(char *path, CTSlice slice)
{
	FILE *fp;
	unsigned short delta[512];
	int i, j, x, left, right, offsetx, offsety;
	
	if ((fp=fopen(path, "w")) == NULL)
		return(0);
	
	/* write header */
	for (i=0; i<2048; i++)
		fputc(255, fp);
	
	if (CTSliceWidth(slice) < 512)
		offsetx = (512-CTSliceWidth(slice))/2;
	else
		offsetx = 512;

	if (CTSliceHeight(slice) < 512)
		offsety = (512-CTSliceHeight(slice))/2;
	else
		offsety = 512;
	
	/* write the delta's */
	memset(delta, '\0', sizeof(short)*512);
	for (x=CTSliceMinX(slice); x<=CTSliceMaxX(slice)&&x-CTSliceMinX(slice)<512;	x++) 
	{
		for (left=CTSliceMinY(slice); left<CTSliceYSize(slice)/2; left++)
		{
			if (CTSliceData(slice, x, left) != 0)
				break;
		}
		
		for (right=CTSliceMaxY(slice); right>CTSliceYSize(slice)/2; right--)
		{
			if (CTSliceData(slice, x, right) != 0)
				break;
		}
		
		if (CTSliceYSize(slice)/2-left > right-CTSliceYSize(slice)/2)
			delta[x] = (CTSliceYSize(slice)/2-left);
		else
			delta[x] = (right-CTSliceYSize(slice)/2);
	}
	fwrite(delta, 2, 512, fp);
	
	/* write out the slice */
	for (i=0; i<offsetx/2; i++)
	{
		for (i=CTSliceMinX(slice); i<=CTSliceMaxX(slice); i++) 
		{
			for (j=CTSliceYSize(slice)/2-delta[i]; j<CTSliceYSize(slice)/2+delta[i];
			j++)
				fwrite(&CTSliceShortData(slice, i, j), sizeof(short), 1, fp);
		}
	}
	fclose(fp);
	return(1);
}

/*-------------------------------------------------------------------
* CTWritePGM - write a slice in portable greymap format
*-------------------------------------------------------------------
*/
static int CTWritePGM(char *path, CTSlice slice)
{
	FILE *fp;
	int i, j;
//	int c, width, height, max;
//	unsigned char *line;
	
	if ((fp=fopen(path, "w")) == NULL)
		return(0);
	
	if (CTSliceMaxD(slice) <= 255) 
	{
		fprintf(fp, "P5\n%d %d\n%d\n", CTSliceXSize(slice), CTSliceYSize(slice),
			CTSliceMaxD(slice));
		if (CTSliceChars(slice)) 
		{
			for (i=CTSliceYSize(slice)-1; i>=0; i--)
			{
				for (j=0; j<CTSliceXSize(slice); j++)
					fputc(CTSliceCharData(slice, j, i), fp);
			}
		}
		else if (CTSliceShorts(slice)) 
		{
			for (i=CTSliceYSize(slice)-1; i>=0; i--)
			{
				for (j=0; j<CTSliceXSize(slice); j++)
					fputc(CTSliceShortData(slice, j, i), fp);
			}
		}
	}
	else
	{
		fprintf(fp, "P2\n%d %d\n%d\n", CTSliceXSize(slice), CTSliceYSize(slice),
			CTSliceMaxD(slice));
		for (i=CTSliceYSize(slice)-1; i>=0; i--)
		{
			for (j=0; j<CTSliceXSize(slice); j++)
				fputc(CTSliceShortData(slice, j, i), fp);
		}
	}
	
	fclose(fp);
	return(1);
}

/*-------------------------------------------------------------------
* CTWriteRAS - write a slice in sun rasterfile format
*-------------------------------------------------------------------
*/
static int CTWriteRAS(char *path, CTSlice slice)
{
	return(0);
}

/*-------------------------------------------------------------------
* CTWriteH3D - write a slice in H3D format
*-------------------------------------------------------------------
*/
static int CTWriteH3D(CTSlice slice)
{
	fwrite(slice->fdata, sizeof(float), CTVolXSize(slice->vol)*
		CTVolYSize(slice->vol), slice->vol->outfp);
	return 0;
}

/*-------------------------------------------------------------------
* CTWriteSLC - write a slice in SLC format
*-------------------------------------------------------------------
*/
static int CTWriteSLC(CTSlice slice)
{
	switch (CTVolVoxelBits(slice->vol)) 
	{
		case CT_8BIT:
			fwrite(slice->cdata, sizeof(char), CTVolXSize(slice->vol)*
				CTVolYSize(slice->vol), slice->vol->outfp);
			break;

		case CT_16BIT:
			fwrite(slice->data, sizeof(short), CTVolXSize(slice->vol)*
				CTVolYSize(slice->vol), slice->vol->outfp);
			break;
	}
	return 0;
}

/*-------------------------------------------------------------------
* CTWriteVOL - write a slice in raw volume format
*-------------------------------------------------------------------
*/
static int CTWriteVOL(CTSlice slice)
{
	return(0);
}

/*-------------------------------------------------------------------
* CTSliceWrite - write a slice file
*-------------------------------------------------------------------
*/
int CTSliceWrite(char *path, CTSlice slice)
{
	switch (slice->vol->outtype) 
	{
		case CTCTF:
			return(CTWriteCTF(path, slice));

		case CTPGM:
			return(CTWritePGM(path, slice));

		case CTRAS:
			return(CTWriteRAS(path, slice));

		case CTH3D:
			return(CTWriteH3D(slice));

		case CTSLC:
			return(CTWriteSLC(slice));

		case CTVOL:
			return(CTWriteVOL(slice));

		default:
			return(0);
	}
}

/*-------------------------------------------------------------------
* CTSliceSetVals - set a block of values in a CT image
*-------------------------------------------------------------------
*/
void CTSliceSetVals(CTSlice slice, int x1, int y1, int x2, int y2,
                    int val, int z)
{
	int x, y;
	
	if (x1 > x2) 
	{
		x = x1;
		x1 = x2;
		x2 = x;
	}
	
	if (y1 > y2) 
	{
		y = y1;
		y1 = y2;
		y2 = y;
	}
	
	x1 = MAX(x1, slice->x1);
	y1 = MAX(y1, slice->y1);
	x2 = MIN(x2, slice->x2);
	y2 = MIN(y2, slice->y2);
	
	if (z) 
	{
		if (x1 <= x2) 
		{
			/* set the first line */
			if (CTSliceShorts(slice)) 
			{
				for (y=y1; y<=y2; y++)
					CTSliceShortData(slice, x1, y) = val;
				for (x1++; x1<=x2; x1++)
					memcpy(&CTSliceShortData(slice, x1, y1),&CTSliceShortData(slice, x1-1, y1),
						sizeof(short)*(y2-y1+1));
			}
			else 
			{
				for (y=y1; y<=y2; y++)
					CTSliceCharData(slice, x1, y) = val;
				for (x1++; x1<=x2; x1++)
					memcpy(&CTSliceCharData(slice, x1, y1),&CTSliceCharData(slice, x1-1, y1),
						sizeof(char)*(y2-y1+1));
			}
		}
	}
	else 
	{
		if (CTSliceShorts(slice)) 
		{
			for (x=x1; x<=x2; x++)
			{
				for (y=y1; y<=y2; y++)
					if (CTSliceShortData(slice, x, y) > 1000)
						CTSliceShortData(slice, x, y) = val;
					else
					CTSliceShortData(slice, x, y) = 0;
			}
		}
		else 
		{
			for (x=x1; x<=x2; x++)
			{
				for (y=y1; y<=y2; y++)
					if (CTSliceCharData(slice, x, y) > 1000)
						CTSliceCharData(slice, x, y) = val;
					else
						CTSliceCharData(slice, x, y) = 0;
			}
		}
	}
}

/*-------------------------------------------------------------------
* CTSliceSetVal - set a value in a CT image
*-------------------------------------------------------------------
*/
void CTSliceSetVal(CTSlice slice, int x, int y, float val)
{
	if (CTSliceShorts(slice))
		CTSliceShortData(slice, x, y) = (unsigned short) val;
	else if (CTSliceChars(slice))
		CTSliceCharData(slice, x, y) = (unsigned char) val;
	else if (CTSliceFloats(slice))
		CTSliceFloatData(slice, x, y) = val;
}

/*-------------------------------------------------------------------
* CTSliceFree - free a CTSlice
*-------------------------------------------------------------------
*/
CTSlice CTSliceFree(CTSlice slice)
{
	if (slice==NULL)
		return(NULL);

	if (slice->data != NULL && CTSliceShorts(slice))
		free(slice->data);

	if (slice->cdata != NULL && CTSliceChars(slice))
		free(slice->cdata);

	if (slice->fdata != NULL && CTSliceFloats(slice))
		free(slice->fdata);

	free(slice);
	return(NULL);
}

/*-------------------------------------------------------------------
* CTSliceCopy - copy a CTSlice
*-------------------------------------------------------------------
*/
CTSlice CTSliceCopy(CTSlice slice, int zero)
{
	CTSlice newSlice;
	
	if (slice==NULL)
		return(NULL);
	newSlice=(CTSlice)malloc(sizeof(CTSliceS));
	newSlice->vol = slice->vol;
	newSlice->num = slice->num;
	newSlice->x1 = slice->x1;
	newSlice->y1 = slice->y1;
	newSlice->x2 = slice->x2;
	newSlice->y2 = slice->y2;
	newSlice->width = slice->width;
	newSlice->height = slice->height;
	newSlice->mind = slice->mind;
	newSlice->maxd = slice->maxd;
	
	if (CTSliceShorts(newSlice)) 
	{
		newSlice->data = (unsigned short *)malloc(
			sizeof(short)*slice->width*slice->height);
		if (zero)
			memset(newSlice->data, '\0', sizeof(short)*slice->width*slice->height);
		else
			memcpy(newSlice->data, slice->data,
			sizeof(short)*slice->width*slice->height);
	}
	else if (CTSliceChars(newSlice)) 
	{
		newSlice->cdata = (unsigned char *)malloc(
			sizeof(char)*slice->width*slice->height);
		if (zero)
			memset(newSlice->cdata, '\0', sizeof(char)*slice->width*slice->height);
		else
			memcpy(newSlice->cdata, slice->cdata,
			sizeof(char)*slice->width*slice->height);
	}
	else 
	{
		newSlice->fdata = (float *)malloc(
			sizeof(float)*slice->width*slice->height);
		if (zero)
			memset(newSlice->fdata, '\0', sizeof(float)*slice->width*slice->height);
		else
			memcpy(newSlice->fdata, slice->fdata,
			sizeof(float)*slice->width*slice->height);
	}
	return(newSlice);
}

/*-------------------------------------------------------------------
* CTSliceCreate - create a CTSlice
*-------------------------------------------------------------------
*/
CTSlice CTSliceCreate(int width, int height, float mind, float maxd,
                      CTVolume vol, int num)
{
	CTSlice newSlice;
	
	newSlice=(CTSlice)malloc(sizeof(CTSliceS));
	newSlice->x1 = 0;
	newSlice->y1 = 0;
	newSlice->x2 = width-1;
	newSlice->y2 = height-1;
	newSlice->width = width;
	newSlice->height = height;
	newSlice->mind = mind;
	newSlice->maxd = maxd;
	newSlice->num = num;
	newSlice->vol = vol;
	
	if (CTSliceShorts(newSlice)) 
	{
		newSlice->data = (unsigned short *)malloc(sizeof(short)*width*height);
		newSlice->cdata = NULL; newSlice->fdata = NULL;
	}
	else if (CTSliceChars(newSlice)) 
	{
		newSlice->cdata = (unsigned char *)malloc(sizeof(char)*width*height);
		newSlice->data = NULL; newSlice->fdata = NULL;
	}
	else if (CTSliceFloats(newSlice)) 
	{
		newSlice->fdata = (float *)malloc(sizeof(float)*width*height);
		newSlice->cdata = NULL; newSlice->data = NULL;
	}
	return(newSlice);
}
