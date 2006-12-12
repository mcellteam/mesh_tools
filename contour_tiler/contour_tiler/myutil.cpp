/*****************************************************************************/
/*****************************************************************************/
/**                                                                         **/
/** This SHASTRA software is not in the Public Domain. It is distributed on **/
/** a person to person basis, solely for educational use and permission is  **/
/** NOT granted for its transfer to anyone or for its use in any commercial **/
/** product.  There is NO warranty on the available software and neither    **/
/** Purdue University nor the Applied Algebra and Geometry group directed   **/
/** by C.  Bajaj accept responsibility for the consequences of its use.     **/
/**                                                                         **/
/*****************************************************************************/
/*****************************************************************************/

/* --------------------------------------------------------------
myutil.c
some utility functions. 
----------------------------------------------------------------- */
/*
#define MAX(x, y) ((x)>(y)?(x):(y))
#define MIN(x, y) ((x)<(y)?(x):(y))
*/


#define MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <memory.h>
#include <string.h>
#include <limits.h>

#include "common.h"
#include "myutil.h"


char *mymalloc(int size)
{
    char *p;
    if(size==0)size=4;

    if((p=(char *)malloc(size))==NULL)
    {
		fprintf(stderr,"malloc error, can not alocate %d bytes\n", size);
		exit(1);
    }
    memset(p,0,size);
    return p;
}

char *mycalloc(int size)
{
    char *p;
    p=(char *)mymalloc(size);
    memset(p,0,size);
    return p;
}

void malloc_copy(char **str1, char *str2)
{
    int len;
    len=strlen(str2)+1;
    
	if((*str1=(char*)malloc(len))==NULL)
    {
        fprintf(stderr,"malloc_copy(), error in malloc\n");
        exit(1);
    }
    strcpy(*str1,str2);
}



void print_slice_structure(CTSlice p)
{
    printf("\nprintf_slice_structure(0x%X):\n",p);
    if(p==NULL) return;
    
	printf("*CTVolume=0X%X, sh *data=0X%X, ch *cdata=0X%X, fl *fdata=0X%X\n",
        (int*)p->vol, (int *)p->data, (int*)p->cdata, (int*)p->fdata);
    
	printf("num=%d, x1=%d, y1=%d, x2=%d, y2=%d, width=%d, height=%d\n",
        p->num, p->x1, p->y1, p->x2, p->y2, p->width, p->height);
    
	printf("mind=%f, maxd=%f, ratio=%f\n\n",p->mind,p->maxd,p->ratio);
}

void print_volume_structure(CTVolume p)
{
    printf("\nprintf_volume_structure(0x%X):\n",p);
    if(p==NULL) return;
    
	printf("CTVolType=%d, prefix=%s, suffix=%s\n",
        p->voltype, p->prefix, p->suffix);
    
	printf("xdata=%d ydata=%d zdata=%d\n", p->xdata, p->ydata, p->zdata);
    
	printf("xunit=%f, yunit=%f zunit=%f\n",p->xunit,p->yunit,p->zunit);
    
	printf("numslice=%d, first=%d, last=%d\n",p->numslice,p->first,p->last);
    
	if(p->byteorder==CT_MSBF)printf("byteorder=CT_MSBF ");
    else if(p->byteorder==CT_LSBF)printf("byteorder=CT_LSBF ");
    else printf("byteorder=unknown ");
    if(p->voxelbits==CT_8BIT)printf("voxelbits=CT_8BIT\n");
    else if(p->voxelbits==CT_16BIT)printf("voxelbits=CT_16BIT\n");
    else if(p->voxelbits==CT_FLOAT)printf("voxelbits=CT_FLOAT\n");
    else printf("voxelbits=unknown\n");
}


CTVolume  InitVolume(int type, char *prefix, char *suffix, 
					 int first, int last, double zunits)
					 /* CTVolType    type;  */
{
    char      file[320];          /* 256+64, ct.h  */
#ifndef i860
    double noise =0.0;

    switch (type) 
	{
		case CTCTF:
			return ( CTVolumeInitCTF(prefix,suffix,first,last,
				zunits, noise, CT_MSBF,CT_16BIT) );
		case CTRAS:
			return ( CTVolumeInitRAS(prefix, suffix, first, last, zunits,noise) );
		case CTPGM:
			return ( CTVolumeInitPGM(prefix, suffix, first, last, zunits, noise) );
		case CTVOL:
			sprintf(file, "%s%s", prefix, suffix);
			return ( CTVolumeInitVOL(file, first, last, zunits, noise, CT_MSBF, CT_16BIT) );
		case CTSLC:
			sprintf(file, "%s%s", prefix, suffix);
			return ( CTVolumeInitSLC(file,noise) );
    }
#else
    switch (type) 
	{
		case CTCTF:
			return ( CTVolumeInitCTF(prefix,suffix,first,last,
				zunits,  CT_MSBF,CT_16BIT) );
			break;
		case CTRAS:
			return ( CTVolumeInitRAS(prefix, suffix, first, last, zunits) );
			break;
		case CTPGM:
			return ( CTVolumeInitPGM(prefix, suffix, first, last, zunits) );
			break;
		case CTVOL:
			sprintf(file, "%s%s", prefix, suffix);
			return ( CTVolumeInitVOL(file, first, last, zunits, CT_MSBF, CT_16BIT) );
			break;
		case CTSLC:
			sprintf(file, "%s%s", prefix, suffix);
			return ( CTVolumeInitSLC(file) );
		break;
    }
#endif
    return NULL;
}


