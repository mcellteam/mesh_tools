#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <malloc.h>

#include "common.h"
#include "myutil.h"
#include "filter.h"

static float Mbuf[10];

extern int STATIC_CHQ;

void my_clear_filter()
{
	int i;

	for (i=0; i<10; i++)
		Mbuf[i] =0.0;
}
static int float_compare(float *i, float *j)
{
    if(*i>*j)return 1;
    if(*i<*j)return -1;
    return 0;
}

static int float_compare2(const void *i, const void *j)
{
	return float_compare((float*) i, (float*) j);
}

void median_filter33(float *buf, int dimx, int dimy)
{
    int i,j,k,k1,area,dimx1,dimy1;
    float *buf1;
	
    area=dimx*dimy;
    dimx1=dimx-1;
    dimy1=dimy-1;
    buf1=(float*)mymalloc(sizeof(*buf1)*area);
    /* set the boundary to be 0 */

    for(i=0; i<dimx; i++)buf1[i]=0.0;
    
	for(i=area-1; i>=area-dimx; i--)buf1[i]=0.0;
    
	for(i=0; i<area; i=i+dimx)buf1[i]=0.0;
    
	for(i=dimx-1; i<area; i=i+dimx)buf1[i]=0.0;
	
    for(j=1; j<dimy1; j++)
	{
		for(i=1; i<dimx1; i++)
		{
			k=j*dimx+i;
			k1=k-dimx-1;
			Mbuf[0]=buf[k1++];
			Mbuf[1]=buf[k1++];
			Mbuf[2]=buf[k1];
			k1=k-1;
			Mbuf[3]=buf[k1++];
			Mbuf[4]=buf[k1++];
			Mbuf[5]=buf[k1];
			k1=k-1+dimx;
			Mbuf[6]=buf[k1++];
			Mbuf[7]=buf[k1++];
			Mbuf[8]=buf[k1];
			qsort(Mbuf, 9, sizeof(float), float_compare2);
			buf1[k]=Mbuf[4];
		}
    }
    memcpy(buf, buf1, area*sizeof(*buf));
    free(buf1);
}

void lowpass_filter33(float *buf, int dimx, int dimy)
{
    int i,j,k,area,dimx1,dimy1;
    float *buf1,sum;
	
    area=dimx*dimy;
    dimx1=dimx-1;
    dimy1=dimy-1;
    buf1=(float*)mymalloc(sizeof(*buf1)*area);
    /* set the boundary to be 0 */

    for(i=0; i<dimx*2; i++)buf1[i]=0.0;
    
	for(i=area-1; i>=area-(2*dimx); i--)buf1[i]=0.0;
    
	for(i=0; i<area; i=i+dimx){buf1[i]=0.0; buf1[i+1]=0.0;}
    
	for(i=dimx-1; i<area; i=i+dimx)buf1[i]=0.0;
	
	/* this is temporary for freddy which has noise on the boundary */
    for(j=2; j<dimy1-1; j++)
	{
		for(i=2; i<dimx1; i++)
		{
			k=j*dimx+i;
			sum=(buf[k]*2)+buf[k-1]+buf[k+1]+buf[k-dimx]+buf[k+dimx];
			/*
			sum=sum+buf[k-dimx-1]+buf[k-dimx+1]+buf[k+dimx+1]+buf[k+dimx-1];
			*/
			buf1[k]=sum/6.0f;
		}
    }
    memcpy(buf, buf1, area*sizeof(*buf));
    free(buf1);
}



