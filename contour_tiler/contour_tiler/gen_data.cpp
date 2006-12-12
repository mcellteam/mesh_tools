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

/* -----------------------------------------------------------------------
gen_data.c
It generates a perfect sphere as the testing data. It replace the
CTSliceRead(), so no real file is written.

-------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include "common.h"
#include "myutil.h"
#include "gen_data.h"

#define XDIM 256 
#define YDIM 256 
#define ZDIM 32 
#define XDIM2 (XDIM/2) 
#define YDIM2 (YDIM/2) 
#define ZDIM2 (ZDIM/2) 
/* generate a sphere as slice data */

extern int MATH_FUNC;
extern int STATIC_CHQ;

void my_clear_gen_data()
{
	specialCTSliceRead(NULL, STATIC_CHQ, STATIC_CHQ, STATIC_CHQ, STATIC_CHQ, STATIC_CHQ);	
}


CTVolume myInitVolume(int type, char *prefix, char *suffix, 
					  int beg, int end, double detz)
{
    CTVolume vol;
    vol=(CTVolume)mymalloc(sizeof(*vol));
    vol->voltype=CTPGM;
    malloc_copy(&(vol->prefix),prefix);
    malloc_copy(&(vol->suffix),suffix);
    vol->xdata=XDIM;
    vol->ydata=YDIM;
    vol->zdata=end-beg+1;
    vol->xunit=1.0;
    vol->yunit=1.0;
    vol->zunit=(double)detz;
    vol->numslice=end-beg+1;
    vol->first=beg;
    vol->last=end;
    vol->byteorder= CT_MSBF;
    vol->voxelbits=CT_8BIT;
    return vol;
}
CTVolume myInitVolume16bits(int type, char *prefix, char *suffix, 
							int beg, int end, double detz)
{
    char *mymalloc(int);
    CTVolume vol;
    vol=(CTVolume)mymalloc(sizeof(*vol));
    vol->voltype=CTPGM;
    malloc_copy(&(vol->prefix),prefix);
    malloc_copy(&(vol->suffix),suffix);
    vol->xdata=XDIM;
    vol->ydata=YDIM;
    vol->zdata=end-beg+1;
    vol->xunit=1.0;
    vol->yunit=1.0;
    vol->zunit=(double)detz;
    vol->numslice=end-beg+1;
    vol->first=beg;
    vol->last=end;
    vol->byteorder= CT_MSBF;
    vol->voxelbits=CT_16BIT;
    return vol;
}

/* generate a sphere */
CTSlice myCTSliceRead(CTVolume vol, int num, int x1, int y1, int x2,
					  int y2)
{
    int i,j,cnt;
    double x,y,z,d;
    CTSlice slice;
    char *mymalloc(int);
    unsigned char *buf=NULL;
    unsigned short *sbuf=NULL;
	
    slice=(CTSlice)mycalloc(sizeof(*slice));
    slice->num = num;
    slice->vol = vol;

    if(x1== -1 || x2 == -1) 
	{
		x1=0; x2=vol->xdata-1; 
	}

    if(y1== -1 || y2 == -1)
	{
		y1=0; y2=vol->ydata-1; 
	}
    
	slice->x1 = x1;
    slice->y1 = y1;
    slice->x2 = x2;
    slice->y2 = y2;
    CTSliceWidth(slice) = (x2-x1+1);
    CTSliceHeight(slice) = (y2-y1+1);
    slice->maxd = 0;
    
	if (CTSliceShorts(slice)) 
	{
		slice->data = 
			(unsigned short *)mycalloc(sizeof(short)*(x2-x1+1)*(y2-y1+1));
		sbuf = slice->data;
    } 
	else 
	{
		slice->cdata = 
			(unsigned char *)malloc(sizeof(unsigned char)*(x2-x1+1)*(y2-y1+1));
		buf = slice->cdata;
    }
    
	cnt=0;
    z=(double)(num- ZDIM2) * vol->zunit;
    z=z*z;
    
	for(i=x1; i<=x2; i++) 
	{
		for(j=y1; j<=y2; j++) 
		{
			x=(i-XDIM2);
			y=(j-YDIM2);
			switch (MATH_FUNC)
			{
				case 1:
					d=x*x + y*y+ z;
					d=sqrt(d);
					break;
				case 2:
					d=x*x + y*y;
					d=sqrt(d);
					break;
				case 3:
					d=x*x + y*y-z;
					if(d<0)d=0;
					d=sqrt(d);
					break;
			}

			if(sbuf) 
			{
				d *= 256.0;
				if(d>65535)sbuf[cnt]=65535;
				else sbuf[cnt]= (unsigned short)d;
			}
			else 
			{
				if(d>255)buf[cnt]=255;
				else buf[cnt]=(unsigned char)d;
			}
			cnt++;
			
		}
    }
	
    return(slice);
}


void generate_ellipse(unsigned short *sbuf, int width, int height, int origx,
					  int origy, double a, double b)
{
    int i,j;
    int cx,cy,val,ref;
    if(width<=10 || height <=10)
		fprintf(stderr,"**Warning, width=%d, height=%d\n",width,height);
    cx=origx+width/2;
    cy=origy+height/2;
    ref= (int) (a*(width*width)/4 + b*(height*height)/4);

    for(j=origy; j<origy+height; j++)
	{
		for(i=origx; i<origx+width; i++)
		{
			val=(int)(a*((i-cx)*(i-cx))+b*((j-cy)*(j-cy)));
			val=ref-val;
			if(val<0)	val=0;
			sbuf[j*XDIM+i]=(unsigned short)val;
		}
    }
}


CTSlice specialCTSliceRead(CTVolume vol, int num, int x1, int y1, int x2,
						   int y2)
{
#define MAX_SLICE 100
    static int cnt=0,pre_num=-1;
    static int index[MAX_SLICE]={5,1};
	
    static double a[MAX_SLICE][10]=
	{
		{0.1, 0.12, 0.20, 0.18, 0.4, 0.0},
		{0.030, 0.0}
	};
		/*
		{0.040, 0.0}};
		*/
	static double b[MAX_SLICE][10]=
	{
		{0.10, 0.20, 0.10, 0.18, 0.6, 0.0},
		{0.05, 0.0}
	};
		
	static int width[MAX_SLICE][10]=
	{
		{80,80,80,80, 30,0},
		{140,0}
	};
	
	static int height[MAX_SLICE][10]=
	{
		{80,80,80,80,40,0},
		{140,0}
	};
			
	static int origx[MAX_SLICE][10]=
	{
		{40,40,140,140,100,0},
		{50,0}
	};

	static int origy[MAX_SLICE][10]=
	{
		{40,140,140, 40, 100,0},
		{30,0}
	};
		
		int i;
		CTSlice slice;
		char *mymalloc(int);
		unsigned short *sbuf=NULL;
		
		if(pre_num==num)cnt--;
	
		if(cnt>=MAX_SLICE)
		{
			fprintf(stderr,"special_CTSliceRead(), only %d slice allowed\n",
				MAX_SLICE);
			exit(1);
		}
		
		slice=(CTSlice)mycalloc(sizeof(*slice));
		slice->num = num;
		slice->vol = vol;
		
		if(x1== -1 || x2 == -1)
		{
			x1=0; x2=vol->xdata-1; 
		}
		
		if(y1== -1 || y2 == -1) 
		{
			y1=0; y2=vol->ydata-1; 
		}
		
		slice->x1 = x1;
		slice->y1 = y1;
		slice->x2 = x2;
		slice->y2 = y2;
		CTSliceWidth(slice) = (x2-x1+1);
		CTSliceHeight(slice) = (y2-y1+1);
		slice->maxd = 0;
		
		if (CTSliceShorts(slice)) 
		{
			slice->data = 
				(unsigned short *)mycalloc(sizeof(short)*(x2-x1+1)*(y2-y1+1));
			sbuf = slice->data;
		}
		else 
		{
			fprintf(stderr,"Internal Error, must use myInitVolume16bits()\n");
			exit(1);
		}
		
		for(i=0; i<index[cnt]; i++)
			generate_ellipse(sbuf, width[cnt][i], height[cnt][i], origx[cnt][i],
			origy[cnt][i], a[cnt][i], b[cnt][i]);
		pre_num=num;
		cnt++;
		return(slice);
}


