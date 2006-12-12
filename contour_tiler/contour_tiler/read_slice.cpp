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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include "common.h"
#include "myutil.h"
#include "gen_data.h"
#include "read_slice.h"

#define MAX_SLICE_NUM 1000
extern int Begin_slice_num;
extern int MATH_FUNC;
extern Interpo_struct Istru;
extern Name_struct Nstru;
static int *Posx_ary,*Posy_ary;
CTSlice Slice_ary[MAX_SLICE_NUM];
CTVolume Vol_in; /* a pointer to a structure */

extern int STATIC_CHQ;

void my_clear_read_slice()
{
	Posx_ary = NULL;
	Posy_ary = NULL;

	slicefunction(STATIC_CHQ, STATIC_CHQ);
}

void make_pos_table(CTSlice slice) /* compute the cube position */
{
    int i,y,dimy,dimx,x;
	
    y = CTSliceMinY(slice);
    x = CTSliceMinX(slice);
    dimx=CTSliceWidth(slice)/Istru.xcube_vox;
    dimy=CTSliceHeight(slice)/Istru.ycube_vox;
    Posy_ary = (int *)mymalloc(dimy*sizeof(int));
    Posx_ary = (int *)mymalloc(dimx*sizeof(int));

    for(i=0; i<dimy; i++) 
	{ 
        Posy_ary[i]=y+i*Istru.ycube_vox;
    }
    
	for(i=0; i<dimx; i++) 
	{ 
        Posx_ary[i]=x+i*Istru.xcube_vox;
    }
}


void pass_volume(CTVolume vol_in)
{
    Vol_in = vol_in;
}

float valuefunction(int x, int y, int z)
{
    float val;
    
	if(Slice_ary[z]==NULL)
	{ /* This slice is not read yet */
		fprintf(stderr,"Error! slice %d is not read yet\n",z);
		exit(1);
    }
    
	val=(float)CTSliceData(Slice_ary[z], Posx_ary[x], Posx_ary[y]);
    return val;
}

void slicefunction(int z, int needs)
{
    static int first=1;
    
	if((z+Begin_slice_num)>MAX_SLICE_NUM) 
	{
		fprintf(stderr,"slice # =%d > maximum # %d\n",(z+Begin_slice_num),
			MAX_SLICE_NUM);
		exit(0);
    }
	/*
    if((z+Begin_slice_num)<Istru.beg || (z+Istru.beg) > Istru.end) {
	fprintf(stderr,"Internal ERROR! slice %d\n",(z+Istru.beg));
	fprintf(stderr,"z=%d, beg=%d, end=%d\n",z,Istru.beg,Istru.end);
	exit(0);
    }
	*/
    
	if(needs==0)
	{
		CTSliceFree(Slice_ary[z]);
		Slice_ary[z] = NULL;
    }
    else 
    {
		if(Slice_ary[z] !=NULL) 
		{
			fprintf(stderr,"Warning: slicefunction() slice %d is read\n");
			return;
		}
        if(MATH_FUNC==0)
		{
            fprintf(stderr,"N reading slice # %d ..",(z+Begin_slice_num));
            //Slice_ary[z] = CTSliceRead(Vol_in, (z+Begin_slice_num), Istru.x1,Istru.y1, Istru.x2, Istru.y2);//KLC
			Slice_ary[z] = myCTSliceRead(Vol_in, (z+Begin_slice_num), Istru.x1,Istru.y1, Istru.x2, Istru.y2); 
        }
        else if(MATH_FUNC==4)
		{
            fprintf(stderr,"generating input slice # %d ..",(z+Begin_slice_num));
            Slice_ary[z] = specialCTSliceRead(Vol_in, (z+Begin_slice_num),
				Istru.x1,Istru.y1, Istru.x2, Istru.y2);
		}
        else 
		{
            fprintf(stderr,"generating input slice # %d ..",(z+Begin_slice_num));
            Slice_ary[z] = myCTSliceRead(Vol_in, (z+Begin_slice_num), Istru.x1, 
				Istru.y1, Istru.x2, Istru.y2);
        }
        fprintf(stderr," Done\n");
		
		if(Slice_ary[z]==NULL)
		{
			fprintf(stderr,"ERROR in reading slice %d\n",(z+Begin_slice_num));
			exit(1);
		}
		if(first)
		{
			first=0;
			make_pos_table(Slice_ary[z]);
		}
    }
}


