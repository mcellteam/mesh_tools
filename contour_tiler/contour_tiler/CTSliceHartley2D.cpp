/*
Function	:CTSliceHartley2D
Parameters	:in - slice to transform
Purpose	:Applies a Fast Hartley Transform to the input slice, in.
A new slice is created that has the same number and is
from the same volume as the input slice.  This transformed
slice is returned.
Author	:Sean Vyain
Algorithm	:Electronic Letters 4th June 1992 Vol. 28 No. 12 pp1077-8
*/
#include <malloc.h>
#include <stdio.h>
#include "ct/ct.h" 

CTSlice
CTSliceHartley2D(CTSlice in)
{
	CTSlice slice;
	int m,n,i,j,M,N;
	float Z,min = 0.0f,max = 0.0f;
	register float *d1, *w1;
	extern float *H2D_Ncas, *H2D_Mcas, *D, *w;
	
	/* Get size of source slice. */
	M = CTSliceWidth(in);
	N = CTSliceHeight(in);
	slice = CTSliceCreate(M,N,min,max,CTSliceVol(in),CTSliceNum(in));
	
	/* Precompute necessary cas function values. */
	CTSliceH2DInit(M,N);
	
	/* Copy input slice to D array for faster reading. */
	d1 = D;

	for (j=0;j<N;j++) 
	{
		for (i=0;i<M;i++) 
		{
			(*d1++) = CTSliceFloatData(in,i,j);
		}
	}
	
	/* Algorithm, O(MMN+MNN). */
	min = 1.0e10;
	max = -1.0e10;
	
	/* Compute the transform of the array column by column. */
	for (i=0;i<M;i++) 
	{		
		/* Compute w(i,n) */
		d1 = D;
		w1 = w;
		for (n=0;n<N;n++) 
		{
			*w1 = 0.0;
			for (m=0;m<M;m++) {
				*w1 += (*d1++)*H2D_Mcas[i*m];
			}
			w1++;
		}
		
		/* Use w(i,n) to compute each of the entries in column i. */
		for (j=0;j<N;j++) 
		{
			Z = 0.0;
			w1 = w;
		
			for (n=0;n<N;n++) 
			{
				Z += (*w1++)*H2D_Ncas[j*n];
			}
			
			if (Z>max) max = Z;
			if (Z<min) min = Z;
			CTSliceSetVal(slice,i,j,Z);
		}
	}
	
	/* Return the transformed slice. */
	CTSliceMinD(slice) = min;
	CTSliceMaxD(slice) = max;
	return(slice);
}
