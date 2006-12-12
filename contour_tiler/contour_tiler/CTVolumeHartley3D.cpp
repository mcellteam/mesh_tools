/*
Function	:CTVolumeHartley3D
Parameters	:vol - volume to transform
Purpose	:Performs a 3D Fast Hartley Transform on the volume, vol.
Returns a new volume of type CTH3D with data values of
type CT_FLOAT.
Author	:Sean Vyain
Algorithm	:Electronic Letters 4th June 1992 Vol. 28 No. 12, pp1077-8
*/
#include <time.h>
#include <malloc.h>
#include <stdio.h>
#include <math.h>

#include "ct/ct.h"

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

CTVolume
CTVolumeHartley3D(CTVolume vol)
{
	CTVolume vol2;
	CTSlice slice;
	int l,m,n,L,M,N,i,j,k;
	
	#ifdef USE_PTIME
		int t1,t2;
	#endif
	
	float *d1,*d2,*X,*D,Z,*w,*w1,*w2;
	float *Mcas, *Lcos, *Lsin, *Ncos, *Nsin, max, min;
	double dL, dM, dN, s, c;
	
	/* Get source volume info. */
	L = CTVolXSize(vol);
	M = CTVolYSize(vol);
	N = CTVolNumSlices(vol);
	dL = (double)L;
	dM = (double)M;
	dN = (double)N;
	
	/* Initialize the new volume. */
	vol2 = CTVolumeCreate(CTH3D,CTVolPrefix(vol),".h3d",CTVolXSize(vol),
		CTVolYSize(vol),CTVolZSize(vol),CTVolXSep(vol),
		CTVolYSep(vol),CTVolZSep(vol),CTVolFirstSlice(vol),
		CTVolLastSlice(vol),1.0e10,CT_MSBF,CT_FLOAT);
	
	/* Open the new volume for writing. */
	CTVolumeOpen(vol2,CTH3D,CTVolPrefix(vol2),".h3d");
	
	/* Read in the volume data. */
	D = (float *)malloc(L*M*N*sizeof(float));

	switch (CTVolVoxelBits(vol)) 
	{
		
		/* Read in and convert CT_8BIT values to CT_FLOAT. */
		case CT_8BIT:
			for (n=0;n<N;n++) 
			{
				slice = CTSliceRead(vol,n+CTVolFirstSlice(vol),-1,-1,-1,-1);
				for (l=0;l<L;l++) 
				{
					for (m=0;m<M;m++) 
					{
						D[n+N*(l+L*m)] = (float)CTSliceCharData(slice,l,m);
					}
				}
				CTSliceFree(slice);
			}
			break;
			
			/* Read in and convert CT_16BIT values to CT_FLOAT. */
		case CT_16BIT:
			for (n=0;n<N;n++) 
			{
				slice = CTSliceRead(vol,n+CTVolFirstSlice(vol),-1,-1,-1,-1);
				for (l=0;l<L;l++) 
				{
					for (m=0;m<M;m++) 
					{
						D[n+N*(l+L*m)] = (float)CTSliceShortData(slice,l,m);
					}
				}
				CTSliceFree(slice);
			}
			break;
			
			/* Read in CT_FLOAT values. */
		case CT_FLOAT:
			for (n=0;n<N;n++) 
			{
				slice = CTSliceRead(vol,n+CTVolFirstSlice(vol),-1,-1,-1,-1);
				for (l=0;l<L;l++) 
				{
					for (m=0;m<M;m++) 
					{
						D[n+N*(l+L*m)] = CTSliceFloatData(slice,l,m);
					}
				}
				CTSliceFree(slice);
			}
			break;
	}
	
	/* Precompute necessary sin, cos, and cas function values. */
	Lcos = (float *)malloc(sizeof(float)*L*L);
	Lsin = (float *)malloc(sizeof(float)*L*L);
	Ncos = (float *)malloc(sizeof(float)*N*N);
	Nsin = (float *)malloc(sizeof(float)*N*N);
	Mcas = (float *)malloc(sizeof(float)*M*M);

	for (i=0;i<L*L;i++) 
	{
		s = sin(2.0*M_PI*(double)i/dL);
		c = cos(2.0*M_PI*(double)i/dL);
		Lcos[i] = (float)c;
		Lsin[i] = (float)s;
	}
	
	for (i=0;i<N*N;i++) 
	{
		s = sin(2.0*M_PI*(double)i/dN);
		c = cos(2.0*M_PI*(double)i/dN);
		Ncos[i] = (float)c;
		Nsin[i] = (float)s;
	}
	
	for (i=0;i<M*M;i++) 
	{
		s = sin(2.0*M_PI*(double)i/dM);
		c = cos(2.0*M_PI*(double)i/dM);
		Mcas[i] = (float)(s+c);
	}
	
	/* Algorithm is O(LLMN+LMMN+LMNN). */
	X = (float *)malloc(L*M*sizeof(float));
	w1 = (float *)malloc(L*M*sizeof(float));
	w2 = (float *)malloc(M*sizeof(float));
	
#ifdef USE_PTIME
	t1 = t2 = time(NULL);
#endif
	
	/*
    Transfrom the original volume and write the new volume one slice at
    a time, in order.
	*/
	for (k=0;k<N;k++) 
	{
		
		/* Initialize max and min values for each slice. */
		min = 1.0e10;
		max = -1.0e10;
		
		/* Compute w1(i,m,n) */
		d1 = D;
		w = w1;
		for (m=0;m<M;m++) 
		{
			for (l=0;l<L;l++) 
			{
				d2 = &D[N*((L-l)%L+L*((M-m)%M))];
				*w = 0.0;
				
				for (n=0;n<N;n++) 
				{
					*w += (*d1++)*Ncos[k*n]+(*d2++)*Nsin[k*n];
				}
				
				w++;
			}
		}
		
		/* Compute w2(i,j,n) */
		for (i=0;i<L;i++) 
		{
			w = w2;
			d1 = w1;
			for (m=0;m<M;m++) 
			{
				*w = 0.0;
				d2 = &w1[L*((M-m)%M)];
				
				for (l=0;l<L;l++) 
				{
					*w += (*d1++)*Lcos[i*l]+(*d2++)*Lsin[i*l];
				}
				
				w++;
			}
			
			/* Compute X(i,j,k) */
			for (j=0;j<M;j++) 
			{
				Z = 0.0;
				w = w2;
				
				for (m=0;m<M;m++) 
				{
					Z += (*w++)*Mcas[j*m];
				}
				
				if (Z>max) max = Z;
				if (Z<min) min = Z;
				X[i+L*j] = Z;
			}
		}
		
		/* Write the transformed slice. */
		slice = CTSliceCreate(L,M,min,max,vol2,k);

		for (i=0;i<L;i++) 
		{
			for (j=0;j<M;j++) 
			{
				CTSliceSetVal(slice,i,j,X[i+L*j]);
			}
		}
		CTSliceWrite(NULL,slice);
		CTSliceFree(slice);
		
#ifdef USE_PTIME
		printf("Slice %d/%d, ",k+1,N);
		ptime((int)((float)(time(NULL)-t1)/(float)(k+1)*(float)(N-k-1)));
		printf("\n");
#endif
		
	}
	
	/* Done writing to the new volume. */
	CTVolumeClose(vol2);
	
#ifdef USE_PTIME
	t2 = time(NULL);
	printf("(3D)Total time: ");
	ptime(t2-t1);
	printf("\n");
#endif
	
	/* Return the transformed volume. */
	return(vol2);
}
