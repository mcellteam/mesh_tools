/*
Function	:CTSliceH2DInit
Parameters	:M - horizontal size
N - vertical size
Purpose	:Pre-computes the values of the Mcas and Ncas arrays.  This
saves time if several calls are made to CTSliceHartley2D
with same-size slices.
Author	:Sean Vyain
*/

#include <stdio.h>
#include <malloc.h>
#include <math.h>

#include "ct/ct.h" //klc

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

int H2D_M = -1;
int H2D_N = -1;
float *H2D_Mcas=NULL;
float *H2D_Ncas=NULL;
float *D=NULL;
float *X=NULL;
float *w=NULL;
extern int STATIC_CHQ;

void my_clear_CTSliceH2DInit()
{
	CTSliceH2DInit(STATIC_CHQ, STATIC_CHQ);

}


void
CTSliceH2DInit(int M, int N)
{
	static int i;
	static double s,c,dM,dN;
	
	/* Don't do any work if the arrays have already been initialized. */
	if ((H2D_M==M) && (H2D_N==N)) return;
	
	/* Initialize the Mcas and Ncas arrays, given M and N. */
	H2D_M = M;
	H2D_N = N;
	dM = (double)M;
	dN = (double)N;
	free(H2D_Mcas);
	free(H2D_Ncas);
	free(D);
	free(X);
	free(w);
	H2D_Mcas = (float *)malloc(sizeof(float)*M*M);
	H2D_Ncas = (float *)malloc(sizeof(float)*N*N);
	D = (float *)malloc(sizeof(float)*M*N);
	X = (float *)malloc(sizeof(float)*M*N);
	w = (float *)malloc(sizeof(float)*N);

	for (i=0;i<M*M;i++) 
	{
		s = sin(2.0*M_PI*(double)i/dM);
		c = cos(2.0*M_PI*(double)i/dM);
		H2D_Mcas[i] = (float)(s+c);
	}
	
	for (i=0;i<N*N;i++) 
	{
		s = sin(2.0*M_PI*(double)i/dN);
		c = cos(2.0*M_PI*(double)i/dN);
		H2D_Ncas[i] = (float)(s+c);
	}
}
