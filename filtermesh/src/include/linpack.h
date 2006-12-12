
#if defined(ANSI)

void sqrdc_(double* a, int* m, int* m2, int* n,
	    double* qraux, int* jpvt, double* work, int* pivot);
/* a_fortran[n][m], qraux[n], jpvt[n], work[n] */

void sqrsl_(double* a, int* m, int* m2, int* n,
	    double* qraux, double* rh,
	    double*, double* qty, double* b, double*, double *,
	    int* mode, int* info);
/* a_fortran[n][m], qraux[n], rh[m], qty[n], b[n]
   mode==100 do linear least squares, return b */

#else

void sqrdc_();
void sqrsl_();

#endif
