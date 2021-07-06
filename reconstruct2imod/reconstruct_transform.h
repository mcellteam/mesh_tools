#ifndef RECONSTRUCT_TRANSFORM_H
#define RECONSTRUCT_TRANSFORM_H

double Xforward( int dim, double *a, double *b, double x, double y );
double Yforward( int dim, double *a, double *b, double x, double y );
void XYinverse( int dim, double *a, double *b, double *x, double *y );

#endif
