/* the file f2c/libF77/r_sign.c */

#include "f2c.h"

double r_sign(a,b)
real *a, *b;
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}
