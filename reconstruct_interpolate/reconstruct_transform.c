#include <math.h>

double Xforward( int dim, double *a, double *b, double x, double y )
{
  double result;
  switch (dim) {
  case 1:                // translation only
    result = a[0] + x;
    break;
  case 2:
    result = a[0] + a[1]*x;
    break;
  case 3:                // affine transformation
    result = a[0] + a[1]*x + a[2]*y;
    break;
  case 4:
    result = a[0] + (a[1] + a[3]*y)*x + a[2]*y;
    break;
  case 5:
    result = a[0] + (a[1] + a[3]*y + a[4]*x)*x + a[2]*y;
    break;
  case 6:
    result = a[0] + (a[1] + a[3]*y + a[4]*x)*x + (a[2] + a[5]*y)*y;
    break;
  default:
    result = x;            // identity
  }
  return result;
}


double Yforward( int dim, double *a, double *b, double x, double y )
{
  double result;
  switch (dim) {
  case 1:                // translation only
    result = b[0] + y;
    break;
  case 2:
    result = b[0] + b[1]*y;
    break;
  case 3:                // affine transformation
    result = b[0] + b[1]*x + b[2]*y;
    break;
  case 4:
    result = b[0] + (b[1] + b[3]*y)*x + b[2]*y;
    break;
  case 5:
    result = b[0] + (b[1] + b[3]*y + b[4]*x)*x + b[2]*y;
    break;
  case 6:
    result = b[0] + (b[1] + b[3]*y + b[4]*x)*x + (b[2] + b[5]*y)*y;
    break;
  default:
    result = y;            // identity
  }
   return result;
}


void XYinverse( int dim, double *a, double *b, double *x, double *y )
{
  int i;
  double epsilon = 5e-10;
  double e,l,m,n,o,p,x0,y0,u0,v0,u,v;

  switch (dim) {
  case 0:
    break;              // no change when identity Nform
  case 1:
    *x = *x - a[0];
    *y = *y - b[0];
    break;
  case 2:
  case 3:
    u = *x - a[0];
    v = *y - b[0];
    p = a[1]*b[2] - a[2]*b[1];
    if ( fabs(p) > epsilon ) {
      *x = (b[2]*u - a[2]*v)/p;       // inverse of rotational part
      *y = (a[1]*v - b[1]*u)/p;
      }
    break;

  case 4:
  case 5:
  case 6:
    u = *x;             // (u,v) for which we want (x,y)
    v = *y;
    x0 = 0.0;           // initial guess of (x,y)
    y0 = 0.0;
    u0 = Xforward(dim, a, b, x0,y0);          //  get forward tform of initial guess
    v0 = Yforward(dim, a, b, x0,y0);
    i = 0;              // allow no more than 10 iterations
    e = 1.0;            // to reduce error to this limit
    while ( (e > epsilon) && (i<10) ) {
      i++;
      l = a[1] + a[3]*y0 + 2.0*a[4]*x0; // compute Jacobian
      m = a[2] + a[3]*x0 + 2.0*a[5]*y0;
      n = b[1] + b[3]*y0 + 2.0*b[4]*x0;
      o = b[2] + b[3]*x0 + 2.0*b[5]*y0;
      p = l*o - m*n;            // determinant for inverse
      if ( fabs(p) > epsilon ) {
        x0 += (o*(u-u0) - m*(v-v0))/p;  // inverse of Jacobian
        y0 += (l*(v-v0) - n*(u-u0))/p;  // and use to increment (x0,y0)
        }
      else {
        x0 += l*(u-u0) + n*(v-v0);    // try Jacobian transpose instead
        y0 += m*(u-u0) + o*(v-v0);
        }
      u0 = Xforward(dim, a, b, x0,y0);            //  get forward tform of current guess
      v0 = Yforward(dim, a, b, x0,y0);
      e = fabs(u-u0) + fabs(v-v0);    // compute closeness to goal
      }
    *x = x0;
    *y = y0;            // return final estimate of (x,y)
  } // end switch
}


