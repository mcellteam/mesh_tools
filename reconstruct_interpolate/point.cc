#include <iostream>
#include <stdlib.h>
#include <string>
#include <string.h>

#include "reconstruct_transform.h"
#include "point.h"

using std::cout;
using std::endl;

Point::Point (void)
  :x(0),y(0)
{
}

Point::Point (double xval, double yval)
  :x(xval),y(yval)
{
}

Point::Point (char const * str, int const dim, double * const t)
  :x(0),y(0)
{
  char val[80];
  char *eptr;
  int i;
  double xval,yval;
  // get past 'points'
  while (strchr(" points=\"\t",*str)!=NULL) {str++;}
  // grab x coordinate
  i=0;
  while (strchr("0123456789+-eE.",*str)!=NULL) { val[i++] = *str++; }
  val[i]=0;
  xval = strtod(val,&eptr);
  if (val==eptr) {
    printf("Error in reading x coordinate\n");
    printf("str =%s\n",str);
    return;
  }
  // grab y coordinate
  while (strchr(" \t,",*str)!=NULL) { str++; }
  i=0;
  while (strchr("0123456789+-eE.",*str)!=NULL) { val[i++] = *str++; }
  val[i]=0;
  yval = strtod(val,&eptr);
  if (val==eptr) {
    printf("Error in reading y coordinate\n");
    return;
  }

// Apply the full inverse transform with method from Reconstruct source code
  XYinverse(dim, &t[0], &t[6], &xval, &yval);
  x = xval;
  y = yval;

// Tom Bartol's initial faulty fix 
// But we have learned that that we need to apply the Inverse Transform
// Only works for dim=0 or dim=1 transform, Duh!
//  x = -t[0] + t[1]*xval + t[2]*yval + t[3]*xval*yval + t[4]*xval*xval + t[5]*yval*yval;
//  y = -t[6] + t[7]*xval + t[8]*yval + t[9]*xval*yval + t[10]*xval*xval + t[11]*yval*yval;

// Justin Kinney's original code for applying the Forward Transform
// But we have learned that that we need to apply the Inverse Transform
// Only works for identity transform, Duh!
//  x = t[0] + t[1]*xval + t[2]*yval + t[3]*xval*yval + t[4]*xval*xval + t[5]*yval*yval;
//  y = t[6] + t[7]*xval + t[8]*yval + t[9]*xval*yval + t[10]*xval*xval + t[11]*yval*yval;
}

