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

/* --------------------------------------------------------------
    quick sort
----------------------------------------------------------------- */


#define MAIN
#include <stdio.h>
#include <stdlib.h>

#include "qsort.h"

static double *Val_ary;
extern int STATIC_CHQ;

void my_clear_qsort()
{
	Val_ary = NULL;
}

int int_compare(int *i, int *j)
{
    double x,y;
    x=Val_ary[*i];
    y=Val_ary[*j];
    if(x>y)return 1;
    if(x<y)return -1;
    return 0;
}

int int_compare2(const void *i, const void *j)
{
	return int_compare((int*) i, (int*) j);
}

int quick_sort(double *val_ary, int *index_ary, int index)
{
    Val_ary=val_ary;
    qsort(index_ary, index, sizeof(*index_ary), int_compare2);
    return 1;
}

