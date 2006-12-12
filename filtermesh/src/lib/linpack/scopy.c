/* scopy.f -- translated by f2c (version of 8 February 1991  13:22:30).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/*<       subroutine scopy(n,sx,incx,sy,incy) >*/
/* Subroutine */ int scopy_(n, sx, incx, sy, incy)
integer *n;
real *sx;
integer *incx;
real *sy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     copies a vector, x, to a vector, y. */
/*     uses unrolled loops for increments equal to 1. */
/*     jack dongarra, linpack, 3/11/78. */

/*<       real sx(1),sy(1) >*/
/*<       integer i,incx,incy,ix,iy,m,mp1,n >*/

/*<       if(n.le.0)return >*/
    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
/*<       if(incx.eq.1.and.incy.eq.1)go to 20 >*/
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

/*<       ix = 1 >*/
    ix = 1;
/*<       iy = 1 >*/
    iy = 1;
/*<       if(incx.lt.0)ix = (-n+1)*incx + 1 >*/
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
/*<       if(incy.lt.0)iy = (-n+1)*incy + 1 >*/
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
/*<       do 10 i = 1,n >*/
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/*<         sy(iy) = sx(ix) >*/
	sy[iy] = sx[ix];
/*<         ix = ix + incx >*/
	ix += *incx;
/*<         iy = iy + incy >*/
	iy += *incy;
/*<    10 continue >*/
/* L10: */
    }
/*<       return >*/
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

/*<    20 m = mod(n,7) >*/
L20:
    m = *n % 7;
/*<       if( m .eq. 0 ) go to 40 >*/
    if (m == 0) {
	goto L40;
    }
/*<       do 30 i = 1,m >*/
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
/*<         sy(i) = sx(i) >*/
	sy[i] = sx[i];
/*<    30 continue >*/
/* L30: */
    }
/*<       if( n .lt. 7 ) return >*/
    if (*n < 7) {
	return 0;
    }
/*<    40 mp1 = m + 1 >*/
L40:
    mp1 = m + 1;
/*<       do 50 i = mp1,n,7 >*/
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 7) {
/*<         sy(i) = sx(i) >*/
	sy[i] = sx[i];
/*<         sy(i + 1) = sx(i + 1) >*/
	sy[i + 1] = sx[i + 1];
/*<         sy(i + 2) = sx(i + 2) >*/
	sy[i + 2] = sx[i + 2];
/*<         sy(i + 3) = sx(i + 3) >*/
	sy[i + 3] = sx[i + 3];
/*<         sy(i + 4) = sx(i + 4) >*/
	sy[i + 4] = sx[i + 4];
/*<         sy(i + 5) = sx(i + 5) >*/
	sy[i + 5] = sx[i + 5];
/*<         sy(i + 6) = sx(i + 6) >*/
	sy[i + 6] = sx[i + 6];
/*<    50 continue >*/
/* L50: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* scopy_ */

