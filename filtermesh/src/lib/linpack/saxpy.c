/* saxpy.f -- translated by f2c (version of 8 February 1991  13:22:30).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/*<       subroutine saxpy(n,sa,sx,incx,sy,incy) >*/
/* Subroutine */ int saxpy_(n, sa, sx, incx, sy, incy)
integer *n;
real *sa, *sx;
integer *incx;
real *sy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i, m, ix, iy, mp1;


/*     constant times a vector plus a vector. */
/*     uses unrolled loop for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */

/*<       real sx(1),sy(1),sa >*/
/*<       integer i,incx,incy,ix,iy,m,mp1,n >*/

/*<       if(n.le.0)return >*/
    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
/*<       if (sa .eq. 0.0) return >*/
    if (*sa == (double)0.) {
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
/*<         sy(iy) = sy(iy) + sa*sx(ix) >*/
	sy[iy] += *sa * sx[ix];
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

/*<    20 m = mod(n,4) >*/
L20:
    m = *n % 4;
/*<       if( m .eq. 0 ) go to 40 >*/
    if (m == 0) {
	goto L40;
    }
/*<       do 30 i = 1,m >*/
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
/*<         sy(i) = sy(i) + sa*sx(i) >*/
	sy[i] += *sa * sx[i];
/*<    30 continue >*/
/* L30: */
    }
/*<       if( n .lt. 4 ) return >*/
    if (*n < 4) {
	return 0;
    }
/*<    40 mp1 = m + 1 >*/
L40:
    mp1 = m + 1;
/*<       do 50 i = mp1,n,4 >*/
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 4) {
/*<         sy(i) = sy(i) + sa*sx(i) >*/
	sy[i] += *sa * sx[i];
/*<         sy(i + 1) = sy(i + 1) + sa*sx(i + 1) >*/
	sy[i + 1] += *sa * sx[i + 1];
/*<         sy(i + 2) = sy(i + 2) + sa*sx(i + 2) >*/
	sy[i + 2] += *sa * sx[i + 2];
/*<         sy(i + 3) = sy(i + 3) + sa*sx(i + 3) >*/
	sy[i + 3] += *sa * sx[i + 3];
/*<    50 continue >*/
/* L50: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* saxpy_ */

