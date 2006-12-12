/* sdot.f -- translated by f2c (version of 8 February 1991  13:22:30).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/*<       real function sdot(n,sx,incx,sy,incy) >*/
doublereal sdot_(n, sx, incx, sy, incy)
integer *n;
real *sx;
integer *incx;
real *sy;
integer *incy;
{
    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer i, m;
    static real stemp;
    static integer ix, iy, mp1;


/*     forms the dot product of two vectors. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */

/*<       real sx(1),sy(1),stemp >*/
/*<       integer i,incx,incy,ix,iy,m,mp1,n >*/

/*<       stemp = 0.0e0 >*/
    /* Parameter adjustments */
    --sy;
    --sx;

    /* Function Body */
    stemp = (double)0.;
/*<       sdot = 0.0e0 >*/
    ret_val = (double)0.;
/*<       if(n.le.0)return >*/
    if (*n <= 0) {
	return ret_val;
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
/*<         stemp = stemp + sx(ix)*sy(iy) >*/
	stemp += sx[ix] * sy[iy];
/*<         ix = ix + incx >*/
	ix += *incx;
/*<         iy = iy + incy >*/
	iy += *incy;
/*<    10 continue >*/
/* L10: */
    }
/*<       sdot = stemp >*/
    ret_val = stemp;
/*<       return >*/
    return ret_val;

/*        code for both increments equal to 1 */


/*        clean-up loop */

/*<    20 m = mod(n,5) >*/
L20:
    m = *n % 5;
/*<       if( m .eq. 0 ) go to 40 >*/
    if (m == 0) {
	goto L40;
    }
/*<       do 30 i = 1,m >*/
    i__1 = m;
    for (i = 1; i <= i__1; ++i) {
/*<         stemp = stemp + sx(i)*sy(i) >*/
	stemp += sx[i] * sy[i];
/*<    30 continue >*/
/* L30: */
    }
/*<       if( n .lt. 5 ) go to 60 >*/
    if (*n < 5) {
	goto L60;
    }
/*<    40 mp1 = m + 1 >*/
L40:
    mp1 = m + 1;
/*<       do 50 i = mp1,n,5 >*/
    i__1 = *n;
    for (i = mp1; i <= i__1; i += 5) {
/*<    >*/
	stemp = stemp + sx[i] * sy[i] + sx[i + 1] * sy[i + 1] + sx[i + 2] * 
		sy[i + 2] + sx[i + 3] * sy[i + 3] + sx[i + 4] * sy[i + 4];
/*<    50 continue >*/
/* L50: */
    }
/*<    60 sdot = stemp >*/
L60:
    ret_val = stemp;
/*<       return >*/
    return ret_val;
/*<       end >*/
} /* sdot_ */

