/* sscal.f -- translated by f2c (version of 8 February 1991  13:22:30).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/*<       subroutine sscal(n,sa,sx,incx) >*/
/* Subroutine */ int sscal_(n, sa, sx, incx)
integer *n;
real *sa, *sx;
integer *incx;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, m, nincx, mp1;


/*     scales a vector by a constant. */
/*     uses unrolled loops for increment equal to 1. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */

/*<       real sa,sx(1) >*/
/*<       integer i,incx,m,mp1,n,nincx >*/

/*<       if( n.le.0 .or. incx.le.0 )return >*/
    /* Parameter adjustments */
    --sx;

    /* Function Body */
    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
/*<       if(incx.eq.1)go to 20 >*/
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

/*<       nincx = n*incx >*/
    nincx = *n * *incx;
/*<       do 10 i = 1,nincx,incx >*/
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2) {
/*<         sx(i) = sa*sx(i) >*/
	sx[i] = *sa * sx[i];
/*<    10 continue >*/
/* L10: */
    }
/*<       return >*/
    return 0;

/*        code for increment equal to 1 */


/*        clean-up loop */

/*<    20 m = mod(n,5) >*/
L20:
    m = *n % 5;
/*<       if( m .eq. 0 ) go to 40 >*/
    if (m == 0) {
	goto L40;
    }
/*<       do 30 i = 1,m >*/
    i__2 = m;
    for (i = 1; i <= i__2; ++i) {
/*<         sx(i) = sa*sx(i) >*/
	sx[i] = *sa * sx[i];
/*<    30 continue >*/
/* L30: */
    }
/*<       if( n .lt. 5 ) return >*/
    if (*n < 5) {
	return 0;
    }
/*<    40 mp1 = m + 1 >*/
L40:
    mp1 = m + 1;
/*<       do 50 i = mp1,n,5 >*/
    i__2 = *n;
    for (i = mp1; i <= i__2; i += 5) {
/*<         sx(i) = sa*sx(i) >*/
	sx[i] = *sa * sx[i];
/*<         sx(i + 1) = sa*sx(i + 1) >*/
	sx[i + 1] = *sa * sx[i + 1];
/*<         sx(i + 2) = sa*sx(i + 2) >*/
	sx[i + 2] = *sa * sx[i + 2];
/*<         sx(i + 3) = sa*sx(i + 3) >*/
	sx[i + 3] = *sa * sx[i + 3];
/*<         sx(i + 4) = sa*sx(i + 4) >*/
	sx[i + 4] = *sa * sx[i + 4];
/*<    50 continue >*/
/* L50: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* sscal_ */

