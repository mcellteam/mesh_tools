/* snrm2.f -- translated by f2c (version of 8 February 1991  13:22:30).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/

#include "f2c.h"

/*<       real function snrm2 ( n, sx, incx) >*/
doublereal snrm2_(n, sx, incx)
integer *n;
real *sx;
integer *incx;
{
    /* Initialized data */

    static real zero = (double)0.;
    static real one = (double)1.;
    static real cutlo = (double)4.441e-16;
    static real cuthi = (double)1.304e19;

    /* Format strings */
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_110[] = "";

    /* System generated locals */
    integer i__1;
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static real xmax;
    static integer next, i, j, ix;
    static real hitest, sum;

    /* Assigned format variables */
    char *next_fmt;

/*<       integer i, incx, ix, j, n, next >*/
/*<       real   sx(1),  cutlo, cuthi, hitest, sum, xmax, zero, one >*/
/*<       data   zero, one /0.0e0, 1.0e0/ >*/
    /* Parameter adjustments */
    --sx;

    /* Function Body */

/*     euclidean norm of the n-vector stored in sx() with storage */
/*     increment incx . */
/*     if    n .le. 0 return with result = 0. */
/*     if n .ge. 1 then incx must be .ge. 1 */

/*           c.l.lawson, 1978 jan 08 */
/*     modified to correct failure to update ix, 1/25/92. */
/*     modified 3/93 to return if incx .le. 0. */

/*     four phase method     using two built-in constants that are */
/*     hopefully applicable to all machines. */
/*         cutlo = maximum of  sqrt(u/eps)  over all known machines. */
/*         cuthi = minimum of  sqrt(v)      over all known machines. */
/*     where */
/*         eps = smallest no. such that eps + 1. .gt. 1. */
/*         u   = smallest positive no.   (underflow limit) */
/*         v   = largest  no.            (overflow  limit) */

/*     brief outline of algorithm.. */

/*     phase 1    scans zero components. */
/*     move to phase 2 when a component is nonzero and .le. cutlo */
/*     move to phase 3 when a component is .gt. cutlo */
/*     move to phase 4 when a component is .ge. cuthi/m */
/*     where m = n for x() real and m = 2*n for complex. */

/*     values for cutlo and cuthi.. */
/*     from the environmental parameters listed in the imsl converter */
/*     document the limiting values are as follows.. */
/*     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are 
*/
/*                   univac and dec at 2**(-103) */
/*                   thus cutlo = 2**(-51) = 4.44089e-16 */
/*     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec. */
/*                   thus cuthi = 2**(63.5) = 1.30438e19 */
/*     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec. */
/*                   thus cutlo = 2**(-33.5) = 8.23181d-11 */
/*     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19 */
/*     data cutlo, cuthi / 8.232d-11,  1.304d19 / */
/*     data cutlo, cuthi / 4.441e-16,  1.304e19 / */
/*<       data cutlo, cuthi / 4.441e-16,  1.304e19 / >*/

/*<       if(n .gt. 0 .and. incx.gt.0) go to 10 >*/
    if (*n > 0 && *incx > 0) {
	goto L10;
    }
/*<          snrm2  = zero >*/
    ret_val = zero;
/*<          go to 300 >*/
    goto L300;

/*<    10 assign 30 to next >*/
L10:
    next = 0;
    next_fmt = fmt_30;
/*<       sum = zero >*/
    sum = zero;
/*<       i = 1 >*/
    i = 1;
/*<       ix = 1 >*/
    ix = 1;
/*                                                 begin main loop */
/*<    20    go to next,(30, 50, 70, 110) >*/
L20:
    switch ((int)next) {
	case 0: goto L30;
	case 1: goto L50;
	case 2: goto L70;
	case 3: goto L110;
    }
/*<    30 if( abs(sx(i)) .gt. cutlo) go to 85 >*/
L30:
    if ((r__1 = sx[i], dabs(r__1)) > cutlo) {
	goto L85;
    }
/*<       assign 50 to next >*/
    next = 1;
    next_fmt = fmt_50;
/*<       xmax = zero >*/
    xmax = zero;

/*                        phase 1.  sum is zero */

/*<    50 if( sx(i) .eq. zero) go to 200 >*/
L50:
    if (sx[i] == zero) {
	goto L200;
    }
/*<       if( abs(sx(i)) .gt. cutlo) go to 85 >*/
    if ((r__1 = sx[i], dabs(r__1)) > cutlo) {
	goto L85;
    }

/*                                prepare for phase 2. */
/*<       assign 70 to next >*/
    next = 2;
    next_fmt = fmt_70;
/*<       go to 105 >*/
    goto L105;

/*                                prepare for phase 4. */

/*<   100 continue >*/
L100:
/*<       ix = j >*/
    ix = j;
/*<       assign 110 to next >*/
    next = 3;
    next_fmt = fmt_110;
/*<       sum = (sum / sx(i)) / sx(i) >*/
    sum = sum / sx[i] / sx[i];
/*<   105 xmax = abs(sx(i)) >*/
L105:
    xmax = (r__1 = sx[i], dabs(r__1));
/*<       go to 115 >*/
    goto L115;

/*                   phase 2.  sum is small. */
/*                             scale to avoid destructive underflow. */

/*<    70 if( abs(sx(i)) .gt. cutlo ) go to 75 >*/
L70:
    if ((r__1 = sx[i], dabs(r__1)) > cutlo) {
	goto L75;
    }

/*                     common code for phases 2 and 4. */
/*                     in phase 4 sum is large.  scale to avoid overflow. 
*/

/*<   110 if( abs(sx(i)) .le. xmax ) go to 115 >*/
L110:
    if ((r__1 = sx[i], dabs(r__1)) <= xmax) {
	goto L115;
    }
/*<          sum = one + sum * (xmax / sx(i))**2 >*/
/* Computing 2nd power */
    r__1 = xmax / sx[i];
    sum = one + sum * (r__1 * r__1);
/*<          xmax = abs(sx(i)) >*/
    xmax = (r__1 = sx[i], dabs(r__1));
/*<          go to 200 >*/
    goto L200;

/*<   115 sum = sum + (sx(i)/xmax)**2 >*/
L115:
/* Computing 2nd power */
    r__1 = sx[i] / xmax;
    sum += r__1 * r__1;
/*<       go to 200 >*/
    goto L200;


/*                  prepare for phase 3. */

/*<    75 sum = (sum * xmax) * xmax >*/
L75:
    sum = sum * xmax * xmax;


/*     for real or d.p. set hitest = cuthi/n */
/*     for complex      set hitest = cuthi/(2*n) */

/*<    85 hitest = cuthi/double( n ) >*/
L85:
    hitest = cuthi / (real) (*n);

/*                   phase 3.  sum is mid-range.  no scaling. */

/*<       do 95 j = ix, n >*/
    i__1 = *n;
    for (j = ix; j <= i__1; ++j) {
/*<          if(abs(sx(i)) .ge. hitest) go to 100 >*/
	if ((r__1 = sx[i], dabs(r__1)) >= hitest) {
	    goto L100;
	}
/*<          sum = sum + sx(i)**2 >*/
/* Computing 2nd power */
	r__1 = sx[i];
	sum += r__1 * r__1;
/*<          i = i + incx >*/
	i += *incx;
/*<    95 continue >*/
/* L95: */
    }
/*<       snrm2 = sqrt( sum ) >*/
    ret_val = sqrt(sum);
/*<       go to 300 >*/
    goto L300;

/*<   200 continue >*/
L200:
/*<       ix = ix + 1 >*/
    ++ix;
/*<       i = i + incx >*/
    i += *incx;
/*<       if( ix .le. n ) go to 20 >*/
    if (ix <= *n) {
	goto L20;
    }

/*              end of main loop. */

/*              compute square root and adjust for scaling. */

/*<       snrm2 = xmax * sqrt(sum) >*/
    ret_val = xmax * sqrt(sum);
/*<   300 continue >*/
L300:
/*<       return >*/
    return ret_val;
/*<       end >*/
} /* snrm2_ */

