/* Splines.f -- translated by f2c (version 19990311).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int bsplvb_(t, jhigh, index, x, left, biatx)
doublereal *t;
integer *jhigh, *index;
doublereal *x;
integer *left;
doublereal *biatx;
{
    /* Initialized data */

    static integer j = 1;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal term;
    static integer i__;
    static doublereal saved, deltal[20], deltar[20];
    static integer jp1;

    /* Parameter adjustments */
    --t;
    --biatx;

    /* Function Body */
    switch ((int)*index) {
	case 1:  goto L10;
	case 2:  goto L20;
    }
L10:
    j = 1;
    biatx[1] = (float)1.;
    if (j >= *jhigh) {
	goto L99;
    }
L20:
    jp1 = j + 1;
    deltar[j - 1] = t[*left + j] - *x;
    deltal[j - 1] = *x - t[*left + 1 - j];
    saved = (float)0.;
    i__1 = j;
    for (i__ = 1; i__ <= i__1; ++i__) {
	term = biatx[i__] / (deltar[i__ - 1] + deltal[jp1 - i__ - 1]);
	biatx[i__] = saved + deltar[i__ - 1] * term;
	saved = deltal[jp1 - i__ - 1] * term;
/* L26: */
    }
    biatx[jp1] = saved;
    j = jp1;
    if (j < *jhigh) {
	goto L20;
    }
L99:
    return 0;
} /* bsplvb_ */

/* Subroutine */ int bsplvd(t, k, x, left, dbiatx, nderiv)
doublereal *t;
integer *k;
doublereal *x;
integer *left;
doublereal *dbiatx;
integer *nderiv;
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer jlow, kp1mm;
    static doublereal a[225]	/* was [15][15] */;
    static integer i__, j, m, mhigh, jp1mid;
    static doublereal fkp1mm;
    static integer il;
    static doublereal factor;
    static integer ideriv;
    extern /* Subroutine */ int bsplvb_();
    static integer ldummy, kp1;
    static doublereal sum;

/* ******  changes from de Boor */
/*  also a missing in arg list */
/* **************** */
    /* Parameter adjustments */
    dbiatx -= 16;
    --t;

    /* Function Body */
    for (i__ = 1; i__ <= 15; ++i__) {
	for (j = 1; j <= 15; ++j) {
	    a[i__ + j * 15 - 16] = (float)0.;
	}
    }
/* Computing MAX */
    i__1 = min(*nderiv,*k);
    mhigh = max(i__1,1);
    kp1 = *k + 1;
    i__1 = kp1 - mhigh;
    bsplvb_(&t[1], &i__1, &c__1, x, left, &dbiatx[16]);
    if (mhigh == 1) {
	goto L99;
    }
    ideriv = mhigh;
    i__1 = mhigh;
    for (m = 2; m <= i__1; ++m) {
	jp1mid = 1;
	i__2 = *k;
	for (j = ideriv; j <= i__2; ++j) {
	    dbiatx[j + ideriv * 15] = dbiatx[jp1mid + 15];
	    ++jp1mid;
/* L11: */
	}
	--ideriv;
	i__2 = kp1 - ideriv;
	bsplvb_(&t[1], &i__2, &c__2, x, left, &dbiatx[16]);
/* L15: */
    }
    jlow = 1;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *k;
	for (j = jlow; j <= i__2; ++j) {
	    a[j + i__ * 15 - 16] = (float)0.;
/* L19: */
	}
	jlow = i__;
	a[i__ + i__ * 15 - 16] = (float)1.;
/* L20: */
    }
    i__1 = mhigh;
    for (m = 2; m <= i__1; ++m) {
	kp1mm = kp1 - m;
	fkp1mm = (doublereal) kp1mm;
	il = *left;
	i__ = *k;
	i__2 = kp1mm;
	for (ldummy = 1; ldummy <= i__2; ++ldummy) {
	    factor = fkp1mm / (t[il + kp1mm] - t[il]);
	    i__3 = i__;
	    for (j = 1; j <= i__3; ++j) {
		a[i__ + j * 15 - 16] = (a[i__ + j * 15 - 16] - a[i__ - 1 + j *
			 15 - 16]) * factor;
/* L24: */
	    }
	    --il;
	    --i__;
/* L25: */
	}
	i__2 = *k;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum = (float)0.;
	    jlow = max(i__,m);
	    i__3 = *k;
	    for (j = jlow; j <= i__3; ++j) {
		sum += a[j + i__ * 15 - 16] * dbiatx[j + m * 15];
/* L35: */
	    }
	    dbiatx[i__ + m * 15] = sum;
/* L36: */
	}
/* L40: */
    }
L99:
    return 0;
} /* bsplvd_ */

doublereal bvalue(t, bcoef, n, k, x, jderiv)
doublereal *t, *bcoef;
integer *n, *k;
doublereal *x;
integer *jderiv;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    static doublereal fkmj;
    static integer i__, j, mflag, jcmin, jcmax;
    static doublereal aj[20];
    static integer jc;
    static doublereal dl[20];
    static integer jj;
    static doublereal dr[20];
    extern /* Subroutine */ int interv_();
    static integer km1, imk, kmj, nmi, ilo;

/*     Function returns the value of the (jderiv)th derivative of */
/*       f(x) = Sum[bcoef(i) * BSpline(i, x), i=1..n] */
/*     at the point x */
    /* Parameter adjustments */
    --t;
    --bcoef;

    /* Function Body */
    ret_val = (float)0.;
    if (*jderiv >= *k) {
	goto L99;
    }
    i__1 = *n + *k;
    interv_(&t[1], &i__1, x, &i__, &mflag);
    if (mflag != 0) {
	goto L99;
    }
    km1 = *k - 1;
    if (km1 > 0) {
	goto L1;
    }
    ret_val = bcoef[i__];
    goto L99;
L1:
    jcmin = 1;
    imk = i__ - *k;
    if (imk >= 0) {
	goto L8;
    }
    jcmin = 1 - imk;
    i__1 = i__;
    for (j = 1; j <= i__1; ++j) {
	dl[j - 1] = *x - t[i__ + 1 - j];
/* L5: */
    }
    i__1 = km1;
    for (j = i__; j <= i__1; ++j) {
	aj[*k - j - 1] = (float)0.;
	dl[j - 1] = dl[i__ - 1];
/* L6: */
    }
    goto L10;
L8:
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	dl[j - 1] = *x - t[i__ + 1 - j];
/* L9: */
    }
L10:
    jcmax = *k;
    nmi = *n - i__;
    if (nmi >= 0) {
	goto L18;
    }
    jcmax = *k + nmi;
    i__1 = jcmax;
    for (j = 1; j <= i__1; ++j) {
	dr[j - 1] = t[i__ + j] - *x;
/* L15: */
    }
    i__1 = km1;
    for (j = jcmax; j <= i__1; ++j) {
	aj[j] = (float)0.;
	dr[j - 1] = dr[jcmax - 1];
/* L16: */
    }
    goto L20;
L18:
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	dr[j - 1] = t[i__ + j] - *x;
/* L19: */
    }
L20:
    i__1 = jcmax;
    for (jc = jcmin; jc <= i__1; ++jc) {
	aj[jc - 1] = bcoef[imk + jc];
/* L21: */
    }
    if (*jderiv == 0) {
	goto L30;
    }
    i__1 = *jderiv;
    for (j = 1; j <= i__1; ++j) {
	kmj = *k - j;
	fkmj = (doublereal) kmj;
	ilo = kmj;
	i__2 = kmj;
	for (jj = 1; jj <= i__2; ++jj) {
	    aj[jj - 1] = (aj[jj] - aj[jj - 1]) / (dl[ilo - 1] + dr[jj - 1]) * 
		    fkmj;
	    --ilo;
/* L22: */
	}
/* L23: */
    }
L30:
    if (*jderiv == km1) {
	goto L39;
    }
    i__1 = km1;
    for (j = *jderiv + 1; j <= i__1; ++j) {
	kmj = *k - j;
	ilo = kmj;
	i__2 = kmj;
	for (jj = 1; jj <= i__2; ++jj) {
	    aj[jj - 1] = (aj[jj] * dl[ilo - 1] + aj[jj - 1] * dr[jj - 1]) / (
		    dl[ilo - 1] + dr[jj - 1]);
	    --ilo;
/* L32: */
	}
/* L33: */
    }
L39:
    ret_val = aj[0];
L99:
    return ret_val;
} /* bvalue_ */

/* Subroutine */ int gauss(n, x, w)
integer *n;
doublereal *x, *w;
{
    /* Initialized data */

    static doublereal eps = 1e-15;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double acos(), cos();

    /* Local variables */
    static integer i__, j, m;
    static doublereal z__, p1, p2, p3, z1, pi, pp;

/* ***************************************************************** */

/*  Gaussian coordinates and weights for the interval [0..1] */
/*        adapted from "setgau" in Numerical Recipes */
/* ****************************************************************** */
    /* Parameter adjustments */
    --w;
    --x;

    /* Function Body */
    pi = acos(-1.);
    m = (*n + 1) / 2;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ = cos(pi * (i__ - .25) / (*n + .5));
L100:
	p1 = 1.;
	p2 = 0.;
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    p3 = p2;
	    p2 = p1;
	    p1 = ((j * 2. - 1.) * z__ * p2 - (j - 1.) * p3) / j;
/* L110: */
	}
	pp = *n * (z__ * p1 - p2) / (z__ * z__ - 1.);
	z1 = z__;
	z__ = z1 - p1 / pp;
	if ((d__1 = z__ - z1, abs(d__1)) > eps) {
	    goto L100;
	}
	x[i__] = (1. - z__) * .5;
	x[*n + 1 - i__] = (z__ + 1.) * .5;
	w[i__] = 1. / ((1. - z__ * z__) * pp * pp);
	w[*n + 1 - i__] = w[i__];
/* L120: */
    }
    return 0;
} /* gauss_ */

/* Subroutine */ int interv_(xt, lxt, x, left, mflag)
doublereal *xt;
integer *lxt;
doublereal *x;
integer *left, *mflag;
{
    /* Initialized data */

    static integer ilo = 1;

    static integer istep, middle, ihi;

    /* Parameter adjustments */
    --xt;

    /* Function Body */
    ihi = ilo + 1;
    if (ihi < *lxt) {
	goto L20;
    }
    if (*x >= xt[*lxt]) {
	goto L110;
    }
    if (*lxt <= 1) {
	goto L90;
    }
    ilo = *lxt - 1;
    ihi = *lxt;
L20:
    if (*x >= xt[ihi]) {
	goto L40;
    }
    if (*x >= xt[ilo]) {
	goto L100;
    }
    istep = 1;
L31:
    ihi = ilo;
    ilo = ihi - istep;
    if (ilo <= 1) {
	goto L35;
    }
    if (*x >= xt[ilo]) {
	goto L50;
    }
    istep <<= 1;
    goto L31;
L35:
    ilo = 1;
    if (*x < xt[1]) {
	goto L90;
    }
    goto L50;
L40:
    istep = 1;
L41:
    ilo = ihi;
    ihi = ilo + istep;
    if (ihi >= *lxt) {
	goto L45;
    }
    if (*x < xt[ihi]) {
	goto L50;
    }
    istep <<= 1;
    goto L41;
L45:
    if (*x >= xt[*lxt]) {
	goto L110;
    }
    ihi = *lxt;
L50:
    middle = (ilo + ihi) / 2;
    if (middle == ilo) {
	goto L100;
    }
    if (*x < xt[middle]) {
	goto L53;
    }
    ilo = middle;
    goto L50;
L53:
    ihi = middle;
    goto L50;
L90:
    *mflag = -1;
    *left = 1;
    return 0;
L100:
    *mflag = 0;
    *left = ilo;
    return 0;
L110:
    *mflag = 1;
    *left = *lxt;
    return 0;
} /* interv_ */

