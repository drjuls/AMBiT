#ifndef SPLINE_H
#define SPLINE_H

extern "C"{
/** Get value of the B-Splines of order k and derivatives at point
 *      x = [t[left], t[left+1]).
 *  where t is the array of knot points.
 *  POST: dbiatx[m, n] = (n-1)th order derivate of B_(left-k+m) at point x.
 *          where m = 1..k  and
 *          n = 1..nderiv
 */
void bsplvd_(const double* t, int* k, double* x, int* left, double* dbiatx, int* nderiv);

/** Function returns the value of the (jderiv)th derivative of
 *      f(x) = Sum[bcoef(i) * BSpline(i, x), i=1..n]
 *  at the point x.
 */
double bvalue_(const double* t, const double* bcoef, int* n, int* k, double* x, int* jderiv);

const int MaximumK = 15;    //Parameter KX in spline.f

/** Gaussian coordinates (x[n]) and weights (w[n]) for the interval [0..1].
 *  Adapted from "setgau" in Numerical Recipes.
 */
void gauss_(int* n, double* x, double* w);
}

#endif
