#ifndef BSPLINE_GRID_H
#define BSPLINE_GRID_H

#include "Universal/Lattice.h"

class BSplineGrid : public Lattice
{
    /** The spline grid is made of knots t[i] interspersed with gaussian points
     *      r[i] = t[i] + (t[i+1] - t[i])*xgauss[m]
     *  where the knots are defined on a logarithmic scale
     *      t[i] = beta * (exp(h*(i-k+1)) - 1), i = k-1..n
     *  subject to constraints
     *      t[0] = t[1] = ... = t[k-1] = 0
     *      t[k] = dr0
     *      t[n] = t[n+1] = ... = t[n+k-1] = rmax
     */
public:
    BSplineGrid(unsigned int n, unsigned int k, double dr0, double rmax);
    virtual ~BSplineGrid(void);

    const double* GetSplineKnots() { return t; }
    unsigned int GetN() { return n; }
    unsigned int GetK() { return k; }

protected:
    /** Calculate the value that r[i] should be. */
    double lattice_to_real(unsigned int i) const;

    /** Calculate the lattice spacing at a point. */
    double calculate_dr(double r_point) const;

    /** Resizes the lattice such that NumPoints > min_size. */
    void ReSize(unsigned int min_size);

    /** Create a grid (knot sequence) for the splines.
     *  PRE: t[n+k]
     *  If there are to be n splines in total then there are n+k knots.
     *      t[0] = t[1] =...= t[k-1] = 0
     *      t[n] = t[n+1] =...= t[n+k-1] = rmax
     *  The sequence itself is logarithmic, constrained by
     *      t[k] - t[k-1] = dr0
     */
    void CreateSplineKnots(double dr0, double rmax);

protected:
    unsigned int n;
    unsigned int k;

    // B-Spline Knots
    double* t;

    // Gaussian points and weights for intervals between knots
    double* xgauss;
    double* wgauss;
};

#endif
