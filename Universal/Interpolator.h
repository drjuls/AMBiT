#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "Lattice.h"

class Interpolator
{
public:
    Interpolator(Lattice* lat): lattice(lat) {}
    ~Interpolator() {}

    /** Interpolate yvalues over the lattice lat, using points preferably
     *  centered around xvalue.
     *  Aitken's method is used to construct a Lagrange interpolating polynomial
     *  using "order" number of points.
     *  Derivative is dy/dr at xvalue.
     */
    void Interpolate(const std::vector<double>& yvalues, double xvalue,
                     double& yvalue, double& derivative, unsigned int order);

    /** Get derivative of y on lattice from Lagrange interpolation.
     *      dy[i] = (dy/dr)[i] * dr[i]
     */
    void GetDerivative(const std::vector<double>& y, std::vector<double>& dy, unsigned int order);

protected:
    Lattice* lattice;
};

#endif