#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "Lattice.h"

class Interpolator
{
    /** Interpolator is constructed with either a lattice or a vector of points for the
        radial grid. The functions Interpolate and GetDerivative assume that the functions
        y or yvalues are defined on this grid.
     */
public:
    Interpolator(Lattice* lat): lattice(lat) {}
    Interpolator(const std::vector<double>& R_points, unsigned int order = 6);
    Interpolator(const std::vector<double>& R_points, const std::vector<double>& dR_points);
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
     */
    void GetDerivative(const std::vector<double>& y, std::vector<double>& dydr, unsigned int order);

protected:
    Lattice* lattice;
    std::vector<double> R, dR;
};

#endif
