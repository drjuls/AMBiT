#ifndef OPERATOR_INTEGRATOR_H
#define OPERATOR_INTEGRATOR_H

#include "Universal/Lattice.h"
#include "Orbital.h"

class OPIntegrator
{
public:
    OPIntegrator(Lattice* lat) : lattice(lat) {}

    virtual double Integrate(const std::vector<double>& integrand) = 0;

    /** < a | b > = Integral (f_a * f_b + g_a * g_b) dr */
    virtual double GetInnerProduct(const SpinorFunction& a, const SpinorFunction& b) = 0;

    Lattice* GetLattice();

protected:
    Lattice* lattice;
};

inline Lattice* OPIntegrator::GetLattice()
{
    return lattice;
}

#endif
