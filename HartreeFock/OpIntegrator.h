#ifndef OPERATOR_INTEGRATOR_H
#define OPERATOR_INTEGRATOR_H

#include "Universal/Lattice.h"
#include "Orbital.h"
#include <boost/shared_ptr.hpp>

class OPIntegrator
{
public:
    OPIntegrator(Lattice* lat): lattice(lat) {}

    virtual double Integrate(const RadialFunction& integrand) const = 0;

    /** < a | b > = Integral (f_a * f_b + g_a * g_b) dr */
    virtual double GetInnerProduct(const SpinorFunction& a, const SpinorFunction& b) const;

    Lattice* GetLattice() { return lattice; }

protected:
    Lattice* lattice;
};

typedef boost::shared_ptr<OPIntegrator> pOPIntegrator;

class SimpsonsIntegrator : public OPIntegrator
{
public:
    SimpsonsIntegrator(Lattice* lat): OPIntegrator(lat) {}

    virtual double Integrate(const RadialFunction& integrand) const;
};

#endif
