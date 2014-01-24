#ifndef OPERATOR_INTEGRATOR_H
#define OPERATOR_INTEGRATOR_H

#include "Universal/Lattice.h"
#include "Orbital.h"
#include <boost/shared_ptr.hpp>

class OPIntegrator
{
public:
    OPIntegrator(pLattice lat): lattice(lat) {}

    virtual double Integrate(const RadialFunction& integrand) const = 0;

    /** < a | b > = Integral (f_a * f_b + g_a * g_b) dr */
    virtual double GetInnerProduct(const SpinorFunction& a, const SpinorFunction& b) const;

    pLattice GetLattice() { return lattice; }

protected:
    pLattice lattice;
};

typedef boost::shared_ptr<OPIntegrator> pOPIntegrator;

class SimpsonsIntegrator : public OPIntegrator
{
public:
    SimpsonsIntegrator(pLattice lat): OPIntegrator(lat) {}

    virtual double Integrate(const RadialFunction& integrand) const;
};

#endif
