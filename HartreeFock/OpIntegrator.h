#ifndef OPERATOR_INTEGRATOR_H
#define OPERATOR_INTEGRATOR_H

#include "Universal/Lattice.h"
#include "Universal/SpinorFunction.h"
#include <boost/shared_ptr.hpp>

class OPIntegrator
{
public:
    OPIntegrator(pLattice lat): lattice(lat) {}

    virtual double Integrate(const RadialFunction& integrand) const = 0;

    /** < a | b > = Integral (f_a * f_b + g_a * g_b) dr */
    virtual double GetInnerProduct(const SpinorFunction& a, const SpinorFunction& b) const;

    /** < a | a > */
    virtual double GetNorm(const SpinorFunction& a) const;

    /** < a | V | b > = Integral (f_a * f_b + g_a * g_b) * V(r) dr */
    virtual double GetPotentialMatrixElement(const SpinorFunction& a, const SpinorFunction& b, const RadialFunction& V) const;

    /** < a | V | a > */
    virtual double GetPotentialMatrixElement(const SpinorFunction& a, const RadialFunction& V) const;

    pLattice GetLattice() { return lattice; }

protected:
    pLattice lattice;
};

typedef boost::shared_ptr<OPIntegrator> pOPIntegrator;

/** SimpsonsIntegrator overrides many functions in OPIntegrator for speed. */
class SimpsonsIntegrator : public OPIntegrator
{
public:
    SimpsonsIntegrator(pLattice lat): OPIntegrator(lat) {}

    virtual double Integrate(const RadialFunction& integrand) const override;

    /** < a | b > = Integral (f_a * f_b + g_a * g_b) dr */
    virtual double GetInnerProduct(const SpinorFunction& a, const SpinorFunction& b) const override;

    /** < a | a > */
    virtual double GetNorm(const SpinorFunction& a) const override;

    /** < a | V | b > = Integral (f_a * f_b + g_a * g_b) * V(r) dr */
    virtual double GetPotentialMatrixElement(const SpinorFunction& a, const SpinorFunction& b, const RadialFunction& V) const override;

    /** < a | V | a > */
    virtual double GetPotentialMatrixElement(const SpinorFunction& a, const RadialFunction& V) const override;

protected:
    template<typename LambdaIntegrand>
    double Integrate(int size, LambdaIntegrand&& integrand) const
    {
        double total = 0.;
        const double* dR = lattice->dR();

        int i;
        for(i = 1; i < size-1; i+=2)
        {
            total += 4. * integrand(i) * dR[i]
                    + 2. * integrand(i+1) * dR[i+1];
        }
        total = total/3.;
        total += integrand(0) * dR[0];

        while(i < size)
        {   total += integrand(i) * dR[i];
            i++;
        }
        
        return total;
    }
};

#endif
