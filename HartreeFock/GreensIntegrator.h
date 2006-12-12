#ifndef GREENS_INTEGRATOR_H
#define GREENS_INTEGRATOR_H

#include "Universal/Integrator.h"
#include "Universal/CoupledFunction.h"

class GreensIntegrator : public Integrator
{
    /** Solves the integral equations
            GreensInfinity(r) = (1/W) Integral  [F(r'). f_inf(r') + alpha^2. G(r'). g_inf(r')]dr'
                                     r->infinity
            GreensOrigin(r)   = (1/W) Integral [F(r'). f_0(r') + alpha^2. G(r'). g_0(r')]dr'
                                        0->r
        W is the Wronskian
            W = f_0(r)*g_inf(r) - f_inf(r)*g_0(r)
        and should be a constant. We include it in the integral just in case it isn't.

        Generally this is used to solve some coupled first order equation
        (namely the Dirac equation)
            df/dr + p(r).f(r) + q(r).g(r) = G(r)
            dg/dr + s(r).f(r) + t(r).g(r) = F(r)
        where p(r) + t(r) = 0.
        The solutions of the homogenous equation [F = G = 0] are the equations
                   s0 = (f_0(r)) should vanish at r=0
                        (g_0(r))
            sInfinity = (f_inf(r)) should vanish as r->(practical)infinity.
                        (g_inf(r))

        Also we group into a CoupledFunction rather arbitrarily:
            GreensIntegrand = (F(r))
                              (G(r))
        Note that F and G will be upper and lower components of the exchange potential.
     */
public:
    GreensIntegrator(Lattice* lat): Integrator(lat), G(NULL), s_0(NULL), s_inf(NULL) {}
    ~GreensIntegrator() {}

    /** Set a solution of the homogenous equation that is regular at r = 0. */
    inline void SetSolutionOrigin(const CoupledFunction *s0)
    {   s_0 = s0;
    }

    /** Set a solution of the homogenous equation that is regular at r -> infinity. */
    inline void SetSolutionInfinity(const CoupledFunction *sInfinity)
    {   s_inf = sInfinity;
    }

    inline void SetGreensIntegrand(const CoupledFunction *FG)
    {   G = FG;
    }

    std::vector<double> GetGreensInfinity();
    std::vector<double> GetGreensOrigin();

protected:
    // Although these are CoupledFunctions, the derivatives are ignored.
    const CoupledFunction *G;       // (F(r'), G(r'))
    const CoupledFunction *s_0;     // (f_0(r'), g_0(r'))
    const CoupledFunction *s_inf;   // (f_inf(r'), g_inf(r'))
};

#endif
