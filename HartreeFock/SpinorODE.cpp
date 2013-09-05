#include "SpinorODE.h"

void ZeroODE::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    w[0] = w[1] = 0.;
}

void ZeroODE::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    w_f[0] = w_f[1] = 0;
    w_g[0] = w_g[1] = 0;
    w_const[0] = w_const[1] = 0;
}

void ZeroODE::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    jacobian[0][0] = 0.;
    jacobian[0][1] = 0.;
    jacobian[1][0] = 0.;
    jacobian[1][1] = 0.;

    dwdr[0] = dwdr[1] = 0.;
}

void ZeroODE::EstimateOrbitalNearOrigin(unsigned int numpoints, SpinorFunction& s) const
{
    unsigned int i;
    for(i = 0; i < numpoints; i++)
    {   s.f[i] = 0.;
        s.g[i] = 0.;
        s.dfdr[i] = 0.;
        s.dgdr[i] = 0.;
    }
}

void ZeroODE::EstimateOrbitalNearInfinity(unsigned int numpoints, Orbital& s) const
{
    unsigned int i;
    for(i = s.Size() - numpoints; i < s.Size(); i++)
    {   s.f[i] = 0.;
        s.g[i] = 0.;
        s.dfdr[i] = 0.;
        s.dgdr[i] = 0.;
    }
}
