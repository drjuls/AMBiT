#include "NonRelativisticSMSOperator.h"
#include "Universal/Interpolator.h"

SpinorFunction NonRelativisticSMSOperator::ApplyTo(const SpinorFunction& a) const
{
    return OneBodyOperatorDecorator::ApplyTo(a);
}

SpinorFunction NonRelativisticSMSOperator::ApplyTo(const SpinorFunction& a, int kappa_b) const
{
    int l_b = kappa_b;
    if(kappa_b < 0)
        l_b = - kappa_b - 1;

    SpinorFunction wrapped_result = OneBodyOperatorDecorator::ApplyTo(a, kappa_b);
    double coeff_fa = 0.;

    // Check angular momenta delta functions
    if(l_b == a.L()+1)
        coeff_fa = -double(l_b);
    else if(l_b == a.L()-1)
        coeff_fa = double(a.L());
    else
        return wrapped_result;

    SpinorFunction p(kappa_b, a.Size());

    pLattice lattice = integrator->GetLattice();
    const double* R = lattice->R();

    unsigned int i = 0;
    for(i = 0; i < a.Size(); i++)
    {   p.f[i] = a.dfdr[i] + coeff_fa/R[i] * a.f[i];
    }

    // Get derivative of p->f
    Interpolator I(lattice);
    I.GetDerivative(a.dfdr, p.dfdr, 6);  // First term
    // Second term
    for(i = 0; i < p.Size(); i++)
        p.dfdr[i] += coeff_fa/R[i] * (a.dfdr[i] - a.f[i]/R[i]);

    p += wrapped_result;
    return p;
}