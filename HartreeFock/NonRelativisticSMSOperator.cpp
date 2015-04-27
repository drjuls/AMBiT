#include "NonRelativisticSMSOperator.h"
#include "Universal/Interpolator.h"

bool NonRelativisticSMSOperator::SetLocalParameters(int new_K, pSpinorFunctionConst new_c, pSpinorFunctionConst new_d)
{
    K = new_K;
    c = new_c;
    d = new_d;

    if(K == 1)
    {
        SpinorFunction td = ApplyOperator(*new_d, new_c->Kappa());
        p_cd = integrator->GetInnerProduct(td, *new_c);
    }
    else
    {   p_cd = 0.0;
    }

    return (p_cd != 0.0);
}

bool NonRelativisticSMSOperator::isZero() const
{
    bool ret = HartreeYDecorator::isZero();
    return ret && (p_cd * lambda == 0);
}

/** < b | t | a > for an operator t. */
double NonRelativisticSMSOperator::GetMatrixElement(const Orbital& b, const Orbital& a, bool reverse) const
{
    if(!integrator)
        throw "NonRelativisticSMSOperator::GetMatrixElement(): no integrator found.";

    SpinorFunction ta = ApplyTo(a, b.Kappa(), reverse);
    return integrator->GetInnerProduct(ta, b);
}

SpinorFunction NonRelativisticSMSOperator::ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse) const
{
    SpinorFunction ret = ApplyOperator(a, kappa_b) * lambda * p_cd;
    if(reverse)
        ret = ret * (-1.0);

    ret += HartreeYDecorator::ApplyTo(a, kappa_b, reverse);

    return ret;
}

SpinorFunction NonRelativisticSMSOperator::ApplyOperator(const SpinorFunction& a, int kappa_b) const
{
    if(!integrator)
        throw "NonRelativisticSMSOperator::GetMatrixElement(): no integrator found.";

    int l_b = kappa_b;
    if(kappa_b < 0)
        l_b = - kappa_b - 1;

    double coeff_fa = 0.;

    // Check angular momenta delta functions
    if(l_b == a.L()+1)
        coeff_fa = -double(l_b);
    else if(l_b == a.L()-1)
        coeff_fa = double(a.L());
    else
        return SpinorFunction(kappa_b);

    SpinorFunction p(kappa_b, a.size());

    pLattice lattice = integrator->GetLattice();
    const double* R = lattice->R();

    unsigned int i = 0;
    for(i = 0; i < a.size(); i++)
    {   p.f[i] = a.dfdr[i] + coeff_fa/R[i] * a.f[i];
    }

    // Get derivative of p->f
    Interpolator I(lattice);
    I.GetDerivative(a.dfdr, p.dfdr, 6);  // First term
    // Second term
    for(i = 0; i < p.size(); i++)
        p.dfdr[i] += coeff_fa/R[i] * (a.dfdr[i] - a.f[i]/R[i]);

    return p;
}
