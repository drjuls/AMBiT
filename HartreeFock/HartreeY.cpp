#include "HartreeY.h"

HartreeY::HartreeY(pOPIntegrator integration_strategy, pCoulombOperator coulomb):
    HartreeYBase(integration_strategy), coulomb(coulomb)
{}

bool HartreeY::SetParameters(int K, const SpinorFunction& c, const SpinorFunction& d)
{
    this->K = K;

    if(((K + c.L() + d.L())%2 == 0) &&
       (abs(c.TwoJ() - d.TwoJ()) <= 2 * K) &&
       (2 * K <= c.TwoJ() + d.TwoJ()))
    {
        RadialFunction density = c.GetDensity(d);
        density.resize(integrator->GetLattice()->Size());

        coulomb->GetPotential(K, density, potential);
        return true;
    }
    else
    {   potential.Clear();
        return false;
    }
}

double HartreeY::GetMatrixElement(const SpinorFunction& b, const SpinorFunction& a, bool reverse) const
{
    if(!isZero() &&
       ((K + a.L() + b.L())%2 == 0) &&
       (abs(a.TwoJ() - b.TwoJ()) <= 2 * K) &&
       (2 * K <= a.TwoJ() + b.TwoJ()))
    {
        return integrator->GetPotentialMatrixElement(a, b, potential);
    }
    else
        return 0.;
}

SpinorFunction HartreeY::ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse) const
{
    SpinorFunction ret(kappa_b);

    if(!isZero() &&
       ((K + a.L() + ret.L())%2 == 0) &&
       (abs(a.TwoJ() - ret.TwoJ()) <= 2 * K) &&
       (2 * K <= a.TwoJ() + ret.TwoJ()))
    {
        // Need to make a temporary to avoid resetting ret.kappa
        SpinorFunction temp = a * potential;
        temp.swap(ret);
    }

    return ret;
}
