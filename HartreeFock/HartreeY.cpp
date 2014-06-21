#include "HartreeY.h"

HartreeY::HartreeY(pOPIntegrator integration_strategy, pCoulombOperator coulomb):
    HartreeYBase(integration_strategy), LatticeObserver(integration_strategy->GetLattice()), coulomb(coulomb)
{
    two_body_reverse_symmetry_exists = true;
}

bool HartreeY::SetParameters(int K, const SpinorFunction& c, const SpinorFunction& d)
{
    this->K = K;

    if(((K + c.L() + d.L())%2 == 0) &&
       (abs(c.TwoJ() - d.TwoJ()) <= 2 * K) &&
       (2 * K <= c.TwoJ() + d.TwoJ()))
    {
        RadialFunction density = c.GetDensity(d);
        density.resize(integrator->GetLattice()->size());

        coulomb->GetPotential(K, density, potential);
        return true;
    }
    else
    {   potential.Clear();
        return false;
    }
}

HartreeY* HartreeY::Clone() const
{
    HartreeY* copy = new HartreeY(integrator, coulomb);
    copy->K = K;
    copy->potential = potential;

    return copy;
}

void HartreeY::Alert()
{
    if(isZero())
        return;

    unsigned int i = potential.size();
    potential.resize(lattice->size());

    // If potential has grown, add points
    if(i < potential.size())
    {
        const double* R_k = lattice->Rpower(K);
        const double* R_kplusone = lattice->Rpower(K+1);
        double charge = potential.f[i-1] * R_k[i-1];

        while(i < potential.size())
        {
            potential.f[i] = charge/R_k[i];
            potential.dfdr[i] = - K * charge/R_kplusone[i];
            i++;
        }
    }
}

double HartreeY::GetMatrixElement(const Orbital& b, const SingleParticleWavefunction& a, bool reverse) const
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
