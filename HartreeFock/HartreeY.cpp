#include "HartreeY.h"

HartreeY::HartreeY(pIntegrator integration_strategy, pCoulombOperator coulomb):
    HartreeYBase(integration_strategy), LatticeObserver(integration_strategy->GetLattice()), coulomb(coulomb)
{
    two_body_reverse_symmetry_exists = true;
}

bool HartreeY::SetParameters(int new_K, pSpinorFunctionConst new_c, pSpinorFunctionConst new_d)
{
    K = new_K;
    c = new_c;
    d = new_d;

    if(((K + c->L() + d->L())%2 == 0) &&
       (abs(c->TwoJ() - d->TwoJ()) <= 2 * K) &&
       (2 * K <= c->TwoJ() + d->TwoJ()))
    {
        RadialFunction density = c->GetDensity(*d);
        density.resize(integrator->GetLattice()->size());

        coulomb->GetPotential(K, density, potential);
        return true;
    }
    else
    {   potential.Clear();
        return false;
    }
}

int HartreeY::SetOrbitals(pSpinorFunctionConst new_c, pSpinorFunctionConst new_d)
{
    c = new_c;
    d = new_d;

    K = GetMinK();
    if((K != -1) && (parent == nullptr))
        SetK(K);

    return K;
}

int HartreeY::GetMinK() const
{
    int Knew = abs(c->TwoJ() - d->TwoJ())/2;
    if((Knew + c->L() + d->L())%2 == 1)
        Knew++;

    return Knew;
}

int HartreeY::NextK()
{
    if(K != -1 && c && d)
    {   K += 2;
        if(SetParameters(K, c, d))
            return K;
    }

    K = -1;
    return K;
}

int HartreeY::GetMaxK() const
{
    int maxK = -1;
    if(c && d)
    {
        maxK = (c->TwoJ() + d->TwoJ())/2;
        if((maxK + c->L() + d->L())%2 == 1)
            maxK--;
    }

    return maxK;
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
        const double* R_kplusone = lattice->Rpower(K+1);
        const double* R_kplustwo = lattice->Rpower(K+2);
        double charge = potential.f[i-1] * R_kplusone[i-1];

        while(i < potential.size())
        {
            potential.f[i] = charge/R_kplusone[i];
            potential.dfdr[i] = - K * charge/R_kplustwo[i];
            i++;
        }
    }
}

double HartreeY::GetMatrixElement(const Orbital& b, const Orbital& a, bool reverse) const
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
