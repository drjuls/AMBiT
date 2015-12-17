#include "Hyperfine.h"
#include "Include.h"
#include "Universal/MathConstant.h"

HyperfineDipoleOperator::HyperfineDipoleOperator(pLattice lattice, pOPIntegrator integration_strategy):
    SpinorOperator(1, integration_strategy), lattice(lattice)
{}

SpinorFunction HyperfineDipoleOperator::ApplyTo(const SpinorFunction& a, int kappa_b) const
{
    SpinorFunction ret(kappa_b, a.size());
    MathConstant* math = MathConstant::Instance();

    // Atomic units
    double prefactor = double(a.Kappa())/(a.J() * (a.J() + 1.)) * math->NuclearMagneton();

    const double* R2 = lattice->Rpower(2);
    const double* R3 = lattice->Rpower(3);

    for(unsigned int i = 0; i < ret.size(); i++)
    {
        ret.f[i] = a.g[i]/R2[i];
        ret.dfdr[i] = a.dgdr[i]/R2[i] - 2. * a.g[i]/R3[i];
        ret.g[i] = a.f[i]/R2[i];
        ret.dgdr[i] = a.dfdr[i]/R2[i] - 2. * a.g[i]/R3[i];
    }

    return ret * prefactor;
}
