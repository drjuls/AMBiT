#include "Hyperfine.h"
#include "Include.h"
#include "Universal/MathConstant.h"

HyperfineDipoleOperator::HyperfineDipoleOperator(pIntegrator integration_strategy, double nuclear_magnetic_radius_fm):
    SpinorOperator(1, Parity::even, integration_strategy), lattice(integration_strategy->GetLattice())
{
    nuclear_radius = nuclear_magnetic_radius_fm/MathConstant::Instance()->BohrRadiusInFermi();
    nuclear_radius_lattice = lattice->real_to_lattice(nuclear_radius);
}

SpinorFunction HyperfineDipoleOperator::ReducedApplyTo(const SpinorFunction& a, int kappa_b) const
{
    SpinorFunction ret(kappa_b, a.size());
    MathConstant* math = MathConstant::Instance();

    double prefactor = -(a.Kappa() + kappa_b) * math->SphericalTensorReducedMatrixElement(-kappa_b, a.Kappa(), 1);
    if(!prefactor)
        return ret;
    prefactor *= math->NuclearMagneton();

    const double* R = lattice->R();
    const double* R2 = lattice->Rpower(2);
    const double* R3 = lattice->Rpower(3);

    double Rn3 = gsl_pow_3(nuclear_radius);
    unsigned int i;
    for(i = 0; i < nuclear_radius_lattice; i++)
    {
        ret.f[i] = R[i] * a.g[i]/Rn3;
        ret.dfdr[i] = (R[i] * a.dgdr[i] + a.g[i])/Rn3;
        ret.g[i] = R[i] * a.f[i]/Rn3;
        ret.dgdr[i] = (R[i] * a.dfdr[i] + a.f[i])/Rn3;
    }
    for(i = nuclear_radius_lattice; i < ret.size(); i++)
    {
        ret.f[i] = a.g[i]/R2[i];
        ret.dfdr[i] = a.dgdr[i]/R2[i] - 2. * a.g[i]/R3[i];
        ret.g[i] = a.f[i]/R2[i];
        ret.dgdr[i] = a.dfdr[i]/R2[i] - 2. * a.f[i]/R3[i];
    }

    return ret * prefactor;
}

HyperfineDipoleCalculator::HyperfineDipoleCalculator(MultirunOptions& user_input, Atom& atom):
    TransitionCalculator(user_input, atom.GetBasis(), atom.GetLevels())
{
    pHFOperatorConst hf = atom.GetHFOperator();

    double magnetic_radius = user_input("NuclearMagneticRadius", 0.0);
    op = std::make_shared<HyperfineDipoleOperator>(hf->GetIntegrator(), magnetic_radius);

    g_I = user_input("gOnI", 1.0);

    if(user_input.search("--rpa"))
        op = MakeStaticRPA(std::static_pointer_cast<HyperfineDipoleOperator>(op), hf, atom.GetHartreeY());
}

void HyperfineDipoleCalculator::PrintHeader() const
{
    *outstream << "Hyperfine dipole matrix elements (stretched states) in MHz: " << std::endl;
}

void HyperfineDipoleCalculator::PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const
{
    MathConstant* math = MathConstant::Instance();
    double value = matrix_element * g_I * math->AtomicFrequencyMHz();
    if(left == right)
        value = value/left.first->GetJ();

    *outstream << "  " << Name(left) << " -> " << Name(right)
               << " = " << std::setprecision(6) << value << std::endl;
}
