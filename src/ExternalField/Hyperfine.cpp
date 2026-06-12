#include "Hyperfine.h"
#include "Include.h"
#include "Specification/TransitionSpecifications.h"
#include "Universal/MathConstant.h"

namespace Ambit
{
HyperfineEJOperator::HyperfineEJOperator(int J, pIntegrator integration_strategy, double nuclear_radius_fm):
    SpinorOperator(J, (J%2? Parity::odd: Parity::even), integration_strategy), lattice(integration_strategy->GetLattice())
{
    nuclear_radius = nuclear_radius_fm/MathConstant::Instance()->BohrRadiusInFermi();
    nuclear_radius_lattice = lattice->real_to_lattice(nuclear_radius);
}

SpinorFunction HyperfineEJOperator::ReducedApplyTo(const SpinorFunction& a, int kappa_b) const
{
    SpinorFunction ret(kappa_b, a.size());
    MathConstant* math = MathConstant::Instance();

    double prefactor = -math->SphericalTensorReducedMatrixElement(kappa_b, a.Kappa(), K);
    if(!prefactor)
        return ret;

    const double* R = lattice->R();
    const double* RKp1 = lattice->Rpower(K+1);
    const double* RKp2 = lattice->Rpower(K+2);

    double RnKp2 = gsl_pow_int(nuclear_radius, K+2);
    unsigned int i;
    for(i = 0; i < nuclear_radius_lattice; i++)
    {
        ret.f[i] = R[i] * a.f[i]/RnKp2;
        ret.dfdr[i] = (R[i] * a.dfdr[i] + a.f[i])/RnKp2;
        ret.g[i] = R[i] * a.g[i]/RnKp2;
        ret.dgdr[i] = (R[i] * a.dgdr[i] + a.g[i])/RnKp2;
    }
    for(i = nuclear_radius_lattice; i < ret.size(); i++)
    {
        ret.f[i] = a.f[i]/RKp1[i];
        ret.dfdr[i] = a.dfdr[i]/RKp1[i] - (K+1) * a.f[i]/RKp2[i];
        ret.g[i] = a.g[i]/RKp1[i];
        ret.dgdr[i] = a.dgdr[i]/RKp1[i] - (K+1) * a.g[i]/RKp2[i];
    }

    return ret * prefactor;
}

HyperfineMJOperator::HyperfineMJOperator(int J, pIntegrator integration_strategy, double nuclear_radius_fm):
    SpinorOperator(J, (J%2? Parity::even: Parity::odd), integration_strategy), lattice(integration_strategy->GetLattice())
{
    nuclear_radius = nuclear_radius_fm/MathConstant::Instance()->BohrRadiusInFermi();
    nuclear_radius_lattice = lattice->real_to_lattice(nuclear_radius);
}

SpinorFunction HyperfineMJOperator::ReducedApplyTo(const SpinorFunction& a, int kappa_b) const
{
    SpinorFunction ret(kappa_b, a.size());
    MathConstant* math = MathConstant::Instance();

    double prefactor = -(a.Kappa() + kappa_b) * math->SphericalTensorReducedMatrixElement(-kappa_b, a.Kappa(), K)/K;
    if(!prefactor)
        return ret;

    const double* R = lattice->R();
    const double* RKp1 = lattice->Rpower(K+1);
    const double* RKp2 = lattice->Rpower(K+2);

    double RnKp2 = gsl_pow_int(nuclear_radius, K+2);
    unsigned int i;
    for(i = 0; i < nuclear_radius_lattice; i++)
    {
        ret.f[i] = R[i] * a.g[i]/RnKp2;
        ret.dfdr[i] = (R[i] * a.dgdr[i] + a.g[i])/RnKp2;
        ret.g[i] = R[i] * a.f[i]/RnKp2;
        ret.dgdr[i] = (R[i] * a.dfdr[i] + a.f[i])/RnKp2;
    }
    for(i = nuclear_radius_lattice; i < ret.size(); i++)
    {
        ret.f[i] = a.g[i]/RKp1[i];
        ret.dfdr[i] = a.dgdr[i]/RKp1[i] - (K+1) * a.g[i]/RKp2[i];
        ret.g[i] = a.f[i]/RKp1[i];
        ret.dgdr[i] = a.dfdr[i]/RKp1[i] - (K+1) * a.f[i]/RKp2[i];
    }

    return ret * prefactor;
}

HyperfineDipoleCalculator::HyperfineDipoleCalculator(Atom& atom, TransitionsHFS1Specification& specification):
    TransitionCalculator(atom.GetBasis(), atom.GetLevels(), specification)
{
    pHFOperatorConst hf = atom.GetHFOperator();

    //double magnetic_radius = user_input("NuclearMagneticRadius", -1.0);
    double magnetic_radius;
    auto radius_opt = specification.nuclear_magnetic_radius;
    if(radius_opt)
    {
        magnetic_radius = radius_opt.value();
    }
    else
    {
        pNucleusDecorator nuc = atom.GetNucleusDecorator();
        if(nuc)
            magnetic_radius = nuc->CalculateNuclearRMSRadius() * sqrt(5./3.);
        else
            magnetic_radius = 0.;
    }
    op = std::make_shared<HyperfineMJOperator>(1, hf->GetIntegrator(), magnetic_radius);

    g_I = specification.g_on_i;

    if(specification.do_rpa)
        op = MakeRPA(std::static_pointer_cast<HyperfineMJOperator>(op), hf, atom.GetHartreeY());
}

void HyperfineDipoleCalculator::PrintHeader() const
{
    *outstream << "Hyperfine dipole matrix elements (stretched states) in MHz: " << std::endl;
}

void HyperfineDipoleCalculator::PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const
{
    MathConstant* math = MathConstant::Instance();
    double value = matrix_element * g_I * math->NuclearMagneton() * math->AtomicFrequencyMHz();
    if(left == right)
        value = value/left.first->GetJ();

    *outstream << "  " << Name(left) << " -> " << Name(right)
               << " = " << std::setprecision(6) << value << std::endl;
}

HyperfineQuadrupoleCalculator::HyperfineQuadrupoleCalculator(Atom& atom, TransitionsHFS2Specification& specification):
    TransitionCalculator(atom.GetBasis(), atom.GetLevels(), specification)
{
    pHFOperatorConst hf = atom.GetHFOperator();

    double magnetic_radius;
    auto radius_opt = specification.nuclear_quadrupole_radius;
    if(radius_opt)
    {
        magnetic_radius = radius_opt.value();
    }
    else
    {
        pNucleusDecorator nuc = atom.GetNucleusDecorator();
        if(nuc)
            magnetic_radius = nuc->CalculateNuclearRMSRadius() * sqrt(5./3.);
        else
            magnetic_radius = 0.;
    }
    op = std::make_shared<HyperfineEJOperator>(2, hf->GetIntegrator(), magnetic_radius);

    Q = specification.q;

    if(specification.do_rpa)
        op = MakeRPA(std::static_pointer_cast<HyperfineEJOperator>(op), hf, atom.GetHartreeY());
}

void HyperfineQuadrupoleCalculator::PrintHeader() const
{
    *outstream << "Hyperfine quadrupole matrix elements (stretched states) in MHz: " << std::endl;
}

void HyperfineQuadrupoleCalculator::PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const
{
    MathConstant* math = MathConstant::Instance();
    double value = 2. * Q * matrix_element * math->Barn() * math->AtomicFrequencyMHz();

    *outstream << "  " << Name(left) << " -> " << Name(right)
               << " = " << std::setprecision(6) << value << std::endl;
}

GeneralisedHyperfineCalculator::GeneralisedHyperfineCalculator(Atom& atom, TransitionsHFISpecification& specification):
    TransitionCalculator(atom.GetBasis(), atom.GetLevels(), specification)
{
    MathConstant* math = MathConstant::Instance();
    pHFOperatorConst hf = atom.GetHFOperator();

    double magnetic_radius;
    auto radius_opt = specification.nuclear_radius;
    if(radius_opt)
    {
        magnetic_radius = radius_opt.value();
    }
    else
    {
        pNucleusDecorator nuc = atom.GetNucleusDecorator();
        if(nuc)
            magnetic_radius = nuc->CalculateNuclearRMSRadius() * sqrt(5./3.);
        else
            magnetic_radius = 0.;
    }

    std::string type = specification.op.value_or("");
    int k = 0;
    if(type.size() >= 2 && (type[0] == 'E' || type[0] == 'M'))
    {
        k = atoi(type.c_str() + 1);
        if(k < 1 || k > 9)
            k = 0;

        if(type[0] == 'E')
            EorM = MultipolarityType::E;
        else
            EorM = MultipolarityType::M;
    }

    if(k == 0)
    {   *errstream << "HFI/Operator: Bad string " << type << " (should be EJ or MJ with J a nonnegative integer)." << std::endl;
        exit(1);
    }

    if(EorM == MultipolarityType::E)
    {
        op = std::make_shared<HyperfineEJOperator>(k, hf->GetIntegrator(), magnetic_radius);
        if(specification.do_rpa)
            op = MakeRPA(std::static_pointer_cast<HyperfineEJOperator>(op), hf, atom.GetHartreeY());
    }
    else
    {
        op = std::make_shared<HyperfineMJOperator>(k, hf->GetIntegrator(), magnetic_radius);
        if(specification.do_rpa)
            op = MakeRPA(std::static_pointer_cast<HyperfineMJOperator>(op), hf, atom.GetHartreeY());
    }
}

void GeneralisedHyperfineCalculator::PrintHeader() const
{
    if(specification.reduced_elements)
        *outstream << "Hyperfine-" << Name(EorM) << op->GetK() << " reduced matrix elements (a.u.): " << std::endl;
    else
        *outstream << "Hyperfine-" << Name(EorM) << op->GetK() << " matrix elements (stretched states) in a.u.: " << std::endl;
}

void GeneralisedHyperfineCalculator::PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const
{
    MathConstant* math = MathConstant::Instance();
    double value = matrix_element;
    if(specification.reduced_elements)
    {   int twoj1 = left.first->GetTwoJ();
        int twoj2 = right.first->GetTwoJ();
        value = value/math->Electron3j(twoj2, twoj1, op->GetK(), twoj2, -twoj1);
    }

    *outstream << "  " << Name(left) << " -> " << Name(right)
               << " = " << std::setprecision(6) << value << std::endl;
}
}
