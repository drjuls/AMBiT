#include "FieldShift.h"
#include "Include.h"
#include "HartreeFock/NucleusDecorator.h"

namespace Ambit
{
FieldShiftDecorator::FieldShiftDecorator(pHFOperator wrapped_hf, double R0, double dR, double t, pIntegrator integration_strategy):
    BaseDecorator(wrapped_hf, integration_strategy)
{
    pODESolver ode_solver = std::make_shared<AdamsSolver>(integrator);
    pCoulombOperator coulomb = std::make_shared<CoulombOperator>(lattice, ode_solver);
    pHFOperator zero = std::make_shared<HFOperatorBase>(Z, core, physicalConstant, integrator);

    NucleusDecorator nuc(zero, coulomb, integrator);
    nuc.SetFermiParameters(R0, t);
    RadialFunction pot0 = nuc.GetDirectPotential();
    double Rrms0 = nuc.CalculateNuclearRMSRadius();

    nuc.SetFermiParameters(R0 + dR, t);
    directPotential = nuc.GetDirectPotential() - pot0;
    dRrms = nuc.CalculateNuclearRMSRadius() - Rrms0;

    *outstream << "Change in nuclear RMS radius = " << std::setprecision(6) << dRrms << std::endl;
}

FieldShiftCalculator::FieldShiftCalculator(MultirunOptions& user_input, Atom& atom):
    TransitionCalculator(user_input, atom)
{
    pHFOperatorConst hf = atom.GetHFOperator();
    pHFOperator zero = std::make_shared<HFOperatorBase>(*hf);

    // Nuclear parameters
    double R0 = 0.;
    double t = 0.;
    pNucleusDecorator nuc = atom.GetNucleusDecorator();
    if(nuc)
    {
        R0 = nuc->GetNuclearRadius();
        t = nuc->GetNuclearThickness();
    }

    // Change in nuclear radius
    double dR = user_input("DeltaNuclearRadius", 0.1);
    std::shared_ptr<FieldShiftDecorator> field_shift = std::make_shared<FieldShiftDecorator>(zero, R0, dR, t);

    double scale = user_input("Scale", 1.);
    field_shift->SetScale(scale);

    if(user_input.search("--rpa"))
        op = MakeRPA(field_shift, hf, atom.GetHartreeY());
    else
        op = field_shift;

    // Get conversion factor
    MathConstant* math = MathConstant::Instance();
    double Rrms = nuc->CalculateNuclearRMSRadius();
    double dRrms = field_shift->GetDeltaRMSRadius();
    MHzOverRsq = math->AtomicFrequencyMHz()/(gsl_pow_2(Rrms + dRrms) - gsl_pow_2(Rrms));
}

void FieldShiftCalculator::PrintHeader() const
{
    *outstream << "Field shift in MHz/fm^2: " << std::endl;
}

void FieldShiftCalculator::PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const
{
    double value = matrix_element * MHzOverRsq;

    *outstream << "  " << Name(left) << " -> " << Name(right)
               << " = " << std::setprecision(6) << value << std::endl;
}
}
