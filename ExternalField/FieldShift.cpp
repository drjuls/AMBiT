#include "FieldShift.h"
#include "Include.h"
#include "HartreeFock/NucleusDecorator.h"

namespace Ambit
{
FieldShiftDecorator::FieldShiftDecorator(pHFOperator wrapped_hf, double R0, double dR, double t, double dt, pIntegrator integration_strategy):
    BaseDecorator(wrapped_hf, integration_strategy)
{
    pODESolver ode_solver = std::make_shared<AdamsSolver>(integrator);
    pCoulombOperator coulomb = std::make_shared<CoulombOperator>(lattice, ode_solver);
    pHFOperator zero = std::make_shared<HFOperatorBase>(Z, core, physicalConstant, integrator);

    NucleusDecorator nuc(zero, coulomb, integrator);
    nuc.SetFermiParameters(R0, t);
    RadialFunction pot0 = nuc.GetDirectPotential();
    double Rrms0 = nuc.CalculateNuclearRMSRadius();
    double R40 = nuc.CalculateNuclearR4();

    nuc.SetFermiParameters(R0 + dR, t + dt);
    directPotential = nuc.GetDirectPotential() - pot0;
    dRrms = nuc.CalculateNuclearRMSRadius() - Rrms0;
    dR4   = nuc.CalculateNuclearR4() - R40;

    *outstream << "Change in nuclear RMS radius = " << std::setprecision(6) << dRrms << std::endl;
    if(dt)
        *outstream << "Change in nuclear R4         = " << std::setprecision(6) << dR4 << std::endl;
}

FieldShiftCalculator::FieldShiftCalculator(MultirunOptions& user_input, Atom& atom):
    TransitionCalculator(user_input, atom)
{
    pHFOperatorConst hf = atom.GetHFOperator();
    pHFOperator zero = std::make_shared<HFOperatorBase>(*hf);

    // Nuclear parameters
    double R0 = 0.;
    double t = 2.3;
    pNucleusDecorator nuc = atom.GetNucleusDecorator();
    if(nuc)
    {
        R0 = nuc->GetNuclearRadius();
        t = nuc->GetNuclearThickness();
    }

    // Change in nuclear radius: default 0.1
    double dR = user_input("DeltaNuclearRadius", 0.1);
    double dt = user_input("DeltaNuclearThickness", 0.0);
    std::shared_ptr<FieldShiftDecorator> field_shift = std::make_shared<FieldShiftDecorator>(zero, R0, dR, t, dt);

    double scale = user_input("Scale", 1.);
    field_shift->SetScale(scale);

    if(user_input.search("--rpa"))
        op = MakeRPA(field_shift, hf, atom.GetHartreeY());
    else
        op = field_shift;

    // Get conversion factor
    double Rrms = nuc->CalculateNuclearRMSRadius();
    double dRrms = field_shift->GetDeltaRMSRadius();
    dR2 = (gsl_pow_2(Rrms + dRrms) - gsl_pow_2(Rrms));
    dR4 = field_shift->GetDeltaR4();
}

void FieldShiftCalculator::PrintHeader() const
{
    *outstream << "Field shift in MHz (dR2 = " << std::setprecision(6) << dR2
               << ", dR4 = " << dR4 << "):\n";
}

void FieldShiftCalculator::PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const
{
    MathConstant* math = MathConstant::Instance();
    double value = matrix_element * math->AtomicFrequencyMHz();

    *outstream << "  " << Name(left) << " -> " << Name(right)
               << " = " << std::setprecision(6) << value << std::endl;
}
}
