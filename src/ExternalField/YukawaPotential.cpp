#include "YukawaPotential.h"

namespace Ambit
{
YukawaDecorator::YukawaDecorator(pHFOperator wrapped_hf, double mass, double scale, pIntegrator integration_strategy):
    BaseDecorator(wrapped_hf, integration_strategy), mass(mass)
{
    SetScale(scale);
    GenerateYukawaPotential(mass);
}

void YukawaDecorator::GenerateYukawaPotential(double new_mass)
{
    mass = new_mass;
    double r_c = physicalConstant->GetAlpha()/mass;
    directPotential.resize(lattice->size());
    auto R = lattice->R();

    for(unsigned int i = 0; i < directPotential.size(); i++)
    {
        directPotential.f[i] = exp(-R[i]/r_c)/R[i];
        directPotential.dfdr[i] = directPotential.f[i] * (-1./R[i] - 1/r_c);
    }
}

void YukawaDecorator::Alert()
{
    GenerateYukawaPotential(mass);
}

YukawaCalculator::YukawaCalculator(MultirunOptions& user_input, Atom& atom):
    TransitionCalculator(user_input, atom)
{
    pHFOperatorConst hf = atom.GetHFOperator();
    pHFOperator zero = std::make_shared<HFOperatorBase>(*hf);

    if(user_input.VariableExists("Mass"))
        mass = user_input("Mass", 1.0);
    else if(user_input.VariableExists("MassEV"))
        mass = user_input("MassEV", 1.0)/MathConstant::Instance()->ElectronMassInEV;
    else if(user_input.VariableExists("Rc"))
        mass = 1./(atom.GetPhysicalConstants()->GetAlpha() * user_input("Rc", 1.0));

    double scale = user_input("Scale", 1.);

    auto yukawa_decorator = std::make_shared<YukawaDecorator>(zero, mass, scale, hf->GetIntegrator());

    if(user_input.search("--rpa"))
        op = MakeRPA(yukawa_decorator, hf, atom.GetHartreeY());
    else
        op = yukawa_decorator;
}

void YukawaCalculator::PrintHeader() const
{
    *outstream << "Yukawa shift (a.u.):" << std::endl;
}

void YukawaCalculator::PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const
{
    double value = matrix_element;

    *outstream << "  " << Name(left) << " -> " << Name(right)
               << " = " << std::setprecision(12) << value << std::endl;
}
}
