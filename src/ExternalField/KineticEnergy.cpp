#include "KineticEnergy.h"

namespace Ambit
{
SpinorFunction KineticEnergyDecorator::CalculateExtraExchange(const SpinorFunction& s) const
{
    SpinorFunction exchange(s.Kappa(), s.size());

    std::vector<double> second_derivative_f(s.size());
    std::vector<double> second_derivative_g(s.size());

    const double* R = lattice->R();
    const double* R2 = lattice->Rpower(2);

    differentiator->GetDerivative(s.dfdr, second_derivative_f);
    differentiator->GetDerivative(s.dgdr, second_derivative_g);

    double c = physicalConstant->GetSpeedOfLight()/2.;

    for(int i = 0; i < s.size(); i++)
    {
        exchange.f[i] += c * (s.dgdr[i] - s.Kappa() * s.g[i]/R[i]);
        exchange.dfdr[i] += c * (second_derivative_g[i] - s.Kappa() * s.dgdr[i]/R[i] + s.Kappa() * s.g[i]/R2[i]);

        exchange.g[i] -= c * (s.dfdr[i] + s.Kappa() * s.f[i]/R[i]);
        exchange.dgdr[i] -= c * (second_derivative_f[i] + s.Kappa() * s.dfdr[i]/R[i] - s.Kappa() * s.f[i]/R2[i]);
    }

    return exchange;
}

KineticEnergyCalculator::KineticEnergyCalculator(MultirunOptions& user_input, Atom& atom):
    TransitionCalculator(user_input, atom)
{
    pHFOperatorConst hf = atom.GetHFOperator();
    pHFOperator zero = std::make_shared<HFOperatorBase>(*hf);
    auto rke = std::make_shared<KineticEnergyDecorator>(zero);

    double scale = user_input("Scale", 1.);
    rke->SetScale(scale);

    if(user_input.search("--rpa"))
        op = MakeRPA(rke, hf, atom.GetHartreeY());
    else
        op = rke;
}

void KineticEnergyCalculator::PrintHeader() const
{
    *outstream << "Kinetic Energy (a.u.): " << std::endl;
}

void KineticEnergyCalculator::PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const
{
    double value = matrix_element;

    *outstream << "  " << Name(left) << " -> " << Name(right)
               << " = " << std::setprecision(6) << value << std::endl;
}

}
