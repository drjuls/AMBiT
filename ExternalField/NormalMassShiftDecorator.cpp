#include "NormalMassShiftDecorator.h"
#include "Include.h"
#include "Universal/MathConstant.h"
#include "Universal/PhysicalConstant.h"

namespace Ambit
{
NormalMassShiftDecorator::NormalMassShiftDecorator(pHFOperator wrapped_hf, bool only_rel_nms, bool nonrel):
    BaseDecorator(wrapped_hf)
{
    scale = 0.0;
    do_nonrel_nms = !only_rel_nms;
    do_rel_nms = !nonrel;
}

SpinorFunction NormalMassShiftDecorator::CalculateExtraExchange(const SpinorFunction& s) const
{
    SpinorFunction exchange(s.Kappa(), s.size());

    std::vector<double> second_derivative_f(s.size());
    std::vector<double> second_derivative_g(s.size());

    const double* R = lattice->R();
    const double* R2 = lattice->Rpower(2);

    differentiator->GetDerivative(s.dfdr, second_derivative_f);
    differentiator->GetDerivative(s.dgdr, second_derivative_g);

    if(do_nonrel_nms)
    {
        std::vector<double> third_derivative_f(s.size());
        std::vector<double> third_derivative_g(s.size());

        differentiator->GetSecondDerivative(s.dfdr, third_derivative_f);
        differentiator->GetSecondDerivative(s.dgdr, third_derivative_g);

        double coeff_f = s.Kappa() * (s.Kappa() + 1);
        double coeff_g = s.Kappa() * (s.Kappa() - 1);

        for(int i = 0; i < s.size(); i++)
        {
            exchange.f[i] += 0.5 * (second_derivative_f[i] - coeff_f * s.f[i]/R2[i]);
            exchange.dfdr[i] += 0.5 * (third_derivative_f[i] - coeff_f * (s.dfdr[i] - 2.*s.f[i]/R[i]) /R2[i]);
            exchange.g[i] += 0.5 * (second_derivative_g[i] - coeff_g * s.g[i]/R2[i]);
            exchange.dgdr[i] += 0.5 * (third_derivative_g[i] - coeff_g * (s.dgdr[i] - 2.*s.g[i]/R[i]) /R2[i]);
        }
    }

    if(do_rel_nms)
    {
        double ZalphaOnTwoM = 0.5 * Z * physicalConstant->GetAlpha();

        for(int i = 0; i < s.size(); i++)
        {
            exchange.f[i] -= ZalphaOnTwoM/R[i] * (2. * s.dgdr[i] - (s.Kappa()+1) * s.g[i]/R[i]);
            exchange.dfdr[i] -= ZalphaOnTwoM/R[i]
                * (2. * second_derivative_g[i] - (s.Kappa()+3) * s.dgdr[i]/R[i] + (2.*s.Kappa()+2) * s.g[i]/R2[i]);

            exchange.g[i] += ZalphaOnTwoM * (2. * s.dfdr[i] + (s.Kappa() - 1) * s.f[i]/R[i])/R[i];
            exchange.dgdr[i] += ZalphaOnTwoM/R[i]
                * (2. * second_derivative_f[i] + (s.Kappa()-3) * s.dfdr[i]/R[i] - (2.*s.Kappa()-2) * s.f[i]/R2[i]);
        }
    }

    return exchange;
}
}
