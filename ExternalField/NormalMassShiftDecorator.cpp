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

void NormalMassShiftDecorator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    BaseDecorator::GetODEFunction(latticepoint, fg, w);

    if(do_nonrel_nms)
    {
        const double alpha = physicalConstant->GetAlpha();
        w[0] += -2.*scale/alpha * fg.g[latticepoint];
    }
}

void NormalMassShiftDecorator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    BaseDecorator::GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);

    if(do_nonrel_nms)
    {
        const double alpha = physicalConstant->GetAlpha();
        w_g[0] += -2.*scale/alpha;
    }
}

void NormalMassShiftDecorator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    BaseDecorator::GetODEJacobian(latticepoint, fg, jacobian, dwdr);

    if(do_nonrel_nms)
    {
        const double alpha = physicalConstant->GetAlpha();
        jacobian[0][1] = -2.*scale/alpha;
    }
}

SpinorFunction NormalMassShiftDecorator::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta = BaseDecorator::ApplyTo(a);

    if(do_nonrel_nms)
    {
        const double alphasquared = physicalConstant->GetAlphaSquared();
        for(unsigned int i = 0; i < ta.size(); i++)
        {
            ta.g[i] += 2.*scale/alphasquared * a.g[i];
            ta.dgdr[i] += 2.*scale/alphasquared * a.dgdr[i];
        }
    }

    return ta;
}

SpinorFunction NormalMassShiftDecorator::CalculateExtraExchange(const SpinorFunction& s) const
{
    SpinorFunction exchange(s.Kappa(), s.size());

    std::vector<double> second_derivative_f(s.size());
    std::vector<double> second_derivative_g(s.size());

    const double* R = lattice->R();
    const double* R2 = lattice->Rpower(2);
    const double alpha = physicalConstant->GetAlpha();
    const double alphasquared = physicalConstant->GetAlphaSquared();

    differentiator->GetDerivative(s.dfdr, second_derivative_f);
    differentiator->GetDerivative(s.dgdr, second_derivative_g);

    if(false)
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
            exchange.g[i] += 0.5 * (second_derivative_g[i] - coeff_g * s.g[i]/R2[i]) + 2./alphasquared * s.g[i];
            exchange.dgdr[i] += 0.5 * (third_derivative_g[i] - coeff_g * (s.dgdr[i] - 2.*s.g[i]/R[i]) /R2[i]) + 2./alphasquared * s.dgdr[i];
        }
    }

    if(do_rel_nms)
    {
        double ZalphaOnTwoM = 0.5 * Z * alpha;

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
