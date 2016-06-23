#include "NormalMassShiftDecorator.h"
#include "Include.h"
#include "Universal/MathConstant.h"
#include "Universal/PhysicalConstant.h"

NormalMassShiftDecorator::NormalMassShiftDecorator(pHFOperator wrapped_hf, bool only_rel_nms, bool nonrel):
    BaseDecorator(wrapped_hf), lambda(0.0)
{
    do_nonrel_nms = !only_rel_nms;
    do_rel_nms = !nonrel;
}

void NormalMassShiftDecorator::Alert()
{
    if(currentExchangePotential.size() > lattice->size())
        currentExchangePotential.resize(lattice->size());
}

/** Set exchange (nonlocal) potential and energy for ODE routines. */
void NormalMassShiftDecorator::SetODEParameters(const Orbital& approximation)
{
    HFOperatorDecorator::SetODEParameters(approximation);
    currentExchangePotential = CalculateExtraExchange(approximation);
}

/** Get exchange (nonlocal) potential. */
SpinorFunction NormalMassShiftDecorator::GetExchange(pOrbitalConst approximation) const
{
    SpinorFunction ret = wrapped->GetExchange(approximation);

    if(approximation == NULL)
        ret += currentExchangePotential;
    else
        ret += CalculateExtraExchange(*approximation);

    return ret;
}

void NormalMassShiftDecorator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    wrapped->GetODEFunction(latticepoint, fg, w);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha();
        w[0] += alpha * currentExchangePotential.g[latticepoint];
        w[1] -= alpha * currentExchangePotential.f[latticepoint];
    }
}
void NormalMassShiftDecorator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    wrapped->GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha();
        w_const[0] += alpha * currentExchangePotential.g[latticepoint];
        w_const[1] -= alpha * currentExchangePotential.f[latticepoint];
    }
}
void NormalMassShiftDecorator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    wrapped->GetODEJacobian(latticepoint, fg, jacobian, dwdr);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha();
        dwdr[0] += alpha * currentExchangePotential.dgdr[latticepoint];
        dwdr[1] -= alpha * currentExchangePotential.dfdr[latticepoint];
    }
}

SpinorFunction NormalMassShiftDecorator::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta = wrapped->ApplyTo(a);
    ta -= CalculateExtraExchange(a);

    return ta;
}

SpinorFunction NormalMassShiftDecorator::CalculateExtraExchange(const SpinorFunction& s) const
{
    SpinorFunction exchange(s.Kappa(), s.size());

    if(lambda == 0)
        return exchange;

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
            exchange.f[i] += 0.5 * lambda * (second_derivative_f[i] - coeff_f * s.f[i]/R2[i]);
            exchange.dfdr[i] += 0.5 * lambda * (third_derivative_f[i] - coeff_f * (s.dfdr[i] - 2.*s.f[i]/R[i]) /R2[i]);
            exchange.g[i] += 0.5 * lambda * (second_derivative_g[i] - coeff_g * s.g[i]/R2[i]);
            exchange.dgdr[i] += 0.5 * lambda * (third_derivative_g[i] - coeff_g * (s.dgdr[i] - 2.*s.g[i]/R[i]) /R2[i]);
        }
    }

    if(do_rel_nms)
    {
        double ZalphaOnTwoM = 0.5 * lambda * Z * physicalConstant->GetAlpha();

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
