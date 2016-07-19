#include "ExchangeDecorator.h"
#include "Include.h"
#include "Universal/MathConstant.h"

void ExchangeDecorator::Alert()
{
    if(currentExchangePotential.size() > lattice->size())
        currentExchangePotential.resize(lattice->size());
}

/** Set exchange (nonlocal) potential and energy for ODE routines. */
void ExchangeDecorator::SetODEParameters(const Orbital& approximation)
{
    HFOperatorDecorator::SetODEParameters(approximation);
    currentExchangePotential = CalculateExtraExchange(approximation);
}

/** Get exchange (nonlocal) potential. */
SpinorFunction ExchangeDecorator::GetExchange(pOrbitalConst approximation) const
{
    SpinorFunction ret = wrapped->GetExchange(approximation);

    if(approximation == nullptr)
        ret += currentExchangePotential * scale;
    else
        ret += CalculateExtraExchange(*approximation) * scale;

    return ret;
}

void ExchangeDecorator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    wrapped->GetODEFunction(latticepoint, fg, w);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha() * scale;
        w[0] += alpha * currentExchangePotential.g[latticepoint];
        w[1] -= alpha * currentExchangePotential.f[latticepoint];
    }
}

void ExchangeDecorator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    wrapped->GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha() * scale;
        w_const[0] += alpha * currentExchangePotential.g[latticepoint];
        w_const[1] -= alpha * currentExchangePotential.f[latticepoint];
    }
}

void ExchangeDecorator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    wrapped->GetODEJacobian(latticepoint, fg, jacobian, dwdr);

    if(include_nonlocal && latticepoint < currentExchangePotential.size())
    {   double alpha = physicalConstant->GetAlpha() * scale;
        dwdr[0] += alpha * currentExchangePotential.dgdr[latticepoint];
        dwdr[1] -= alpha * currentExchangePotential.dfdr[latticepoint];
    }
}

SpinorFunction ExchangeDecorator::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta = wrapped->ApplyTo(a);
    ta -= CalculateExtraExchange(a) * scale;

    return ta;
}
