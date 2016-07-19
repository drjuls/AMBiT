#ifndef EXCHANGE_DECORATOR_H
#define EXCHANGE_DECORATOR_H

#include "HartreeFock/HFOperator.h"

/** Add an extra exchange potential to a HF operator.
    Extra exchange potential is stored in currentExchangePotential.
    Note sign is for normal electrostatic potential: V_nucleus(r) > 0.
    Subclasses must define CalculateExtraExchange().
 */
class ExchangeDecorator : public HFOperatorDecorator<HFBasicDecorator, ExchangeDecorator>
{
public:
    ExchangeDecorator(pHFOperator wrapped_hf, pIntegrator integration_strategy = nullptr):
        BaseDecorator(wrapped_hf, integration_strategy), scale(1.), currentExchangePotential(-1)
    {}

    void SetScale(double factor) { scale = factor; }
    double GetScale() const { return scale; }

public:
    virtual void Alert() override;

    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(const Orbital& approximation) override;
    
    /** Get exchange (nonlocal) potential. */
    virtual SpinorFunction GetExchange(pOrbitalConst approximation) const override;

    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const override;
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const override;
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const override;

public:
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;

protected:
    /** Overwrite this in subclasses. */
    virtual SpinorFunction CalculateExtraExchange(const SpinorFunction& s) const
    {   return SpinorFunction(s.Kappa());
    };

protected:
    double scale;
    SpinorFunction currentExchangePotential;
};

#endif
