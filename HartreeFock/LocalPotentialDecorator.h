#ifndef LOCAL_POTENTIAL_DECORATOR
#define LOCAL_POTENTIAL_DECORATOR

#include "SpinorODE.h"
#include "Operator.h"

/** Generic (abstract-ish) class to add an extra local potential to a HF operator.
    Local member extraLocalPotential must be set by subclasses via SetCore() and/or SetODEParameters().
 */
class LocalPotentialDecorator : public OneBodyOperatorDecorator, public SpinorODEDecorator
{
public:
    LocalPotentialDecorator(OneBodyOperator* wrapped_OBO, SpinorODE* wrapped_ODE, OPIntegrator* integration_strategy = NULL);

public:
    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const;
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const;
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const;

public:
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const;

protected:
    RadialFunction extraLocalPotential;
};

/** Decorate HF operator with a local approximation to the exchange potential. Good for first approximations. */
class LocalExchangeApproximation : public LocalPotentialDecorator
{
public:
    LocalExchangeApproximation(OneBodyOperator* wrapped_OBO, SpinorODE* wrapped_ODE);

    virtual void SetCore(const Core* hf_core);
};

#endif
