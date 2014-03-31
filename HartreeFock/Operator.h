#ifndef OPERATOR_H
#define OPERATOR_H

#include "Orbital.h"
#include "OpIntegrator.h"

/** OneBodyOperator is an abstract class for calculating radial matrix elements with an orbital basis.
    It follows the Decorator (Wrapper) pattern, so it is recursively extensive.
    Operators are also components of the Strategy pattern; they are initialised with an Integrator that the client can choose.
    Subclasses may choose their default integrator.
 */
class OneBodyOperator
{
public:
    OneBodyOperator(pOPIntegrator integration_strategy = pOPIntegrator()): integrator(integration_strategy) {}

    /** < b | t | a > for an operator t. */
    virtual double GetMatrixElement(const SpinorFunction& b, const SpinorFunction& a) const;
    
    /** Potential = t | a > for an operator t such that the resulting Potential has the same angular symmetry as a.
        i.e. t | a > has kappa == kappa_a.
     */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const = 0;

    /** Potential = t | a > for an operator t such that the resulting Potential.Kappa() == kappa_b.
        i.e. t | a > has kappa == kappa_b.
        For many operators the initial and final kappas must be the same.
        So we provide a default for this function that is a zero spinor unless a.Kappa() == kappa_b.
     */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const;

    virtual pOPIntegrator GetOPIntegrator() const
    {   return integrator;
    }

protected:
    pOPIntegrator integrator;
};

typedef boost::shared_ptr<OneBodyOperator> pOneBodyOperator;
typedef boost::shared_ptr<const OneBodyOperator> pOneBodyOperatorConst;

/** OneBodyOperatorDecorator is for adding extra terms to an existing operator.
    The Decorator pattern allows nesting of additional terms in any order.
    When using, remember that the Decorator wraps objects, not classes.
    Wrap around a zero operator object if only the extra term is required.
 */
class OneBodyOperatorDecorator : public OneBodyOperator
{
public:
    OneBodyOperatorDecorator(pOneBodyOperator wrapped, pOPIntegrator integration_strategy = pOPIntegrator()): OneBodyOperator(integration_strategy)
    {   component = wrapped;
        if(!integrator)
            integrator = component->GetOPIntegrator();
    }

    /** Subclasses should call this inherited function. */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override
    {   return component->ApplyTo(a);
    }

    /** Subclasses may call this inherited function. */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const override
    {   return component->ApplyTo(a, kappa_b);
    }

protected:
    pOneBodyOperator component;
};

/** < b | 0 | a > = 0.
    Used to wrap Decorators around when you want just the extra bit.
 */
class ZeroOperator : public OneBodyOperator
{
public:
    ZeroOperator(): OneBodyOperator() {}

public:
    virtual double GetMatrixElement(const SpinorFunction& a, const SpinorFunction& b) const override;
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const override;
};

#endif
