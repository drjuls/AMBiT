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
    
    /** Potential = t | a > for an operator t. */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const = 0;

    virtual pOPIntegrator GetOPIntegrator() const
    {   return integrator;
    }

protected:
    pOPIntegrator integrator;
};

/** OneBodyOperatorDecorator is for adding extra terms to an existing operator.
    The Decorator pattern allows nesting of additional terms in any order.
    When using, remember that the Decorator wraps objects, not classes.
    Wrap around a zero operator object if only the extra term is required.
 */
class OneBodyOperatorDecorator : public OneBodyOperator
{
public:
    OneBodyOperatorDecorator(OneBodyOperator* wrapped, pOPIntegrator integration_strategy = pOPIntegrator()): OneBodyOperator(integration_strategy)
    {   component = wrapped;
        if(!integrator)
            integrator = component->GetOPIntegrator();
    }

    /** Subclasses should call this inherited function. */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const
    {   return component->ApplyTo(a);
    }

protected:
    OneBodyOperator* component;
};

/** < b | 0 | a > = 0.
    Used to wrap Decorators around when you want just the extra bit.
 */
class ZeroOperator : public OneBodyOperator
{
public:
    ZeroOperator(): OneBodyOperator() {}

public:
    double GetMatrixElement(const SpinorFunction& a, const SpinorFunction& b) const;
    SpinorFunction ApplyTo(const SpinorFunction& a) const;
};

#endif
