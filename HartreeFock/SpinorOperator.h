#ifndef ONE_BODY_OPERATOR_H
#define ONE_BODY_OPERATOR_H

#include "Orbital.h"
#include "OpIntegrator.h"
#include "Universal/SpinorMatrixElement.h"

/** SpinorOperator is a base class for calculating radial matrix elements with an orbital basis.
    It follows the Decorator (Wrapper) pattern, so it is recursively extensive.
    Operators are also components of the Strategy pattern; they are initialised with an Integrator that the client can choose.
    Subclasses may choose their default integrator.
    SpinorOperator is itself a valid zero operator, < b | 0 | a > = 0,
    useful to wrap Decorators around when you want just the extra bit.
 */
class SpinorOperator : public SpinorMatrixElement
{
public:
    SpinorOperator(pOPIntegrator integration_strategy = nullptr): SpinorMatrixElement(0, integration_strategy) {}

    /** < b | t | a > for an operator t.
        Default behaviour: take ApplyTo(a, b.Kappa) and integrate with b.
     */
    virtual double GetMatrixElement(const Orbital& b, const SingleParticleWavefunction& a) const override;

    /** Potential = t | a > for an operator t such that the resulting Potential has the same angular symmetry as a.
        i.e. t | a > has kappa == kappa_a.
     */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const
    {   return ApplyTo(a, a.Kappa());
    }

    /** Potential = t | a > for an operator t such that the resulting Potential.Kappa() == kappa_b.
        i.e. t | a > has kappa == kappa_b.
     */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const
    {   return SpinorFunction(kappa_b);
    }
};

typedef std::shared_ptr<SpinorOperator> pSpinorOperator;
typedef std::shared_ptr<const SpinorOperator> pSpinorOperatorConst;

/** SpinorOperatorDecorator is for adding extra terms to an existing operator.
    The Decorator pattern allows nesting of additional terms in any order.
    When using, remember that the Decorator wraps objects, not classes.
    Wrap around a zero operator object if only the extra term is required.
 */
class SpinorOperatorDecorator : public SpinorOperator
{
public:
    SpinorOperatorDecorator(pSpinorOperator wrapped, pOPIntegrator integration_strategy = nullptr): SpinorOperator(integration_strategy)
    {   component = wrapped;
        if(!integrator)
            integrator = component->GetOPIntegrator();
    }

    /** Subclasses may call this inherited function. */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override
    {   return ApplyTo(a, a.Kappa());
    }

    /** Subclasses should call this inherited function. */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const override
    {   return component->ApplyTo(a, kappa_b);
    }

    /** Get maximum multipolarity K for this operator and all wrapped (added) operators. */
    virtual int GetMaxK() const override { return mmax(K, component->GetMaxK()); }

protected:
    pSpinorOperator component;
};

inline double SpinorOperator::GetMatrixElement(const Orbital& b, const SingleParticleWavefunction& a) const
{
    SpinorFunction ta = this->ApplyTo(a, b.Kappa());
    if(ta.size())
    {   if(!integrator)
        throw "SpinorOperator::GetMatrixElement(): no integrator found.";
        return integrator->GetInnerProduct(ta, b);
    }
    else
        return 0.0;
}

#endif
