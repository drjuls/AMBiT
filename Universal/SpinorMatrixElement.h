#ifndef SPINOR_MATRIX_ELEMENT_H
#define SPINOR_MATRIX_ELEMENT_H

#include "Universal/SpinorFunction.h"
#include "HartreeFock/Orbital.h"
#include "Include.h"

/** SpinorMatrixElement is a base class for calculating radial matrix elements with an orbital basis.
    It follows the Decorator (Wrapper) pattern, so it is recursively extensive.
    SpinorMatrixElement is itself a valid zero operator, < b || 0 || a > = 0,
    useful to wrap Decorators around when you want just the extra bit.
 */
class SpinorMatrixElement : public std::enable_shared_from_this<SpinorMatrixElement>
{
public:
    SpinorMatrixElement(pIntegrator integration_strategy = nullptr): K(0), P(Parity::even), integrator(integration_strategy) {}
    SpinorMatrixElement(int K, Parity P, pIntegrator integration_strategy = nullptr): K(K), P(P), integrator(integration_strategy) {}
    /** Initialise with P = (-1)^K */
    SpinorMatrixElement(int K, pIntegrator integration_strategy = nullptr):
        integrator(integration_strategy), K(K)
    {   P = (K%2? Parity::odd: Parity::even);
    }

    /** Reduced matrix element < b || t(K) || a > for our operator t(K).
        Usually the operation t|a> makes sense for any Orbital, but in order to have a
        matrix element one of the wavefunctions must be bounded, and hence an "Orbital".
     */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a) const
    {   return 0.;
    }

    /** Get multipolarity K for this operator. */
    virtual int GetK() const
    {   return K;
    }

    /** Get parity of this operator. */
    virtual Parity GetParity() const
    {   return P;
    }

    virtual pIntegrator GetIntegrator() const
    {   return integrator;
    }

protected:
    pIntegrator integrator;
    int K;
    Parity P;
};

typedef std::shared_ptr<SpinorMatrixElement> pSpinorMatrixElement;
typedef std::shared_ptr<const SpinorMatrixElement> pSpinorMatrixElementConst;

/** SpinorMatrixElementDecorator is for adding extra terms to an existing operator.
    The Decorator pattern allows nesting of additional terms in any order.
    When using, remember that the Decorator wraps objects, not classes.
 */
class SpinorMatrixElementDecorator : public SpinorMatrixElement
{
public:
    SpinorMatrixElementDecorator(int K, pSpinorMatrixElement wrapped, pIntegrator integration_strategy = nullptr): SpinorMatrixElement(K, integration_strategy)
    {   component = wrapped;
        if(!integrator)
            integrator = component->GetIntegrator();
    }

    /** Reduced matrix element < b || t(K) || a > for our operator t(K). */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a) const override
    {   return component->GetMatrixElement(b, a);
    }

protected:
    pSpinorMatrixElement component;
};

#endif
