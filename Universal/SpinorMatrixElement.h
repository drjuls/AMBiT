#ifndef SPINOR_MATRIX_ELEMENT_H
#define SPINOR_MATRIX_ELEMENT_H

#include "Universal/SpinorFunction.h"
#include "HartreeFock/Orbital.h"
#include "Include.h"

/** SpinorMatrixElement is a base class for calculating radial matrix elements with an orbital basis.
    It follows the Decorator (Wrapper) pattern, so it is recursively extensive.
    SpinorMatrixElement is itself a valid zero operator, < b | 0 | a > = 0,
    useful to wrap Decorators around when you want just the extra bit.
 */
class SpinorMatrixElement : public std::enable_shared_from_this<SpinorMatrixElement>
{
public:
    SpinorMatrixElement(int K = 0, pOPIntegrator integration_strategy = nullptr): integrator(integration_strategy), K(K) {}

    /** Reduced matrix element < b || t(K) || a > for our operator t(K).
        Usually the operation t|a> makes sense for any Orbital, but in order to have a
        matrix element one of the wavefunctions must be bounded, and hence an "Orbital".
     */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a) const
    {   return 0.;
    }

    /** Get multipolarity K for this operator and all wrapped (added) operators. */
    virtual int GetK() const
    {   return K;
    }

    virtual pOPIntegrator GetOPIntegrator() const
    {   return integrator;
    }

protected:
    pOPIntegrator integrator;
    int K;
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
    SpinorMatrixElementDecorator(int K, pSpinorMatrixElement wrapped, pOPIntegrator integration_strategy = nullptr): SpinorMatrixElement(K, integration_strategy)
    {   component = wrapped;
        if(!integrator)
            integrator = component->GetOPIntegrator();
    }

    /** Reduced matrix element < b || t(K) || a > for our operator t(K). */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a) const override
    {   return component->GetMatrixElement(b, a);
    }

protected:
    pSpinorMatrixElement component;
};

#endif
