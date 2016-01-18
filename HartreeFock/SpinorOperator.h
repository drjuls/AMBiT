#ifndef ONE_BODY_OPERATOR_H
#define ONE_BODY_OPERATOR_H

#include "Orbital.h"
#include "Integrator.h"
#include "Universal/SpinorMatrixElement.h"
#include "Universal/MathConstant.h"

/** SpinorOperator is a base class for calculating radial matrix elements with an orbital basis.
    Operators are also components of the Strategy pattern; they are initialised with an Integrator that the client can choose.
    Subclasses may choose their default integrator.
    At minimum, derived classes should implement ReducedApplyTo(const SpinorFunction& a, int kappa_b);
    the base class defines ApplyTo(), GetMatrixElement(), and GetReducedMatrixElement() from this,
    although of course these can be overridden.
 */
class SpinorOperator : public SpinorMatrixElement
{
public:
    using SpinorMatrixElement::SpinorMatrixElement;

    /** < b | t | a > for stretched states a and b.
        Default behaviour: take ApplyTo(a, b.Kappa) and integrate with b.
     */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a) const override;

    /** < b || t || a > for an operator t.
        Default behaviour: take ApplyTo(a, b.Kappa) and integrate with b.
     */
    virtual double GetReducedMatrixElement(const Orbital& b, const Orbital& a) const override;

    /** Potential = t | a > for an operator t such that the resulting Potential has the same angular symmetry as a.
        i.e. t | a > has kappa == kappa_a.
        Assumes stretched state: a.M() = a.J().
     */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const
    {   return ApplyTo(a, a.Kappa());
    }

    /** Potential = t | a > for an operator t such that the resulting Potential.Kappa() == kappa_b.
        i.e. t | a > has kappa == kappa_b.
        Assumes stretched states for a and b.
     */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const;

    /** Potential = t || a > for an operator t such that the resulting Potential has the same angular symmetry as a.
        i.e. t || a > has kappa == kappa_a.
     */
    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a) const
    {   return ReducedApplyTo(a, a.Kappa());
    }

    /** Take reduced matrix element of angular part and apply radial part of t to |a>.
        Potential = < Kappa_b || t^(ang) || Kappa_a > . t^(rad) | a >
        The resulting Potential.Kappa() == kappa_b.
     */
    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a, int kappa_b) const
    {   return SpinorFunction(kappa_b);
    }
};

typedef std::shared_ptr<SpinorOperator> pSpinorOperator;
typedef std::shared_ptr<const SpinorOperator> pSpinorOperatorConst;

inline double SpinorOperator::GetMatrixElement(const Orbital& b, const Orbital& a) const
{
    SpinorFunction ta = this->ReducedApplyTo(a, b.Kappa());
    if(!integrator)
    {   *errstream << "SpinorOperator::GetMatrixElement(): no integrator found." << std::endl;
        exit(1);
    }

    double stretched = MathConstant::Instance()->Electron3j(a.TwoJ(), b.TwoJ(), K, a.TwoJ(), -b.TwoJ());
    return integrator->GetInnerProduct(ta, b) * stretched;
}

inline double SpinorOperator::GetReducedMatrixElement(const Orbital& b, const Orbital& a) const
{
    SpinorFunction ta = this->ReducedApplyTo(a, b.Kappa());
    if(!integrator)
    {   *errstream << "SpinorOperator::GetMatrixElement(): no integrator found." << std::endl;
        exit(1);
    }

    return integrator->GetInnerProduct(ta, b);
}

inline SpinorFunction SpinorOperator::ApplyTo(const SpinorFunction& a, int kappa_b) const
{
    int b_TwoJ = 2 * abs(kappa_b) - 1;
    double stretched = MathConstant::Instance()->Electron3j(a.TwoJ(), b_TwoJ, K, a.TwoJ(), -b_TwoJ);
    return ReducedApplyTo(a, kappa_b) * stretched;
}

#endif
