#ifndef HARTREE_Y_FUNCTION_H
#define HARTREE_Y_FUNCTION_H

#include "SpinorOperator.h"
#include "CoulombOperator.h"

/** HartreeY is a one-body operator defined by
    \f[
        Y^k_{cd}(r) = \int \frac{r_<^k}{r_>^{k+1}} \psi_c^\dag (r') \psi_d (r') dr'
                      . \xi(k + c.L() + d.L()) . \Delta(k, c.J(), d.J())
    \f]
    Operated on with functions \f$ \psi_a \f$ and \f$ \psi_b \f$ it yields the Coulomb radial matrix element
    \f[
        R^k (ac, bd) = \left< a \right| Y^k_{cd} \left| b \right>
    \f]
    This operator is also extensible via wrapped decorator objects.

    Sometimes we may want to store a lot of these operators, so we want them to be lightweight.
    Therefore we provide a base interface (HartreeYBase) from which decorators and HartreeY itself inherit.
    This allows us to remove potentials, etc. from decorator objects.
    HartreeYBase itself is the zero operator \f$ Y^k_{cd} = 0 \f$,
    used to wrap Decorators around when you want just the extra bit.
 */
class HartreeYBase : public SpinorOperator
{
public:
    HartreeYBase(pOPIntegrator integration_strategy = nullptr): SpinorOperator(integration_strategy) {}
    virtual ~HartreeYBase() {}

    /** Set k and SpinorFunctions c and d according to definition \f$ Y^k_{cd} \f$.
        Return false if resulting HartreeY function is zero (i.e. isZero() returns true).
     */
    virtual bool SetParameters(int new_K, const SpinorFunction& c, const SpinorFunction& d) { return false; }

    /** Check whether the potential Y^k_{cd} is zero. (e.g. angular momentum conditions not satisfied). */
    virtual bool isZero() const { return true; }

    /** < b | t | a > for an operator t. */
    virtual double GetMatrixElement(const SpinorFunction& b, const SpinorFunction& a) const override
    {   return GetMatrixElement(b, a, false);
    }

    virtual double GetMatrixElement(const SpinorFunction& b, const SpinorFunction& a, bool reverse) const
    {   SpinorFunction ta = ApplyTo(a, b.Kappa(), reverse);
        if(ta.size())
        {   if(!integrator)
                throw "HartreeYBase::GetMatrixElement(): no integrator found.";
            return integrator->GetInnerProduct(ta, b);
        }
        else
            return 0.0;
    }

    /** Potential = t | a > for an operator t such that the resulting Potential has the same angular symmetry as a.
     i.e. t | a > has kappa == kappa_a.
     */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override { return ApplyTo(a, a.Kappa(), false); }
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, bool reverse) const { return ApplyTo(a, a.Kappa(), reverse); }

    /** Potential = t | a > for an operator t such that the resulting Potential.Kappa() == kappa_b.
     i.e. t | a > has kappa == kappa_b.
     */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const { return ApplyTo(a, kappa_b, false); }
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse = false) const
    {   return SpinorFunction(kappa_b);
    }
};

typedef boost::shared_ptr<HartreeYBase> pHartreeY;
typedef boost::shared_ptr<const HartreeYBase> pHartreeYConst;

/** HartreeY is the one-body operator defined by
    \f[
        Y^k_{cd}(r) = \int \frac{r_<^k}{r_>^{k+1}} \psi_c^\dag (r') \psi_d (r') dr'
                      . \xi(k + c.L() + d.L()) . \Delta(k, c.J(), d.J())
    \f]
    Operated on with functions \f$ \psi_a \f$ and \f$ \psi_b \f$ it yields the Coulomb radial matrix element
    \f[
        R^k (ac, bd) = \left< a \right| Y^k_{cd} \left| b \right>
    \f]
    This operator is also extensible via wrapped decorator objects.

    A huge saving can often be made when reversing the spinors c and d,
    e.g. HartreeY itself obeys \f$ Y^k_{cd} = Y^k_{dc} \f$.
    Therefore a boolean argument to the usual one-body operator functions provide for reversed versions.
    (HartreeY itself ignores the reverse boolean.)
 */
class HartreeY : public HartreeYBase
{
public:
    HartreeY(pOPIntegrator integration_strategy, pCoulombOperator coulomb);
    virtual ~HartreeY() {}

    /** Set k and SpinorFunctions c and d according to definition \f$ Y^k_{cd} \f$.
        Return false if resulting HartreeY function is zero (i.e. isZero() returns true).
     */
    virtual bool SetParameters(int new_K, const SpinorFunction& c, const SpinorFunction& d) override;

    /** Check whether the potential Y^k_{cd} is zero. (e.g. angular momentum conditions not satisfied). */
    virtual bool isZero() const override { return (potential.size() == 0); }

    /** < b | t | a > for an operator t. */
    virtual double GetMatrixElement(const SpinorFunction& b, const SpinorFunction& a, bool reverse) const override;

    /** Potential = t | a > for an operator t such that the resulting Potential.Kappa() == kappa_b.
        i.e. t | a > has kappa == kappa_b.
     */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse = false) const override;

protected:
    pCoulombOperator coulomb;
    RadialFunction potential;
};

class HartreeYDecorator : public HartreeYBase
{
public:
    HartreeYDecorator(pHartreeY wrapped, pOPIntegrator integration_strategy = nullptr):
        HartreeYBase(integration_strategy)
    {   component = wrapped;
        if(!integrator)
            integrator = component->GetOPIntegrator();
    }

    virtual bool SetParameters(int new_K, const SpinorFunction& c, const SpinorFunction& d) override
    {   return component->SetParameters(new_K, c, d);
    }

    /** Should be overwritten to take boolean AND of all wrapped objects. */
    virtual bool isZero() const override
    {   return component->isZero();
    }

    /** Get maximum multipolarity K for this operator and all wrapped (added) operators. */
    virtual int GetMaxK() const override { return mmax(K, component->GetMaxK()); }

    /** < b | t | a > for an operator t. */
    virtual double GetMatrixElement(const SpinorFunction& b, const SpinorFunction& a, bool reverse) const override
    {   return component->GetMatrixElement(b, a, reverse);
    }

    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse) const override
    {   return component->ApplyTo(a, kappa_b, reverse);
    }

protected:
    pHartreeY component;
};

#endif
