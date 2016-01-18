#ifndef HARTREE_Y_FUNCTION_H
#define HARTREE_Y_FUNCTION_H

#include "SpinorOperator.h"
#include "CoulombOperator.h"

/** HartreeY is a radial one-body operator defined by
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
 
    Each HartreeY operator needs to know whether it is "reversible" in the sense
    \f$ Y^k_{cd} = Y^k_{dc} \f$, i.e. SetParameters(k, c, d) == SetParameters(k, d, c).
    The boolean two_body_reverse_symmetry_exists contains this information, which is
    accessible via the function ReverseSymmetryExists().
 */
class HartreeYBase
{
public:
    HartreeYBase(pIntegrator integration_strategy = nullptr): K(-1), integrator(integration_strategy), two_body_reverse_symmetry_exists(true) {}
    virtual ~HartreeYBase() {}

    /** Set k and SpinorFunctions c and d according to definition \f$ Y^k_{cd} \f$.
        PRE: new_c and new_d must point to valid objects.
        Return false if resulting HartreeY function is zero (i.e. isZero() returns true).
     */
    virtual bool SetParameters(int new_K, pSpinorFunctionConst new_c, pSpinorFunctionConst new_d)
    {   K = new_K;
        return false;
    }

    /** Check whether the potential Y^k_{cd} is zero. (e.g. angular momentum conditions not satisfied). */
    virtual bool isZero() const { return true; }

    /** Whether this is "reversible" in the sense \f$ Y^k_{cd} = Y^k_{dc} \f$,
        i.e. SetParameters(k, c, d) == SetParameters(k, d, c).
     */
    virtual bool ReverseSymmetryExists() const { return two_body_reverse_symmetry_exists; }

    /** Call SetParameters(k, c, d) with k set to its smallest value with non-zero potential.
        PRE: new_c and new_d must point to valid objects.
        Return resulting k, or -1 if there is no such k.
     */
    virtual int SetOrbitals(pSpinorFunctionConst new_c, pSpinorFunctionConst new_d) { return -1; }

    /** Call SetParameters(k, c, d) with previously set orbitals.
        Return false if resulting HartreeY function is zero (i.e. isZero() returns true).
     */
    virtual bool SetK(int new_K)
    {   K = new_K;
        return false;
    }

    /** Using previously set orbitals, set k to its smallest value with non-zero potential.
        Return resulting k, or -1 if there is no such k.
     */
    virtual int SetMinK() { return -1; }

    /** Increment k to next value with non-zero potential using current orbitals.
        Return resulting k, or -1 if there is no such k.
     */
    virtual int NextK() { return false; }

    /** Return largest possible value of k giving a non-zero potential with current orbitals, or -1 if there is no such k. */
    virtual int GetMaxK() const { return -1; }

    /** Get multipolarity K for this operator. */
    virtual int GetK() const
    {   return K;
    }

    /** Get parity of this operator. */
    virtual Parity GetParity() const
    {   return (K%2? Parity::odd: Parity::even);
    }

    virtual pIntegrator GetIntegrator() const
    {   return integrator;
    }

    /** Deep copy of the HartreeY object, particularly including wrapped objects.
        The caller must take responsibility for deallocating the clone. A typical idiom would be
        to wrap the pointer immediately with a shared pointer, e.g.:
            pHartreeY cloned(old_hartreeY.Clone())
     */
    virtual HartreeYBase* Clone() const { return new HartreeYBase(integrator); }

    /** < b | t | a > for an operator t. */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a) const
    {   return GetMatrixElement(b, a, false);
    }

    virtual double GetMatrixElement(const Orbital& b, const Orbital& a, bool reverse) const
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
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const { return ApplyTo(a, a.Kappa(), false); }
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, bool reverse) const { return ApplyTo(a, a.Kappa(), reverse); }

    /** Potential = t | a > for an operator t such that the resulting Potential.Kappa() == kappa_b.
        i.e. t | a > has kappa == kappa_b.
     */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const { return ApplyTo(a, kappa_b, false); }
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse) const
    {   return SpinorFunction(kappa_b);
    }

protected:
    pIntegrator integrator;
    int K;
    bool two_body_reverse_symmetry_exists;
    pSpinorFunctionConst c;
    pSpinorFunctionConst d;
};

typedef std::shared_ptr<HartreeYBase> pHartreeY;
typedef std::shared_ptr<const HartreeYBase> pHartreeYConst;

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
class HartreeY : public HartreeYBase, public LatticeObserver
{
public:
    HartreeY(pIntegrator integration_strategy, pCoulombOperator coulomb);
    virtual ~HartreeY() {}

    /** Resize potential to match lattice. */
    virtual void Alert() override;

    /** Set k and SpinorFunctions c and d according to definition \f$ Y^k_{cd} \f$.
        PRE: new_c and new_d must point to valid objects.
        Return false if resulting HartreeY function is zero (i.e. isZero() returns true).
     */
    virtual bool SetParameters(int new_K, pSpinorFunctionConst new_c, pSpinorFunctionConst new_d) override;

    /** Check whether the potential Y^k_{cd} is zero. (e.g. angular momentum conditions not satisfied). */
    virtual bool isZero() const override { return (potential.size() == 0); }

    /** Call SetParameters(k, c, d) with k set to its smallest value with non-zero potential.
        PRE: new_c and new_d must point to valid objects.
        Return resulting k, or -1 if there is no such k.
     */
    virtual int SetOrbitals(pSpinorFunctionConst new_c, pSpinorFunctionConst new_d) override;

    /** Call SetParameters(k, c, d) with previously set orbitals.
        Return false if resulting HartreeY function is zero (i.e. isZero() returns true).
     */
    virtual bool SetK(int new_K) override
    {   return SetParameters(new_K, c, d);
    }

    /** Using previously set orbitals, set k to its smallest value with non-zero potential.
        Return resulting k, or -1 if there is no such k.
     */
    virtual int SetMinK() override { return SetOrbitals(c, d); }

    /** Increment k to next value with non-zero potential using current orbitals.
        Return resulting k, or -1 if there is no such k.
     */
    virtual int NextK() override;

    /** Return largest possible value of k giving a non-zero potential with current orbitals, or -1 if there is no such k. */
    virtual int GetMaxK() const override;

    /** Deep copy of this HartreeY object.
        NB: uses the same integrator and coulomb operator.
     */
    virtual HartreeY* Clone() const override;

    /** < b | t | a > for an operator t. */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a, bool reverse) const override;

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
    HartreeYDecorator(pHartreeY wrapped, pIntegrator integration_strategy = nullptr):
        HartreeYBase(integration_strategy)
    {   component = wrapped;
        if(!integrator)
            integrator = component->GetIntegrator();
    }

    virtual bool SetParameters(int new_K, pSpinorFunctionConst new_c, pSpinorFunctionConst new_d) override
    {   bool comp_ret = component->SetParameters(new_K, new_c, new_d);
        bool ret = SetLocalParameters(new_K, new_c, new_d);
        return comp_ret || ret;
    }

    /** Should be overwritten to take boolean AND of all wrapped objects. */
    virtual bool isZero() const override
    {   return component->isZero();
    }

    /** Whether this is "reversible" in the sense \f$ Y^k_{cd} = Y^k_{dc} \f$,
        i.e. SetParameters(k, c, d) == SetParameters(k, d, c).
     */
    virtual bool ReverseSymmetryExists() const override final
    {   return two_body_reverse_symmetry_exists && component->ReverseSymmetryExists();
    }

    virtual int SetOrbitals(pSpinorFunctionConst new_c, pSpinorFunctionConst new_d) override
    {   int new_K = component->SetOrbitals(new_c, new_d);
        SetLocalParameters(new_K, new_c, new_d);
        return K;
    }

    virtual bool SetK(int new_K) override
    {   bool ret = SetLocalParameters(K, c, d);
        return component->SetK(new_K) || ret;
    }

    virtual int SetMinK() override { return SetOrbitals(c, d); }

    virtual int NextK() override
    {   int ret = component->NextK();
        if(ret != -1)
            SetLocalParameters(component->GetK(), c, d);
        else
            K = -1;

        return ret;
    }

    virtual int GetMaxK() const override
    {   return component->GetMaxK();
    }
    
    /** Deep copy of the HartreeY object, including wrapped objects. */
    virtual HartreeYDecorator* Clone() const override
    {   pHartreeY wrapped_clone(component->Clone());
        return new HartreeYDecorator(wrapped_clone);
    }

    /** < b | t | a > for an operator t. */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a, bool reverse) const override
    {   return component->GetMatrixElement(b, a, reverse);
    }

    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse) const override
    {   return component->ApplyTo(a, kappa_b, reverse);
    }

protected:
    /** SetParameters for this object only. */
    virtual bool SetLocalParameters(int new_K, pSpinorFunctionConst new_c, pSpinorFunctionConst new_d)
    {   K = new_K; c = new_c; d = new_d;
        return false;
    }

    pHartreeY component;
};

#endif
