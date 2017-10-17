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
class HartreeYBase : public std::enable_shared_from_this<HartreeYBase>
{
    friend class HartreeYBasicDecorator;
public:
    HartreeYBase(pIntegrator integration_strategy = nullptr):
        K(-1), integrator(integration_strategy), two_body_reverse_symmetry_exists(true), off_parity_exists(false), parent(nullptr), lightweight_mode(false)
    {}
    virtual ~HartreeYBase() {}

    /** Whether this is "reversible" in the sense \f$ Y^k_{cd} = Y^k_{dc} \f$,
        i.e. SetParameters(k, c, d) == SetParameters(k, d, c).
     */
    virtual bool ReverseSymmetryExists() const { return two_body_reverse_symmetry_exists; }

    /** Whether any off-parity radial integrals are non-zero, i.e. \f$ (l_1 + l_3 + k) \f$ is odd.
        The total parity must still hold: \f$(l_1 + l_2 + l_3 + l_4)\f$ is even.
     */
    virtual bool OffParityExists() const { return off_parity_exists; }

    /** Set lightweight_mode, used for getting K values and --check-sizes.
        Warning: It prevents long calculations, but if set then GetMatrixElement() and ApplyTo() will not function.
     */
    virtual void SetLightWeightMode(bool setting) { lightweight_mode = setting; }
    bool GetLightWeightMode() const { return lightweight_mode; }

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

    /** Return smallest possible value of k giving a non-zero potential with current orbitals, or -1 if there is no such k. */
    virtual int GetMinK() const { return -1; }

    /** Return largest possible value of k giving a non-zero potential with current orbitals, or -1 if there is no such k. */
    virtual int GetMaxK() const { return -1; };

    /** Increment k to next value with non-zero potential using current orbitals.
        Return resulting k, or -1 if there is no such k.
     */
    virtual int NextK() { return -1; }

    /** Get multipolarity K for this operator. */
    int GetK() const
    {   return K;
    }

    /** Get parity of this operator. */
    virtual Parity GetParity() const
    {   return (K%2? Parity::odd: Parity::even);
    }

    pIntegrator GetIntegrator() const
    {   return integrator;
    }

    /** Deep copy of the HartreeY object, particularly including wrapped objects. */
    virtual std::shared_ptr<HartreeYBase> Clone() const { return std::make_shared<HartreeYBase>(*this); }

    /** < b | t | a > for an operator t. */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a, bool reverse = false) const
    {   SpinorFunction ta = ApplyTo(a, b.Kappa(), reverse);
        if(ta.size())
        {   if(!integrator)
            {   *errstream << "HartreeYBase::GetMatrixElement(): no integrator found." << std::endl;
                exit(1);
            }
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
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse = false) const
    {   return SpinorFunction(kappa_b);
    }

    void SetParent(std::shared_ptr<HartreeYBase> new_parent) { parent = new_parent.get(); }

protected:
    pIntegrator integrator;
    bool two_body_reverse_symmetry_exists;
    bool off_parity_exists;
    HartreeYBase* parent;   //!< Pointer to parent decorator, or null if this is the top-level HartreeY.

    int K;
    pSpinorFunctionConst c;
    pSpinorFunctionConst d;

    bool lightweight_mode;  //!< Prevents long calculations, only used for getting K
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

    /** Return smallest possible value of k giving a non-zero potential with current orbitals, or -1 if there is no such k. */
    virtual int GetMinK() const override;

    /** Return largest possible value of k giving a non-zero potential with current orbitals, or -1 if there is no such k. */
    virtual int GetMaxK() const override;

    /** Increment k to next value with non-zero potential using current orbitals.
        Return resulting k, or -1 if there is no such k.
     */
    virtual int NextK() override;

    /** Deep copy of this HartreeY object. */
    virtual pHartreeY Clone() const override;

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

/** HartreeYDecorators add additional terms to \f$ Y^k_{cd} \f$.
    Derive classes from it using the HartreeYDecorator template below:
        class MyDecorator : public HartreeYDecorator<HartreeYBasicDecorator, MyDecorator>

    Decorator writers should set two_body_reverse_symmetry_exists to the appropriate value for that operator.
    Functions that need to be overridden:
        SetLocalParameters(k, c, d)
        isZeroLocal()    - for current K, c, d
        GetLocalMinK()   - for current c, d
        GetLocalMaxK()   - for current c, d
 */
class HartreeYBasicDecorator : public HartreeYBase
{
public:
    HartreeYBasicDecorator(pHartreeY wrapped, pIntegrator integration_strategy = nullptr):
        HartreeYBase(integration_strategy)
    {   component = wrapped;
        component->parent = this;
        if(!integrator)
            integrator = component->GetIntegrator();
    }

    virtual ~HartreeYBasicDecorator()
    {   component->parent = nullptr;
    }

    virtual void SetLightWeightMode(bool setting) override final
    {   lightweight_mode = setting;
        component->SetLightWeightMode(setting);
    }

    virtual bool SetParameters(int new_K, pSpinorFunctionConst new_c, pSpinorFunctionConst new_d) override
    {
        bool comp_ret = component->SetParameters(new_K, new_c, new_d);
        bool ret = SetLocalParameters(new_K, new_c, new_d);
        return comp_ret || ret;
    }

    /** Should be overwritten to take boolean AND of all wrapped objects. */
    virtual bool isZero() const override final
    {   return isZeroLocal() && component->isZero();
    }

    /** Return smallest possible value of k giving a non-zero potential with current orbitals, or -1 if there is no such k. */
    virtual int GetMinK() const override final
    {
        if(c && d)
        {   int Kmin = GetLocalMinK();
            int otherKmin = component->GetMinK();
            if(otherKmin == -1)
                return Kmin;
            else
                return mmin(Kmin, otherKmin);
        }
        else
            return component->GetMinK();
    }

    /** Return largest possible value of k giving a non-zero potential with current orbitals, or -1 if there is no such k. */
    virtual int GetMaxK() const override final
    {
        if(c && d)
        {   int Kmax = GetLocalMaxK();
            int otherKmax = component->GetMaxK();
            if(otherKmax == -1)
                return Kmax;
            else
                return mmax(Kmax, otherKmax);
        }
        else
            return component->GetMaxK();
    }

    /** Whether this is "reversible" in the sense \f$ Y^k_{cd} = Y^k_{dc} \f$,
        i.e. SetParameters(k, c, d) == SetParameters(k, d, c).
     */
    virtual bool ReverseSymmetryExists() const override final
    {   return two_body_reverse_symmetry_exists && component->ReverseSymmetryExists();
    }

    virtual bool OffParityExists() const override final
    {   return off_parity_exists || component->OffParityExists();
    }

    virtual int SetOrbitals(pSpinorFunctionConst new_c, pSpinorFunctionConst new_d) override final
    {
        c = new_c;
        d = new_d;
        int myKmin = GetLocalMinK();
        int belowKmin = component->SetOrbitals(new_c, new_d);

        if(myKmin == -1)
            K = belowKmin;
        else if(belowKmin == -1)
            K = myKmin;
        else
            K = mmin(myKmin, belowKmin);

        if((K != -1) && (parent == nullptr))
            SetK(K);

        return K;
    }

    virtual bool SetK(int new_K) override
    {   bool ret = SetLocalParameters(K, c, d);
        return component->SetK(new_K) || ret;
    }

    virtual int NextK() override final
    {
        // Only the top-level decorator does this.
        if((K != -1) && (parent == nullptr))
        {
            int maxK = GetMaxK();

            bool nonzero = true;
            do
            {   K++;
                if(K > maxK)
                {   K = -1;
                }
                else
                {   nonzero = SetK(K);
                }
            }while(K != -1 && !nonzero);
        }

        return K;
    }

    /** Deep copy of the HartreeY object, including wrapped objects. */
    virtual pHartreeY Clone() const override
    {
        std::shared_ptr<HartreeYBasicDecorator> ret(std::make_shared<HartreeYBasicDecorator>(*this));
        ret->component = component->Clone();
        ret->component->SetParent(ret);
        return ret;
    }

    /** < b | t | a > for an operator t. */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a, bool reverse = false) const override
    {   return component->GetMatrixElement(b, a, reverse);
    }

    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse = false) const override
    {   return component->ApplyTo(a, kappa_b, reverse);
    }

protected:
    /** SetParameters for this object only. */
    virtual bool SetLocalParameters(int new_K, pSpinorFunctionConst new_c, pSpinorFunctionConst new_d)
    {   K = new_K; c = new_c; d = new_d;
        return false;
    }

    /** Should be overwritten to take boolean AND of all wrapped objects. */
    virtual bool isZeroLocal() const
    {   return true;
    }

    /** Return smallest possible value of k giving a non-zero potential with current orbitals, or -1 if there is no such k. */
    virtual int GetLocalMinK() const
    {   return -1;
    }

    /** Return largest possible value of k giving a non-zero potential with current orbitals, or -1 if there is no such k. */
    virtual int GetLocalMaxK() const
    {   return -1;
    }

    pHartreeY component;
};

/** HartreeYDecorator is a base class for HartreeY decorators.
    It follows the curiously recurring template pattern to provide Clone() function.
    When deriving from it use the syntax
        class MyDecorator : public HartreeYDecorator<BaseDecorator, MyDecorator>
 */
template <typename Base, typename Derived>
class HartreeYDecorator : public Base
{
public:
    using Base::Base;

    HartreeYDecorator(const Derived& other):
        Base(other)
    {}

    virtual pHartreeY Clone() const override
    {
        std::shared_ptr<Derived> ret = std::make_shared<Derived>(static_cast<Derived const &>(*this));
        ret->component = this->component->Clone();
        ret->component->SetParent(ret);
        return ret;
    }

protected:
    typedef HartreeYDecorator<Base, Derived> BaseDecorator;
};

#endif
