#ifndef HF_OPERATOR_H
#define HF_OPERATOR_H

#include "SpinorOperator.h"
#include "SpinorODE.h"
#include "Core.h"
#include "CoulombOperator.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/FornbergDifferentiator.h"

/** The relativistic Hartree-Fock (Dirac-Fock) operator:
    \f[
    h = \left( \begin{array}{cc}
                  -V             & (-d/dr + \kappa/r)/\alpha \\
        (d/dr + \kappa/r)/\alpha &     -2/\alpha^2 - V
        \end{array} \right)
    \f]
    where V is the electrostatic potential (V > 0) which leads to -V for electrons.
 */
class HFOperator : public SpinorOperator, public SpinorODE
{
public:
    HFOperator(double Z, pCoreConst hf_core, pPhysicalConstant physical_constant, pIntegrator integration_strategy, pCoulombOperator coulomb);
    virtual ~HFOperator();

protected:
    /** Protected version, for use by Decorators and Clone(), does not run SetCore(). */
    HFOperator(double Z, const HFOperator& other);

public:
    /** Set/reset the Hartree-Fock core, from which the potential is derived. */
    virtual void SetCore(pCoreConst hf_core);
    virtual pCoreConst GetCore() const;

    virtual RadialFunction GetDirectPotential() const;  //!< Get the direct potential.
    virtual double GetZ() const { return Z; }           //!< Get nuclear charge.
    virtual double GetCharge() const { return charge; } //!< Get ion charge.

    virtual pPhysicalConstant GetPhysicalConstant() const { return physicalConstant; } //!< Get physical constants.
    virtual pFornbergDifferentiator GetDifferentiator() const { return differentiator; }

    /** Deep copy of the HFOperator object, particularly including wrapped objects (but not the core or physical constants). */
    virtual std::shared_ptr<HFOperator> Clone() const { return std::shared_ptr<HFOperator>(new HFOperator(Z, *this)); }

public:
    /** Extend/reduce direct potential to match lattice size. */
    virtual void Alert() override;

    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(int kappa, double energy, const SpinorFunction* exchange = NULL) override;
    
    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(const Orbital& approximation) override;

    /** Get exchange (nonlocal) potential. */
    virtual SpinorFunction GetExchange(pOrbitalConst approximation = pOrbitalConst()) const override;

    /** Get df/dr = w[0] and dg/dr = w[1] given point r, (f, g).
        PRE: w should be an allocated 2 dimensional array;
             latticepoint < size().
     */
    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const override;

    /** Get numerical coefficients of the ODE at the point r, (f,g).
        PRE: w_f, w_g, and w_const should be allocated 2 dimensional arrays;
             latticepoint < size().
     */
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const override;

    /** Get Jacobian (dw[i]/df and dw[i]/dg), and dw[i]/dr at a point r, (f, g).
        PRE: jacobian should be an allocated 2x2 matrix;
             dwdr should be an allocated 2 dimensional array;
             latticepoint < size().
     */
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const override;

    /** Get approximation to eigenfunction for first numpoints near the origin. */
    virtual void EstimateOrbitalNearOrigin(unsigned int numpoints, SpinorFunction& s) const override;

    /** Get approximation to eigenfunction for last numpoints far from the origin.
        This routine can change the size of the orbital.
     */
    virtual void EstimateOrbitalNearInfinity(unsigned int numpoints, Orbital& s) const override;

public:
    /** I prefer to work with normal (not reduced) operators for Hartree-Fock, since in this case
            h | a > = E_a | a >
        For operators with K = 0, hence for non-zero operators j_a = j_b and m_a = m_b, these are related by
            < b || h || a > = sqrt(2j+1) < b | h | a >
     */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;

    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const override
    {
        if(a.Kappa() == kappa_b)
            return ApplyTo(a);
        else
            return SpinorFunction(kappa_b);
    }

    virtual double GetMatrixElement(const Orbital& b, const Orbital& a) const override
    {
        SpinorFunction ta = ApplyTo(a, b.Kappa());
        if(!integrator)
        {   *errstream << "HFOperator::GetMatrixElement(): no integrator found." << std::endl;
            exit(1);
        }
        return integrator->GetInnerProduct(ta, b);
    }

    /** Overriden so multiplication by sqrt(2j+1) is applied to matrix element, not entire orbital. */
    virtual double GetReducedMatrixElement(const Orbital& b, const Orbital& a) const override
    {
        SpinorFunction ta = this->ApplyTo(a, b.Kappa());
        if(!integrator)
        {   *errstream << "HFOperator::GetMatrixElement(): no integrator found." << std::endl;
            exit(1);
        }
        return integrator->GetInnerProduct(ta, b) * sqrt(a.TwoJ() + 1);
    }

    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a) const override
    {   return ApplyTo(a) * sqrt(a.TwoJ() + 1);
    }

    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a, int kappa_b) const override
    {   return ApplyTo(a, kappa_b) * sqrt(a.TwoJ() + 1);
    }

protected:
    virtual SpinorFunction CalculateExchange(const SpinorFunction& s) const;

protected:
    double Z;
    double charge;
    pCoreConst core;
    pCoulombOperator coulombSolver;
    pPhysicalConstant physicalConstant;
    pFornbergDifferentiator differentiator;

    RadialFunction directPotential;
    SpinorFunction currentExchangePotential;
    double currentEnergy;
    int currentKappa;
};

typedef std::shared_ptr<HFOperator> pHFOperator;
typedef std::shared_ptr<const HFOperator> pHFOperatorConst;

/** HFBasicDecorator is a simple base class for HFOperator decorators.
    Derive classes from it using the HFOperatorDecorator template below:
        class MyDecorator : public HFOperatorDecorator<HFBasicDecorator, MyDecorator>
 */
class HFBasicDecorator : public HFOperator
{
public:
    /** If integration_strategy is null, take from decorated_object. */
    HFBasicDecorator(pHFOperator decorated_object, pIntegrator integration_strategy = pIntegrator()):
        HFOperator(decorated_object->GetZ(), *decorated_object), wrapped(decorated_object)
    {
        if(integration_strategy != nullptr)
            integrator = integration_strategy;

        core = decorated_object->GetCore();
        charge = Z - double(core->NumElectrons());
        differentiator = wrapped->GetDifferentiator();
        IncludeExchange(wrapped->IncludeExchange());
    }
    HFBasicDecorator(const HFBasicDecorator& other):
        HFOperator(other.GetZ(), other), wrapped(other.wrapped)
    {}
    virtual ~HFBasicDecorator() {}

    /** Set/reset the Hartree-Fock core, from which the potential is derived. */
    virtual void SetCore(pCoreConst hf_core) override
    {   wrapped->SetCore(hf_core);
        core = hf_core;
        charge = Z - double(core->NumElectrons());
    }

    virtual RadialFunction GetDirectPotential() const override
    {   return wrapped->GetDirectPotential();
    }

    /** Deep copy of the HFOperator object, including wrapped objects. */
    virtual pHFOperator Clone() const override
    {
        std::shared_ptr<HFBasicDecorator> ret(std::make_shared<HFBasicDecorator>(*this));
        ret->wrapped = wrapped->Clone();
        return ret;
    }

    /** Extend/reduce direct potential to match lattice size.
        Note that wrapped object will itself be registered, so there is no need to pass this along.
     */
    virtual void Alert() override {}

    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(int kappa, double energy, const SpinorFunction* exchange = NULL) override
    {   wrapped->SetODEParameters(kappa, energy, exchange);
        currentEnergy = energy;
        currentKappa = kappa;
    }

    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(const Orbital& approximation) override
    {   wrapped->SetODEParameters(approximation);
        currentEnergy = approximation.Energy();
        currentKappa = approximation.Kappa();
    }

    /** Get exchange (nonlocal) potential. */
    virtual SpinorFunction GetExchange(pOrbitalConst approximation = pOrbitalConst()) const override
    {   return wrapped->GetExchange(approximation);
    }

    /** Tell SpinorODE whether to include the nonlocal (w_const) terms in GetODEFunction, GetODECoefficients, and GetODEJacobian. */
    virtual void IncludeExchange(bool include_exchange) override
    {   include_nonlocal = include_exchange;
        wrapped->IncludeExchange(include_exchange);
    }

    /** Get df/dr = w[0] and dg/dr = w[1] given point r, (f, g).
        PRE: w should be an allocated 2 dimensional array.
     */
    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const override
    {   wrapped->GetODEFunction(latticepoint, fg, w);
    }

    /** Get numerical coefficients of the ODE at the point r, (f,g).
        PRE: w_f, w_g, and w_const should be allocated 2 dimensional arrays.
     */
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const override
    {   wrapped->GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);
    }

    /** Get Jacobian (dw[i]/df and dw[i]/dg), and dw[i]/dr at a point r, (f, g).
        PRE: jacobian should be an allocated 2x2 matrix,
        dwdr should be an allocated 2 dimensional array.
     */
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const override
    {   wrapped->GetODEJacobian(latticepoint, fg, jacobian, dwdr);
    }

    /** Get approximation to eigenfunction for first numpoints near the origin. */
    virtual void EstimateOrbitalNearOrigin(unsigned int numpoints, SpinorFunction& s) const override
    {   wrapped->EstimateOrbitalNearOrigin(numpoints, s);
    }

    /** Get approximation to eigenfunction for last numpoints far from the origin.
        This routine can change the size of the orbital.
     */
    virtual void EstimateOrbitalNearInfinity(unsigned int numpoints, Orbital& s) const override
    {   wrapped->EstimateOrbitalNearInfinity(numpoints, s);
    }

public:
    /** Potential = t | a > for an operator t. */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override
    {   return wrapped->ApplyTo(a);
    }

protected:
    pHFOperator wrapped;
};

/** HFOperatorDecorator is a base class for HFOperator decorators.
    It follows the curiously recurring template pattern to provide Clone() function.
    When deriving from it use the syntax
    class MyDecorator : public HFOperatorDecorator<BaseDecorator, MyDecorator>
 */
template <typename Base, typename Derived>
class HFOperatorDecorator : public Base
{
public:
    using Base::Base;

    HFOperatorDecorator(const Derived& other):
        Base(other)
    {}

    virtual pHFOperator Clone() const override
    {
        std::shared_ptr<Derived> ret = std::make_shared<Derived>(static_cast<Derived const &>(*this));
        ret->wrapped = this->wrapped->Clone();
        return ret;
    }

protected:
    typedef HFOperatorDecorator<Base, Derived> BaseDecorator;
};

#endif
