#ifndef HF_OPERATOR_H
#define HF_OPERATOR_H

#include "OneBodyOperator.h"
#include "SpinorODE.h"
#include "Core.h"
#include "CoulombOperator.h"

/** The relativistic Hartree-Fock (Dirac-Fock) operator:
    \f[
    t = \left( \begin{array}{cc}
                  -V             & (-d/dr + \kappa/r)/\alpha \\
        (d/dr + \kappa/r)/\alpha &     -2/\alpha^2 - V
        \end{array} \right)
    \f]
    where V is the electrostatic potential (V > 0) which leads to -V for electrons.
 */
class HFOperator : public OneBodyOperator, public SpinorODE
{
public:
    HFOperator(double Z, pCoreConst hf_core, pOPIntegrator integration_strategy, pCoulombOperator coulomb);
    HFOperator(const HFOperator& other);
    virtual ~HFOperator();

    /** Set/reset the Hartree-Fock core, from which the potential is derived. */
    virtual void SetCore(pCoreConst hf_core);
    virtual pCoreConst GetCore() const;

    virtual RadialFunction GetDirectPotential() const;  //!< Get the direct potential.
    virtual double GetZ() const { return Z; }           //!< Get nuclear charge.
    virtual double GetCharge() const { return charge; } //!< Get ion charge.

    /** Get size of valid latticepoints. */
    virtual unsigned int size() const;

    /** Extend direct potential to match lattice size. */
    virtual void ExtendPotential();

    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(int kappa, double energy, SpinorFunction* exchange = NULL) override;
    
    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(const SingleParticleWavefunction& approximation) override;

    /** Get exchange (nonlocal) potential. */
    virtual SpinorFunction GetExchange(pSingleParticleWavefunctionConst approximation = pSingleParticleWavefunctionConst()) const override;

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
    /** Potential = t | a > for an operator t. */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const override;

protected:
    virtual SpinorFunction CalculateExchange(const SpinorFunction& s) const;

protected:
    double Z;
    double charge;
    pCoreConst core;
    pCoulombOperator coulombSolver;

    RadialFunction directPotential;
    SpinorFunction currentExchangePotential;
    double currentEnergy;
    int currentKappa;
};

typedef boost::shared_ptr<HFOperator> pHFOperator;
typedef boost::shared_ptr<const HFOperator> pHFOperatorConst;

class HFOperatorDecorator : public HFOperator
{
public:
    HFOperatorDecorator(pHFOperator decorated_object, pOPIntegrator integration_strategy = pOPIntegrator()):
        HFOperator(*decorated_object), wrapped(decorated_object)
    {   if(integration_strategy != NULL)
            integrator = integration_strategy;
    }
    virtual ~HFOperatorDecorator() {}

    /** Set/reset the Hartree-Fock core, from which the potential is derived. */
    virtual void SetCore(pCoreConst hf_core) override
    {   wrapped->SetCore(hf_core);
        core = hf_core;
        directPotential.Clear();
    }

    virtual RadialFunction GetDirectPotential() const override
    {   return wrapped->GetDirectPotential();
    }

    /** Get size of valid latticepoints. */
    virtual unsigned int size() const override
    {   return wrapped->size();
    }

    /** Extend direct potential to match lattice size. */
    virtual void ExtendPotential() override
    {   wrapped->ExtendPotential();
    }

    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(int kappa, double energy, SpinorFunction* exchange = NULL) override
    {   wrapped->SetODEParameters(kappa, energy, exchange);
        currentEnergy = energy;
        currentKappa = kappa;
        currentExchangePotential.Clear();
    }

    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(const SingleParticleWavefunction& approximation) override
    {   wrapped->SetODEParameters(approximation);
        currentEnergy = approximation.Energy();
        currentKappa = approximation.Kappa();
        currentExchangePotential.Clear();
    }

    /** Get exchange (nonlocal) potential. */
    virtual SpinorFunction GetExchange(pSingleParticleWavefunctionConst approximation = pSingleParticleWavefunctionConst()) const override
    {   return wrapped->GetExchange(approximation);
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

#endif
