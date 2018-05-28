#ifndef HF_OPERATOR_H
#define HF_OPERATOR_H

#include "HFOperatorBase.h"
#include "CoulombOperator.h"

namespace Ambit
{
/** The relativistic Hartree-Fock (Dirac-Fock) operator:
    \f[
    h = \left( \begin{array}{cc}
                  -V             & (-d/dr + \kappa/r)/\alpha \\
        (d/dr + \kappa/r)/\alpha &     -2/\alpha^2 - V
        \end{array} \right)
    \f]
    where V is the electrostatic potential (V > 0) which leads to -V for electrons.
 */
class HFOperator : public HFOperatorBase
{
public:
    HFOperator(double Z, pCoreConst hf_core, pPhysicalConstant physical_constant, pIntegrator integration_strategy, pCoulombOperator coulomb);
    virtual ~HFOperator() {}

public:
    /** Set/reset the Hartree-Fock core, from which the potential is derived. */
    virtual void SetCore(pCoreConst hf_core) override;

    virtual RadialFunction GetDirectPotential() const override; //!< Get the direct potential.

    /** Deep copy of the HFOperator object, particularly including wrapped objects (but not the core or physical constants). */
    virtual pHFOperator Clone() const override { return std::make_shared<HFOperator>(*this); }

public:
    /** Extend/reduce direct potential to match lattice size. */
    virtual void Alert() override;

    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(int kappa, double energy, const SpinorFunction* exchange = NULL) override;
    
    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(const Orbital& approximation) override;

    /** Get exchange (nonlocal) potential. */
    virtual SpinorFunction GetExchange(pOrbitalConst approximation = nullptr) const override;

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

protected:
    SpinorFunction CalculateExchange(const SpinorFunction& s) const;

protected:
    pCoulombOperator coulombSolver;
    RadialFunction directPotential;
    SpinorFunction currentExchangePotential;
};

}
#endif
