#ifndef BREIT_ZERO_H
#define BREIT_ZERO_H

#include "HartreeFock/HartreeY.h"

/** Reduced two-body operator for frequency-independent Breit interaction, following Johnson book p. 203.
 */
class BreitZero : public HartreeYDecorator
{
public:
    BreitZero(pHartreeY wrapped, pOPIntegrator integration_strategy, pCoulombOperator coulomb):
        HartreeYDecorator(wrapped, integration_strategy), coulomb(coulomb)
    {   two_body_reverse_symmetry_exists = false;
    }

    virtual bool SetParameters(int new_K, const SpinorFunction& c, const SpinorFunction& d) override;

    /** Should be overwritten to take boolean AND of all wrapped objects. */
    virtual bool isZero() const override;

    /** Deep copy of the HartreeY object, including wrapped objects. */
    virtual HartreeYDecorator* Clone() const override
    {   pHartreeY wrapped_clone(component->Clone());
        return new BreitZero(wrapped_clone, integrator, coulomb);
    }

    /** < b | t | a > for an operator t. */
    virtual double GetMatrixElement(const Orbital& b, const SingleParticleWavefunction& a, bool reverse) const override;

    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse) const override;

protected:
    /** P_{ik}(r) for current value of K. */
    SpinorFunction P(const SpinorFunction& k, int kappa_i) const;
    RadialFunction P(const SpinorFunction& i, const SpinorFunction& k) const;

    /** Q_{ik}(r) for current value of K. */
    SpinorFunction Q(const SpinorFunction& k, int kappa_i) const;
    RadialFunction Q(const SpinorFunction& i, const SpinorFunction& k) const;

protected:
    pCoulombOperator coulomb;

    /* Integral of the form
                             r< ^L
          I(r) = Integral[ ---------.P_cd(r').dr' ] . <kappa_c || C^L || kappa_d>
                           r> ^(L+1)
       with L = K+1. Others are similar.
     */
    RadialFunction pot_P_Kplus;
    RadialFunction pot_P_Kminus;
    RadialFunction pot_Q_Kplus;
    RadialFunction pot_Q_Kminus;
    RadialFunction pot_V_K;

    bool jc_plus_jd_even;
};

#endif
