#ifndef BREIT_ZERO_H
#define BREIT_ZERO_H

#include "HartreeFock/HartreeY.h"

/** Two-body operator for frequency-independent Breit interaction, following Johnson book p. 203.
    Note that here the reduced spherical tensor matrix elements (prefactors in Johnson's book) are
    removed: we are calculating the equivalent of HartreeY operators.
 */
class BreitZero : public HartreeYDecorator<HartreeYBasicDecorator, BreitZero>
{
public:
    BreitZero(pHartreeY wrapped, pIntegrator integration_strategy, pCoulombOperator coulomb):
        BaseDecorator(wrapped, integration_strategy), coulomb(coulomb)
    {   two_body_reverse_symmetry_exists = false;
        off_parity_exists = true;
    }

    /** < b | t | a > for an operator t. */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a, bool reverse = false) const override;

    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse = false) const override;

protected:
    virtual bool SetLocalParameters(int new_K, pSpinorFunctionConst new_c, pSpinorFunctionConst new_d) override;
    virtual bool isZeroLocal() const override;
    virtual int GetLocalMinK() const override;
    virtual int GetLocalMaxK() const override;

    /** P_{ik}(r) for current value of K.
        PRE: K > 0.
     */
    SpinorFunction P(const SpinorFunction& k, int kappa_i) const;
    RadialFunction P(const SpinorFunction& i, const SpinorFunction& k) const;

    /** Q_{ik}(r) for current value of K. */
    SpinorFunction Q(const SpinorFunction& k, int kappa_i) const;
    RadialFunction Q(const SpinorFunction& i, const SpinorFunction& k) const;

protected:
    pCoulombOperator coulomb;

    /* Integral of the form
                             r< ^L
          I(r) = Integral[ ---------.P_cd(r').dr' ]
                           r> ^(L+1)
       with L = K+1. Others are similar.
     */
    RadialFunction pot_P_Kminus;
    RadialFunction pot_Q_Kplus;
    RadialFunction pot_V_K;

    /* Integral of the form
                          r< ^(K-1)   r< ^(K+1)
        I1(r) = Integral[ --------- - --------- .P_cd(r').dr', {r', 0, r} ]
                            r> ^K     r> ^(K+2)
     */
    RadialFunction pot_P_fwd;
    /* Integral of the form
                          r< ^(K-1)   r< ^(K+1)
        I2(r) = Integral[ --------- - --------- .Q_cd(r').dr', {r', r, Infinity} ]
                            r> ^K     r> ^(K+2)
     */
    RadialFunction pot_Q_back;
};

#endif
