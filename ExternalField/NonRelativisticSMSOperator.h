#ifndef NONRELATIVISTIC_SMS_OPERATOR_H
#define NONRELATIVISTIC_SMS_OPERATOR_H

#include "HartreeFock/HartreeY.h"
#include "Universal/FornbergDifferentiator.h"

/** Non-relativistic specific mass shift operator [Berengut et al. PRA 68, 022502 (2003)].
    \f[
        P^1 (ac, bd) = \lambda p_{ab} p_{cd}
    \f]
    where
    \f[
        p_{ab} = \int f_a \left( \frac{d}{dr} - \frac{l_a}{r} \right) f_b dr . \,\delta_{l_a, l_b+1}
                + \int f_a \left( \frac{d}{dr} + \frac{l_b}{r} \right) f_b dr . \,\delta_{l_a, l_b-1}
               = - p_{ba}
    \f]
    Multipolarity K == 1.
 */
class NonRelativisticSMSOperator : public HartreeYDecorator<HartreeYBasicDecorator, NonRelativisticSMSOperator>
{
public:
    NonRelativisticSMSOperator(pHartreeY wrapped, pIntegrator integration_strategy = nullptr):
        BaseDecorator(wrapped, integration_strategy), lambda(0.0), p_cd(0.0), interp(integrator->GetLattice())
    {   two_body_reverse_symmetry_exists = false;
    }

    /** Set the inverse nuclear mass: 1/M. */
    void SetInverseMass(double InverseNuclearMass) { lambda = InverseNuclearMass; }
    double GetInverseMass() const { return lambda; }

    /** < b | t | a > for an operator t. */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a, bool reverse = false) const override;

    /** Potential = t | a > for an operator t such that the resulting Potential.Kappa() == kappa_b.
        i.e. t | a > has kappa == kappa_b.
     */
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b, bool reverse = false) const override;

protected:
    virtual bool SetLocalParameters(int new_K, pSpinorFunctionConst new_c, pSpinorFunctionConst new_d) override;
    virtual bool isZeroLocal() const override;
    virtual int GetLocalMinK() const override;
    virtual int GetLocalMaxK() const override;

    /** Returns \f$ \hat p f_a \f$. */
    SpinorFunction ApplyOperator(const SpinorFunction& a, int kappa_b) const;

protected:
    double lambda;
    double p_cd;
    FornbergDifferentiator interp;
};

typedef std::shared_ptr<NonRelativisticSMSOperator> pSMSOperator;

#endif
