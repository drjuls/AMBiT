#ifndef RELATIVISTIC_SMS_OPERATOR_H
#define RELATIVISTIC_SMS_OPERATOR_H

#include "NonRelativisticSMSOperator.h"

/** Relativistic specific mass shift operator
    \f[
        P^1 (ac, bd) = \lambda p_{ab}^{(r)} p_{cd}^{(r)}
    \f]
    where
    \f[
        p^{(r)} = \vec{p} - \frac{\alpha Z}{2r} \left( \vec{\alpha} + (\vec{\alpha}\cdot\hat{r})\hat{r} \right)
    \f]
    Multipolarity K == 1.
 */
class RelativisticSMSOperator : public HartreeYDecorator<NonRelativisticSMSOperator, RelativisticSMSOperator>
{
public:
    RelativisticSMSOperator(pHartreeY wrapped, double Zalpha, pIntegrator integration_strategy = nullptr):
        BaseDecorator(wrapped, integration_strategy), Zalpha(Zalpha)
    {}

    /** < b | t | a > for an operator t. */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a, bool reverse = false) const override;

protected:
    /** Returns \f$ p^{(r)} \psi_a \f$. */
    virtual SpinorFunction ApplyOperator(const SpinorFunction& a, int kappa_b) const override;

    double Zalpha;
};

#endif
