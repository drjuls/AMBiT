#ifndef NONRELATIVISTIC_SMS_OPERATOR_H
#define NONRELATIVISTIC_SMS_OPERATOR_H

#include "Operator.h"

/** Non-relativistic specific mass shift operator [Berengut et al. PRA 68, 022502 (2003)].
    \f[
        p_{ab} = \int f_a \left( \frac{d}{dr} - \frac{l_a}{r} \right) f_b dr . \delta_{l_a, l_b+1}
                + \int f_a \left( \frac{d}{dr} + \frac{l_b}{r} \right) f_b dr . \delta_{l_a, l_b-1}
               = - p_{ba}
    \f]
 */
class NonRelativisticSMSOperator : public OneBodyOperatorDecorator
{
public:
    NonRelativisticSMSOperator(pOneBodyOperator wrapped, pOPIntegrator integration_strategy):
        OneBodyOperatorDecorator(wrapped, integration_strategy)
    {}

    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;
    virtual SpinorFunction ApplyTo(const SpinorFunction& a, int kappa_b) const override;
};


#endif
