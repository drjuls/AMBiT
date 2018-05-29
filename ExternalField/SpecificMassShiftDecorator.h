#ifndef SPECIFIC_MASS_SHIFT_DECORATOR_H
#define SPECIFIC_MASS_SHIFT_DECORATOR_H

#include "HartreeFock/ExchangeDecorator.h"
#include "TwoBodySMSOperator.h"

namespace Ambit
{
/** Add specific mass shift to exchange part of operator [Berengut et al. PRA 68, 022502 (2003)].
    \f{eqnarray*}{
        P &=& p_{an}p_{na} \quad \textrm{where} \\
        p_{ab} &=& \int f_a \left( \frac{d}{dr} - \frac{l_a}{r} \right) f_b dr . \delta_{l_a, l_b+1}
                   + \int f_a \left( \frac{d}{dr} + \frac{l_b}{r} \right) f_b dr . \delta_{l_a, l_b-1} \\
               &=& - p_{ba}
    \f}
    \f$p_{ab}\f$ is calculated using TwoBodySMSOperator.
    The extra exchange is stored in currentExchangePotential, inherited from HFDecorator.
    Typical values of inverse mass (1/M) are of order 0.001.
 */
class SpecificMassShiftDecorator : public HFOperatorDecorator<ExchangeDecorator, SpecificMassShiftDecorator>
{
public:
    SpecificMassShiftDecorator(pHFOperator wrapped_hf, bool nonrel = false, bool nonrel_include_lower = true);

    /** Set the inverse nuclear mass: 1/M. */
    void SetInverseMass(double InverseNuclearMass)
    {   sms_operator->SetInverseMass(InverseNuclearMass);
    }

    double GetInverseMass() const { return sms_operator->GetInverseMass(); }

protected:
    virtual SpinorFunction CalculateExtraExchange(const SpinorFunction& s) const override;

protected:
    pSMSOperator sms_operator;
};

typedef std::shared_ptr<SpecificMassShiftDecorator> pSpecificMassShiftDecorator;
typedef std::shared_ptr<const SpecificMassShiftDecorator> pSpecificMassShiftDecoratorConst;

}
#endif
