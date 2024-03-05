#ifndef BREIT_HF_DECORATOR_H
#define BREIT_HF_DECORATOR_H

#include "HartreeFock/ExchangeDecorator.h"
#include "BreitZero.h"

namespace Ambit
{
/** As defined by Johnson
        \f$ (B_HF)_{ij} = Sum_b \left( b_{ibjb} - b_{ibbj} \right) \f$,
    where \f$ b_{ijkl} \f$ is the Breit operator.
    The direct part \f$ b_{ibjb} \f$ is zero, so only the exchange part contributes.
 */
class BreitHFDecorator : public HFOperatorDecorator<ExchangeDecorator, BreitHFDecorator>
{
public:
    BreitHFDecorator(pHFOperator wrapped_hf, pHartreeY breit_operator):
        BaseDecorator(wrapped_hf), breit_operator(breit_operator)
    {}

public:
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a) const override;

    /** Deep copy including breit_operator. */
    virtual pHFOperator Clone() const override
    {
        auto ret(std::make_shared<BreitHFDecorator>(*this));
        ret->wrapped = wrapped->Clone();
        ret->breit_operator = breit_operator->Clone();
        return ret;
    }

protected:
    virtual SpinorFunction CalculateExtraExchange(const SpinorFunction& s) const override;
    pHartreeY breit_operator;
};

}
#endif
