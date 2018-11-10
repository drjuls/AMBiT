#ifndef KINETIC_ENERGY_H
#define KINETIC_ENERGY_H

#include "HartreeFock/ExchangeDecorator.h"
#include "Transitions.h"

namespace Ambit
{
/** Add additional kinetic energy term to HF operator.
 */
class KineticEnergyDecorator : public HFOperatorDecorator<ExchangeDecorator, KineticEnergyDecorator>
{
public:
    KineticEnergyDecorator(pHFOperator wrapped_hf, pIntegrator integration_strategy = nullptr):
        BaseDecorator(wrapped_hf, integration_strategy)
    {}

    virtual SpinorFunction CalculateExtraExchange(const SpinorFunction& s) const;
};

class KineticEnergyCalculator : public TransitionCalculator
{
public:
    KineticEnergyCalculator(MultirunOptions& user_input, Atom& atom);

    virtual void PrintHeader() const override;
    virtual void PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const override;
};

}

#endif
