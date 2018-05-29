#ifndef YUKAWA_POTENTIAL_H
#define YUKAWA_POTENTIAL_H

#include "HartreeFock/LocalPotentialDecorator.h"
#include "Transitions.h"

namespace Ambit
{
/** Yukawa decorator adds the potential
        V(r) = alpha.exp(-mr)/r
    where alpha and m are constants.
    m is in atomic units (m_e = 1)
 */
class YukawaDecorator: public HFOperatorDecorator<LocalPotentialDecorator, YukawaDecorator>
{
public:
    YukawaDecorator(pHFOperator wrapped_hf, double mass, double scale = 1., pIntegrator integration_strategy = pIntegrator());

    /** Generate Yukawa potential for given mass m. */
    void GenerateYukawaPotential(double new_mass);

protected:
    double mass; // mass (units of m_e)
    virtual void Alert() override;
};

class YukawaCalculator : public TransitionCalculator
{
public:
    YukawaCalculator(MultirunOptions& user_input, Atom& atom);

    virtual void PrintHeader() const override;
    virtual void PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const override;

protected:
    double mass;
};

}
#endif
