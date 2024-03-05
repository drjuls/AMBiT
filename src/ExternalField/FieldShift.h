#ifndef FIELD_SHIFT_H
#define FIELD_SHIFT_H

#include "HartreeFock/LocalPotentialDecorator.h"
#include "Transitions.h"

namespace Ambit
{
class FieldShiftDecorator: public HFOperatorDecorator<LocalPotentialDecorator, FieldShiftDecorator>
{
public:
    /** Initialise potential with nuclear radius R0 (fm) and thickness parameter t, and changes dR and dt. */
    FieldShiftDecorator(pHFOperator wrapped_hf, double R0, double dR, double t = 2.3, double dt = 0.0, pIntegrator integration_strategy = nullptr);

    double GetDeltaRMSRadius() const { return dRrms; }
    double GetDeltaR4() const { return dR4; }

protected:
    double dRrms;
    double dR4;
};

class FieldShiftCalculator : public TransitionCalculator
{
public:
    FieldShiftCalculator(MultirunOptions& user_input, Atom& atom);

    virtual void PrintHeader() const override;
    virtual void PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const override;

protected:
    double dR2;
    double dR4;
};

}
#endif
