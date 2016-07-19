#ifndef FIELD_SHIFT_H
#define FIELD_SHIFT_H

#include "HartreeFock/LocalPotentialDecorator.h"
#include "Transitions.h"

class FieldShiftDecorator: public HFOperatorDecorator<LocalPotentialDecorator, FieldShiftDecorator>
{
public:
    /** Initialise with nuclear radius R0 (fm), change in nuclear radius dR, and thickness parameter t. */
    FieldShiftDecorator(pHFOperator wrapped_hf, double R0, double dR, double t = 2.3, pIntegrator integration_strategy = nullptr);

    double GetDeltaRMSRadius() const { return dRrms; }

protected:
    double dRrms;
};

class FieldShiftCalculator : public TransitionCalculator
{
public:
    FieldShiftCalculator(MultirunOptions& user_input, Atom& atom);

    virtual void PrintHeader() const override;
    virtual void PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const override;

protected:
    double MHzOverRsq;
};

#endif
