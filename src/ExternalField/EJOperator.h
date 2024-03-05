#ifndef EJ_OPERATOR_H
#define EJ_OPERATOR_H

#include "TimeDependentSpinorOperator.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/Enums.h"
#include "Transitions.h"

namespace Ambit
{
class EJOperator : public TimeDependentSpinorOperator
{
public:
    EJOperator(int J, pIntegrator integration_strategy, TransitionGauge gauge = TransitionGauge::Length):
        TimeDependentSpinorOperator(J, (J%2? Parity::odd: Parity::even), integration_strategy), gauge(gauge)
    {}

    void SetGauge(TransitionGauge gauge_type) { gauge = gauge_type; }
    TransitionGauge GetGauge() const { return gauge; }

    /** Reduced matrix element < b || E(J) || a > for our operator E(J). */
//    virtual double GetReducedMatrixElement(const Orbital& b, const Orbital& a) const override;

    /** E(J) || a > for our operator E(J). */
    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a, int kappa_b) const override
    {
        return ReducedApplyTo(a, kappa_b, false);
    }

    /** E(J)^{\dagger} || a > for our operator E(J). */
    virtual SpinorFunction ConjugateReducedApplyTo(const SpinorFunction& a, int kappa_b) const override
    {
        return ReducedApplyTo(a, kappa_b, true);
    }

protected:
    /** E(J) || a > for our operator E(J). */
    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a, int kappa_b, bool conjugate) const;

    TransitionGauge gauge;
};

class MJOperator : public TimeDependentSpinorOperator
{
public:
    MJOperator(int J, pIntegrator integration_strategy):
        TimeDependentSpinorOperator(J, (J%2? Parity::even: Parity::odd), integration_strategy)
    {}

    /** M(J) || a > for our operator M(J). */
    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a, int kappa_b) const override;
};

class EMCalculator : public TransitionCalculator
{
public:
    EMCalculator(MultipolarityType type, int J, MultirunOptions& user_input, Atom& atom);

    virtual void PrintHeader() const override;
    virtual void PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const override;

    virtual std::string Name() const;

protected:
    MultipolarityType type;
    int J;
};

}
#endif
