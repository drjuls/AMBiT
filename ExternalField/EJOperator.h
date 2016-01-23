#ifndef EJ_OPERATOR_H
#define EJ_OPERATOR_H

#include "Universal/SpinorMatrixElement.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/Enums.h"
#include "Atom/Transitions.h"

class EJOperator : public SpinorMatrixElement
{
public:
    EJOperator(int J, pIntegrator integration_strategy, TransitionGauge gauge = TransitionGauge::Length):
        SpinorMatrixElement(J, (J%2? Parity::odd: Parity::even), integration_strategy), gauge(gauge)
    {}

    virtual void SetGauge(TransitionGauge gauge_type) { gauge = gauge_type; }

    /** Reduced matrix element < b || E(J) || a > for our operator E(J). */
    virtual double GetReducedMatrixElement(const Orbital& b, const Orbital& a) const override;

protected:
    TransitionGauge gauge;
};

class MJOperator : public SpinorMatrixElement
{
public:
    MJOperator(int J, pIntegrator integration_strategy):
        SpinorMatrixElement(J, (J%2? Parity::even: Parity::odd), integration_strategy)
    {}

    /** Reduced matrix element < b || E(J) || a > for our operator E(J). */
    virtual double GetReducedMatrixElement(const Orbital& b, const Orbital& a) const override;
};

class EMCalculator : public TransitionCalculator
{
public:
    EMCalculator(MultipolarityType type, int J, MultirunOptions& user_input, pOrbitalManagerConst orbitals, pLevelStore levels, pIntegrator integrator);

    virtual void PrintHeader() const override;
    virtual void PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const override;

protected:
    MultipolarityType type;
    int J;
};

#endif
