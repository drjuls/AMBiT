#ifndef EJ_OPERATOR_H
#define EJ_OPERATOR_H

#include "Universal/SpinorMatrixElement.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/Enums.h"

class EJOperator : public SpinorMatrixElement
{
public:
    EJOperator(pPhysicalConstant constants, int J, pIntegrator integration_strategy, TransitionGauge gauge = TransitionGauge::Length):
        SpinorMatrixElement(J, integration_strategy), constants(constants), gauge(gauge)
    {}

    virtual void SetGauge(TransitionGauge gauge_type) { gauge = gauge_type; }

    /** Reduced matrix element < b || E(J) || a > for our operator E(J). */
    virtual double GetReducedMatrixElement(const Orbital& b, const Orbital& a) const override;

protected:
    TransitionGauge gauge;
    pPhysicalConstant constants;
};

class MJOperator : public SpinorMatrixElement
{
public:
    MJOperator(pPhysicalConstant constants, int J, pIntegrator integration_strategy):
        SpinorMatrixElement(J, integration_strategy), constants(constants)
    {}

    /** Reduced matrix element < b || E(J) || a > for our operator E(J). */
    virtual double GetReducedMatrixElement(const Orbital& b, const Orbital& a) const override;

protected:
    pPhysicalConstant constants;
};

#endif
