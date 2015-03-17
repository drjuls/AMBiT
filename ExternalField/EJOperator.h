#ifndef EJ_OPERATOR_H
#define EJ_OPERATOR_H

#include "Universal/SpinorMatrixElement.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/Enums.h"

class EJOperator : public SpinorMatrixElementDecorator
{
public:
    EJOperator(pPhysicalConstant constants, int J, pOPIntegrator integration_strategy, TransitionGauge gauge = TransitionGauge::Length):
        SpinorMatrixElementDecorator(J, pSpinorMatrixElement(new SpinorMatrixElement()), integration_strategy),
        constants(constants), gauge(gauge)
    {}

    virtual void SetGauge(TransitionGauge gauge_type) { gauge = gauge_type; }

    /** Reduced matrix element < b || E(J) || a > for our operator E(J). */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a) const override;

protected:
    TransitionGauge gauge;
    pPhysicalConstant constants;
};

class MJOperator : public SpinorMatrixElementDecorator
{
public:
    MJOperator(pPhysicalConstant constants, int J, pOPIntegrator integration_strategy):
        SpinorMatrixElementDecorator(J, pSpinorMatrixElement(new SpinorMatrixElement()), integration_strategy),
        constants(constants)
    {}

    /** Reduced matrix element < b || E(J) || a > for our operator E(J). */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a) const override;

protected:
    pPhysicalConstant constants;
};

#endif
