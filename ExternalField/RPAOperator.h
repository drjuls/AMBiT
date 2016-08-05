#ifndef RPA_OPERATOR_H
#define RPA_OPERATOR_H

#include "TimeDependentSpinorOperator.h"
#include "HartreeFock/HartreeY.h"
#include "RPAOrbital.h"
#include "RPASolver.h"

/** Add RPA corrections to external SpinorOperator. */
class RPAOperator : public TimeDependentSpinorOperator
{
public:
    RPAOperator(pSpinorOperator external, pHFOperatorConst hf, pHartreeY hartreeY, pRPASolver rpa_solver);

    /** Solve the RPA equations using solver. */
    void SolveRPA();

    /** Reset frequency and solve RPA core.  */
    virtual void SetFrequency(double frequency) override;

    /** Get RPA core. */
    pCore GetRPACore() { return core; }

    /** Set RPA core. */
    void SetRPACore(pCore rpa_core) { core = rpa_core; }

    /** Return (f + deltaVhf)||a> */
    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a, int kappa_b) const override;

protected:
    pSpinorOperator external;
    bool static_rpa;

    pCore core;
    pHFOperatorConst hf;
    pHartreeY hartreeY;
    pRPASolver solver;
};

typedef std::shared_ptr<RPAOperator> pRPAOperator;
typedef std::shared_ptr<const RPAOperator> pRPAOperatorConst;

#endif
