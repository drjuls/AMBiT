#ifndef RPA_OPERATOR_H
#define RPA_OPERATOR_H

#include "TimeDependentSpinorOperator.h"
#include "HartreeFock/HartreeY.h"
#include "RPAOrbital.h"
#include "RPASolver.h"

namespace Ambit
{
/** Add RPA corrections to external SpinorOperator. */
class RPAOperator : public TimeDependentSpinorOperator
{
public:
    RPAOperator(pSpinorOperator external, pHFOperatorConst hf, pHartreeY hartreeY, pRPASolver rpa_solver);

    /** Returns true if the external operator is frequency-independent. */
    bool IsStaticRPA() const { return static_rpa; }

    /** Solve the RPA equations using solver. */
    void SolveRPA();

    virtual void SetFrequency(double frequency) override;

    void SetScale(double factor) { scale = factor; }
    double GetScale() const { return scale; }

    /** Get RPA core. */
    pCore GetRPACore() { return core; }

    /** Set RPA core. */
    void SetRPACore(pCore rpa_core) { core = rpa_core; }

    /** Clear all deltaOrbitals in the RPA core. */
    void ClearRPACore();

    /** Return (f + deltaVhf)||a> */
    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a, int kappa_b) const override;

    /** Return (f + deltaVhf)^{\dagger}||a> */
    SpinorFunction ConjugateReducedApplyTo(const SpinorFunction& a, int kappa_b) const override;

    /** Return direct expectation value <deltaVhf> */
    RadialFunction GetRPAField() const;

protected:
    SpinorFunction ReducedApplyTo(const SpinorFunction& a, int kappa_b, bool conjugate) const;

protected:
    pSpinorOperator external;
    bool static_rpa;
    double scale;

    pCore core;
    pHFOperatorConst hf;
    pHartreeY hartreeY;
    pRPASolver solver;
};

typedef std::shared_ptr<RPAOperator> pRPAOperator;
typedef std::shared_ptr<const RPAOperator> pRPAOperatorConst;

}
#endif
