#ifndef TIME_DEPENDENT_SPINOR_OPERATOR_H
#define TIME_DEPENDENT_SPINOR_OPERATOR_H

#include "HartreeFock/SpinorOperator.h"

/** Time-dependent external operators oscillate with a frequency omega:
    \f[
        \hat h = f \exp(-i \omega t) + f^\dagger \exp(i \omega t)
    \f]
    The frequency is set in atomic energy units (Hartree).
 */
class TimeDependentSpinorOperator : public SpinorOperator
{
public:
    TimeDependentSpinorOperator(int K, Parity P, pIntegrator integration_strategy = nullptr):
        SpinorOperator(K, P, integration_strategy), omega(0.0)
    {}
    TimeDependentSpinorOperator(int K, pIntegrator integration_strategy = nullptr):
        SpinorOperator(K, integration_strategy), omega(0.0)
    {}

    virtual void SetFrequency(double frequency) { omega = frequency; }
    double GetFrequency() const { return omega; }

    /** Return f^{\dagger}||a> */
    virtual SpinorFunction ConjugateReducedApplyTo(const SpinorFunction& a, int kappa_b) const
    {
        // Default is that the operator is Hermitian.
        return ReducedApplyTo(a, kappa_b);
    }

protected:
    double omega;   //!< Current value of transition frequency (energy in Hartree units) being calculated
};

typedef std::shared_ptr<TimeDependentSpinorOperator> pTimeDependentSpinorOperator;
typedef std::shared_ptr<const TimeDependentSpinorOperator> pTimeDependentSpinorOperatorConst;

#endif
