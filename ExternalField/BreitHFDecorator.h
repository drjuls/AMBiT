#ifndef BREIT_HF_DECORATOR_H
#define BREIT_HF_DECORATOR_H

#include "HartreeFock/HFOperator.h"
#include "BreitZero.h"

/** As defined by Johnson
        \f$ (B_HF)_{ij} = Sum_b \left( b_{ibjb} - b_{ibbj} \right) \f$,
    where \f$ b_{ijkl} \f$ is the Breit operator.
    The direct part \f$ b_{ibjb} \f$ is zero, so only the exchange part contributes.
 */
class BreitHFDecorator : public HFOperatorDecorator<HFBasicDecorator, BreitHFDecorator>
{
public:
    BreitHFDecorator(pHFOperator wrapped_hf, pHartreeY breit_operator):
        BaseDecorator(wrapped_hf), breit_operator(breit_operator)
    {}
    BreitHFDecorator(const BreitHFDecorator& other):
        HFOperatorDecorator(other), breit_operator(other.breit_operator)
    {}

    virtual void Alert() override;

    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(const Orbital& approximation) override;

    /** Get exchange (nonlocal) potential. */
    virtual SpinorFunction GetExchange(pOrbitalConst approximation) const override;

    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const override;
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const override;
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const override;

public:
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a) const override;
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;

protected:
    virtual SpinorFunction CalculateExtraExchange(const SpinorFunction& s) const;
    pHartreeY breit_operator;
};

#endif
