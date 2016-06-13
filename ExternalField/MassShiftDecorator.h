#ifndef MASS_SHIFT_DECORATOR
#define MASS_SHIFT_DECORATOR

#include "HartreeFock/HFOperator.h"
#include "NonRelativisticSMSOperator.h"

/** Add non-relativistic specific mass shift to exchange part of operator [Berengut et al. PRA 68, 022502 (2003)].
    \f{eqnarray*}{
        P &=& p_{an}p_{na} \quad \textrm{where} \\
        p_{ab} &=& \int f_a \left( \frac{d}{dr} - \frac{l_a}{r} \right) f_b dr . \delta_{l_a, l_b+1}
                   + \int f_a \left( \frac{d}{dr} + \frac{l_b}{r} \right) f_b dr . \delta_{l_a, l_b-1} \\
               &=& - p_{ba}
    \f}
    \f$p_{ab}\f$ is calculated using NonRelativisticSMSOperator.
    The extra exchange is stored in currentExchangePotential, inherited from HFDecorator.
    Typical values of inverse mass (1/M) are of order 0.001.
 */
class MassShiftDecorator : public HFOperatorDecorator<HFBasicDecorator, MassShiftDecorator>
{
public:
    MassShiftDecorator(pHFOperator wrapped_hf, bool relativistic_version = false, bool include_NMS = false, bool include_SMS = true);

    /** Set the inverse nuclear mass: 1/M. */
    void SetInverseMass(double InverseNuclearMass)
    {   lambda = InverseNuclearMass;
        sms_operator->SetInverseMass(InverseNuclearMass);
    }

    double GetInverseMass() const { return lambda; }

public:
    virtual void Alert() override;

    /** Set exchange (nonlocal) potential and energy for ODE routines. */
    virtual void SetODEParameters(const Orbital& approximation) override;
    
    /** Get exchange (nonlocal) potential. */
    virtual SpinorFunction GetExchange(pOrbitalConst approximation) const override;

    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const override;
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const override;
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const override;

public:
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;

protected:
    virtual SpinorFunction CalculateExtraExchange(const SpinorFunction& s) const;

protected:
    double lambda;
    bool include_nms;
    bool include_relativistic_nms;
    pSMSOperator sms_operator;
};

typedef std::shared_ptr<MassShiftDecorator> pMassShiftDecorator;
typedef std::shared_ptr<const MassShiftDecorator> pMassShiftDecoratorConst;

#endif
