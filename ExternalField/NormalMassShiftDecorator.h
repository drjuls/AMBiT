#ifndef NORMAL_MASS_SHIFT_DECORATOR
#define NORMAL_MASS_SHIFT_DECORATOR

#include "HartreeFock/HFOperator.h"

/** Add normal mass shift to exchange part of operator.
    The extra exchange is stored in currentExchangePotential, inherited from HFDecorator.
    Typical values of inverse mass (1/M) are of order 0.00001 or even smaller.
 */
class NormalMassShiftDecorator : public HFOperatorDecorator<HFBasicDecorator, NormalMassShiftDecorator>
{
public:
    NormalMassShiftDecorator(pHFOperator wrapped_hf, bool only_rel_nms = false, bool nonrel = false);

    /** Set the inverse nuclear mass: 1/M. */
    void SetInverseMass(double InverseNuclearMass) { lambda = InverseNuclearMass; }
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
    bool do_rel_nms;
    bool do_nonrel_nms;
};

typedef std::shared_ptr<NormalMassShiftDecorator> pNormalMassShiftDecorator;
typedef std::shared_ptr<const NormalMassShiftDecorator> pNormalMassShiftDecoratorConst;

#endif
