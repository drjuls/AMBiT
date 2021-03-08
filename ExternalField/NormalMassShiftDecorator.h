#ifndef NORMAL_MASS_SHIFT_DECORATOR_H
#define NORMAL_MASS_SHIFT_DECORATOR_H

#include "HartreeFock/ExchangeDecorator.h"
#include "ExternalField/Transitions.h"

namespace Ambit
{
/** Add normal mass shift to exchange part of operator.
    The extra exchange is stored in currentExchangePotential, inherited from HFDecorator.
    Typical values of inverse mass (1/M) are of order 0.00001 or even smaller.
 */
class NormalMassShiftDecorator : public HFOperatorDecorator<ExchangeDecorator, NormalMassShiftDecorator>
{
public:
    NormalMassShiftDecorator(pHFOperator wrapped_hf, bool only_rel_nms = false, bool nonrel = false);

    /** Set the inverse nuclear mass: 1/M. */
    inline void SetInverseMass(double InverseNuclearMass)
    {   SetScale(InverseNuclearMass);
    }

    inline double GetInverseMass() const
    {   return GetScale();
    }

    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const override;
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const override;
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const override;

    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;

protected:
    virtual SpinorFunction CalculateExtraExchange(const SpinorFunction& s) const override;

protected:
    bool do_rel_nms;
    bool do_nonrel_nms;
};

typedef std::shared_ptr<NormalMassShiftDecorator> pNormalMassShiftDecorator;
typedef std::shared_ptr<const NormalMassShiftDecorator> pNormalMassShiftDecoratorConst;

class NormalMassShiftOperator : public SpinorOperator
{
public:
    NormalMassShiftOperator(pHFOperator hf, bool only_rel_nms = false, bool nonrel = true);

    virtual SpinorFunction ReducedApplyTo(const SpinorFunction& a, int kappa_b) const override;

protected:
    pHFOperator hf;
    bool do_rel_nms;
    bool do_nonrel_nms;
};

class NormalMassShiftCalculator : public TransitionCalculator
{
public:
    NormalMassShiftCalculator(MultirunOptions& user_input, Atom& atom);

    virtual void PrintHeader() const override;
    virtual void PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const override;
};

}
#endif
