#ifndef NORMAL_MASS_SHIFT_DECORATOR
#define NORMAL_MASS_SHIFT_DECORATOR

#include "HartreeFock/ExchangeDecorator.h"

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

protected:
    virtual SpinorFunction CalculateExtraExchange(const SpinorFunction& s) const;

protected:
    bool do_rel_nms;
    bool do_nonrel_nms;
};

typedef std::shared_ptr<NormalMassShiftDecorator> pNormalMassShiftDecorator;
typedef std::shared_ptr<const NormalMassShiftDecorator> pNormalMassShiftDecoratorConst;

#endif
