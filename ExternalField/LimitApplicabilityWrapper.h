#ifndef LIMIT_APPLICABILITY_WRAPPER_H
#define LIMIT_APPLICABILITY_WRAPPER_H

#include "HartreeFock/HFOperator.h"

/** The idea of this class is to wrap any HFOperatorDecorator such that its applicability can be limited.
    The limitation of a method comes via the ChangeThisSpinorFunction():
        if(ChangeThisSpinorFunction(s) == true)
            call the Decorated version of the method.
        else
            call the undecorated version via wrapped, bypassing the Decorator method.
    Generally the limiting will be done in a bespoke sense (changing the code rather than input files),
    so just modify ChangeThisSpinorFunction() as required.
    To use the class, simply wrap the HFOperatorDecorator, so
        pExampleDecorator t = std::make_shared<ExampleDecorator>(constructor_arguments...)
    becomes
        pExampleDecorator t = std::make_shared<LimitApplicabilityWrapper<ExampleDecorator>>(constructor_arguments...)
 */
template <class DecoratorType>
class LimitApplicabilityWrapper : public HFOperatorDecorator<DecoratorType, LimitApplicabilityWrapper<DecoratorType>>
{
public:
    template <class... DecoratorTypeArgs>
    LimitApplicabilityWrapper(DecoratorTypeArgs... args):
        HFOperatorDecorator<DecoratorType, LimitApplicabilityWrapper<DecoratorType>>(args...)
    {}

    virtual SpinorFunction GetExchange(pOrbitalConst approximation = pOrbitalConst()) const override
    {
        if(ChangeThisSpinorFunction(*approximation.get()))
            return DecoratorType::GetExchange(approximation);
        else
            return DecoratorType::wrapped->GetExchange(approximation);
    }

    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const override
    {
        if(ChangeThisSpinorFunction(fg))
            DecoratorType::GetODEFunction(latticepoint, fg, w);
        else
            DecoratorType::wrapped->GetODEFunction(latticepoint, fg, w);
    }

    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const override
    {
        if(ChangeThisSpinorFunction(fg))
            DecoratorType::GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);
        else
            DecoratorType::wrapped->GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);
    }

    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const override
    {
        if(ChangeThisSpinorFunction(fg))
            DecoratorType::GetODEJacobian(latticepoint, fg, jacobian, dwdr);
        else
            DecoratorType::wrapped->GetODEJacobian(latticepoint, fg, jacobian, dwdr);
    }

    virtual void EstimateOrbitalNearOrigin(unsigned int numpoints, SpinorFunction& s) const override
    {
        if(ChangeThisSpinorFunction(s))
            DecoratorType::EstimateOrbitalNearOrigin(numpoints, s);
        else
            DecoratorType::wrapped->EstimateOrbitalNearOrigin(numpoints, s);
    }

    virtual void EstimateOrbitalNearInfinity(unsigned int numpoints, Orbital& s) const override
    {
        if(ChangeThisSpinorFunction(s))
            DecoratorType::EstimateOrbitalNearInfinity(numpoints, s);
        else
            DecoratorType::wrapped->EstimateOrbitalNearInfinity(numpoints, s);
    }

    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override
    {
        if(ChangeThisSpinorFunction(a))
            return DecoratorType::ApplyTo(a);
        else
            return DecoratorType::wrapped->ApplyTo(a);
    }

protected:
    virtual bool ChangeThisSpinorFunction(const SpinorFunction& fg) const
    {
        // Example: only shift s-wave
        return (fg.Kappa() == -1);

        // To limit a particular orbital, cast with
        //   const Orbital* = dynamic_cast<const Orbital*>(&fg);
    }
};

#endif
