#ifndef RADIATIVE_POTENTIAL_H
#define RADIATIVE_POTENTIAL_H

#include "HartreeFock/HFOperator.h"
#include "HartreeFock/LocalPotentialDecorator.h"
#include "LimitApplicabilityWrapper.h"

/** Uehling potential decorator.
    Uses step potential from Ginges & Berengut.
 */
class UehlingDecorator: public HFOperatorDecorator<LocalPotentialDecorator, UehlingDecorator>
{
public:
    /** Generate Uehling potential using step potential with radius nuclear_rms_radius (fm). */
    UehlingDecorator(pHFOperator wrapped_hf, double nuclear_rms_radius, pIntegrator integration_strategy = pIntegrator());
    UehlingDecorator(pHFOperator wrapped_hf, const RadialFunction& density, pIntegrator integration_strategy = pIntegrator());

    /** Set local direct_potential to Uehling potential with nuclear_rms_radius (fm). */
    void GenerateStepUehling(double nuclear_rms_radius);

    /** Generate Uehling from given density. */
    void GenerateUehling(const RadialFunction& density);
};

/** The Flambaum-Ginges radiative potentials are only applicable to s- and p- waves.
    They are not fitted to d-waves, and vastly overestimate their effect.
    We use a LimitApplicabilityWrapper to achieve this.
 */
template <class DecoratorType>
class LimitToSandP : public HFOperatorDecorator<LimitApplicabilityWrapper<DecoratorType>, LimitToSandP<DecoratorType>>
{
public:
    template <class... DecoratorTypeArgs>
    LimitToSandP(DecoratorTypeArgs... args):
        HFOperatorDecorator<LimitApplicabilityWrapper<DecoratorType>, LimitToSandP<DecoratorType>>(args...)
    {}

protected:
    virtual bool ChangeThisSpinorFunction(const SpinorFunction& fg) const override
    {
        return (fg.L() < 2);
    }
};

/** Radiative potential decorator for magnetic part of self-energy.
    Includes potential for step density and variable density.
 */
class MagneticSelfEnergyDecorator: public HFOperatorDecorator<HFBasicDecorator, MagneticSelfEnergyDecorator>
{
public:
    MagneticSelfEnergyDecorator(pHFOperator wrapped_hf, double nuclear_rms_radius, pIntegrator integration_strategy = pIntegrator());
    MagneticSelfEnergyDecorator(pHFOperator wrapped_hf, const RadialFunction& density, pIntegrator integration_strategy = pIntegrator());

public:
    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const override;
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const override;
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const override;

    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;

public:
    /** Generate magnetic part of self-energy with nuclear_rms_radius (fm). */
    void GenerateStepMagnetic(double nuclear_rms_radius);

    /** Generate magnetic part from given density. */
    void GenerateMagnetic(const RadialFunction& density);

protected:
    /** Magnetic potential is off-diagonal in upper and lower Dirac components. */
    RadialFunction magnetic;
};

/** Radiative potential decorator for electric part of self-energy (high and low frequency).
    Includes potential for step density and variable density.
 */
class ElectricSelfEnergyDecorator: public HFOperatorDecorator<LimitToSandP<LocalPotentialDecorator>, ElectricSelfEnergyDecorator>
{
public:
    ElectricSelfEnergyDecorator(pHFOperator wrapped_hf, double nuclear_rms_radius, bool integrate_offmass_term = true, pIntegrator integration_strategy = pIntegrator());
    ElectricSelfEnergyDecorator(pHFOperator wrapped_hf, const RadialFunction& density, pIntegrator integration_strategy = pIntegrator());

    /** Generate high-frequency electric part of self-energy with nuclear_rms_radius (fm). */
    void GenerateStepEhigh(double nuclear_rms_radius);

    /** Generate high-frequency electric part from given density. */
    void GenerateEhigh(const RadialFunction& density);

    /** Generate low-frequency electric part of self-energy with nuclear_rms_radius (fm). */
    void GenerateStepElow(double nuclear_rms_radius);

    /** Generate low-frequency electric part from given density. */
    void GenerateElow(const RadialFunction& density);

protected:
    /** Initialise Afit and Ra. */
    void Initialize();

protected:
    double Afit;    //!< Fitting function for Ehigh term.
    double Bfit;    //!< Fitting function for Elow term.
    double Ra;      //!< Cutoff radius for offmass term.
    bool integrate_offmass_term;    //!< Integrate offmass term in GenerateStepEhigh.
};

typedef std::shared_ptr<UehlingDecorator> pUehlingDecorator;
typedef std::shared_ptr<MagneticSelfEnergyDecorator> pMagneticSelfEnergyDecorator;
typedef std::shared_ptr<ElectricSelfEnergyDecorator> pElectricSelfEnergyDecorator;

#endif
