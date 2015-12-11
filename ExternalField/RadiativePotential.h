#ifndef RADIATIVE_POTENTIAL_H
#define RADIATIVE_POTENTIAL_H

#include "HartreeFock/HFOperator.h"
#include "HartreeFock/LocalPotentialDecorator.h"

/** Uehling potential decorator.
    Uses step potential from Ginges & Berengut.
 */
class UehlingDecorator: public LocalPotentialDecorator<UehlingDecorator>
{
public:
    /** Generate Uehling potential using step potential with radius nuclear_rms_radius (fm). */
    UehlingDecorator(pHFOperator wrapped_hf, double nuclear_rms_radius, pOPIntegrator integration_strategy = pOPIntegrator());
    UehlingDecorator(pHFOperator wrapped_hf, const RadialFunction& density, pOPIntegrator integration_strategy = pOPIntegrator());

protected:
    /** Set local direct_potential to Uehling potential with nuclear_rms_radius (fm). */
    void GenerateStepUehling(double nuclear_rms_radius);

    /** Generate Uehling from given density. */
    void GenerateUehling(const RadialFunction& density);
};

/** Radiative potential decorator for magnetic part of self-energy.
    Includes potential for step density and variable density.
 */
class MagneticSelfEnergyDecorator: public HFOperatorDecorator<MagneticSelfEnergyDecorator>
{
public:
    MagneticSelfEnergyDecorator(pHFOperator wrapped_hf, double nuclear_rms_radius, pOPIntegrator integration_strategy = pOPIntegrator());
    MagneticSelfEnergyDecorator(pHFOperator wrapped_hf, const RadialFunction& density, pOPIntegrator integration_strategy = pOPIntegrator());

public:
    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const override;
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const override;
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const override;

public:
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;

protected:
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
class ElectricSelfEnergyDecorator: public LocalPotentialDecorator<ElectricSelfEnergyDecorator>
{
public:
    ElectricSelfEnergyDecorator(pHFOperator wrapped_hf, double nuclear_rms_radius, bool integrate_offmass_term = true, pOPIntegrator integration_strategy = pOPIntegrator());
    ElectricSelfEnergyDecorator(pHFOperator wrapped_hf, const RadialFunction& density, pOPIntegrator integration_strategy = pOPIntegrator());

protected:
    /** Generate high-frequency electric part of self-energy with nuclear_rms_radius (fm). */
    void GenerateStepEhigh(double nuclear_rms_radius);

    /** Generate high-frequency electric part from given density. */
    void GenerateEhigh(const RadialFunction& density);

    /** Generate low-frequency electric part of self-energy with nuclear_rms_radius (fm). */
    void GenerateStepElow(double nuclear_rms_radius);

    /** Generate low-frequency electric part from given density. */
    void GenerateElow(const RadialFunction& density);

    /** Initialise Afit and Ra. */
    void Initialize();

protected:
    double Afit;    //!< Fitting function for Ehigh term.
    double Bfit;    //!< Fitting function for Elow term.
    double Ra;      //!< Cutoff radius for offmass term.
    bool integrate_offmass_term;    //!< Integrate offmass term in GenerateStepEhigh.
};

#endif
