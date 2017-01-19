#ifndef RADIATIVE_POTENTIAL_H
#define RADIATIVE_POTENTIAL_H

#include "HartreeFock/HFOperator.h"
#include "HartreeFock/LocalPotentialDecorator.h"
#include "LimitApplicabilityWrapper.h"
#include "Transitions.h"

/** Uehling potential decorator.
    Uses step potential from Ginges & Berengut, J. Phys. B 49, 095001 (2016).
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

/** Radiative potential decorator for magnetic part of self-energy.
    Includes potential for step density and variable density.
    [Ginges & Berengut, Phys. Rev. A 93, 052509 (2016)]
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
    Implementation: s and p-wave potential is stored in directPotential (following LocalPotentialDecorator),
                    d-wave potential is stored in potDWave.
    [Ginges & Berengut, Phys. Rev. A 93, 052509 (2016)]
 */
class ElectricSelfEnergyDecorator: public HFOperatorDecorator<HFBasicDecorator, ElectricSelfEnergyDecorator>
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

public:
    virtual RadialFunction GetDirectPotential() const override;
    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const override;
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const override;
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const override;
    virtual void EstimateOrbitalNearOrigin(unsigned int numpoints, SpinorFunction& s) const override;

public:
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;

protected:
    /** Initialise Afit and Ra. */
    void Initialize();

protected:
    double AfitSP;  //!< Fitting function for Ehigh term for s and p-waves.
    double BfitSP;  //!< Fitting function for Elow term for s and p-waves.
    double BfitD;   //!< Fitting function for Elow term for d-waves.
    double Ra;      //!< Cutoff radius for offmass term.
    bool integrate_offmass_term;    //!< Integrate offmass term in GenerateStepEhigh.

    RadialFunction directPotential; //!< Potential for s and p waves
    RadialFunction potDWave;        //!< Potential for d waves
};

typedef std::shared_ptr<UehlingDecorator> pUehlingDecorator;
typedef std::shared_ptr<MagneticSelfEnergyDecorator> pMagneticSelfEnergyDecorator;
typedef std::shared_ptr<ElectricSelfEnergyDecorator> pElectricSelfEnergyDecorator;

class QEDCalculator : public TransitionCalculator
{
public:
    QEDCalculator(MultirunOptions& user_input, Atom& atom);

    virtual void PrintHeader() const override;
    virtual void PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const override;
};

#endif
