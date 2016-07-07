#ifndef LOCAL_POTENTIAL_DECORATOR
#define LOCAL_POTENTIAL_DECORATOR

#include "HFOperator.h"

/** Add an extra local potential to a HF operator.
    Extra local potential should be stored in directPotential and
    must be set by subclasses via SetCore() and/or SetODEParameters().
    Note sign is for normal electrostatic potential: V_nucleus(r) > 0.
 */
class LocalPotentialDecorator : public HFOperatorDecorator<HFBasicDecorator, LocalPotentialDecorator>
{
public:
    LocalPotentialDecorator(pHFOperator wrapped_hf, pIntegrator integration_strategy = nullptr):
        BaseDecorator(wrapped_hf, integration_strategy), scale(1.)
    {}

    void SetScale(double factor) { scale = factor; }
    double GetScale() const { return scale; }

public:
    virtual RadialFunction GetDirectPotential() const override;
    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const override;
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const override;
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const override;
    virtual void EstimateOrbitalNearOrigin(unsigned int numpoints, SpinorFunction& s) const override;

public:
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;

protected:
    double scale;
    RadialFunction directPotential; // Additional direct potential
};

/** Import a potential from a file.
    Imported file should be a text list of radial points and corresponding potentials
        R   V(R)
    with one point per line.
 */
class ImportedPotentialDecorator : public HFOperatorDecorator<LocalPotentialDecorator, ImportedPotentialDecorator>
{
public:
    ImportedPotentialDecorator(pHFOperator wrapped_hf, const std::string& filename, pIntegrator integration_strategy = pIntegrator());
};

typedef std::shared_ptr<ImportedPotentialDecorator> pImportedPotentialDecorator;
typedef std::shared_ptr<const ImportedPotentialDecorator> pImportedPotentialDecoratorConst;

/** Decorate HF operator with a local approximation to the exchange potential. Good for first approximations.
    The Latter correction (Zeff -> 1 for neutral atoms) is enforced [Phys. Rev. 99, 510 (1955)].
    The exchange potential is given by
        \f[ x_{\alpha} \left( \frac{81}{32 \pi^2} \frac{\rho(r)}{r^2} \right)^{1/3} \f]
    where \f$\rho(r)\f$ is the density of electrons. The prefactor
        \f$x_\alpha = 1\f$ corresponds to the Dirac-Slater potential, while
        \f$x_\alpha = 2/3\f$ gives the Kohn-Sham potential and
        \f$x_\alpha = 0\f$ is Core-Hartree.
    Xalpha is not quite the same as scaling, since the Latter correction is not affected.
    Remember to turn off the exchange potential using wrapped_hf->Include.
 */
class LocalExchangeApproximation : public HFOperatorDecorator<LocalPotentialDecorator, LocalExchangeApproximation>
{
public:
    LocalExchangeApproximation(pHFOperator wrapped_hf, pCoulombOperator coulomb, double Xalpha = 1.0, pIntegrator integration_strategy = nullptr);

    void SetXalpha(double x_alpha) { Xalpha = x_alpha; }
    double GetXalpha() const { return Xalpha; }

    virtual void SetCore(pCoreConst hf_core) override;
    virtual void Alert() override;

protected:
    double Xalpha;
    pCoulombOperator coulombSolver;
};

typedef std::shared_ptr<LocalExchangeApproximation> pLocalExchangeApproximation;
typedef std::shared_ptr<const LocalExchangeApproximation> pLocalExchangeApproximationConst;

#endif
