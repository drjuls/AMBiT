#ifndef LOCAL_POTENTIAL_DECORATOR
#define LOCAL_POTENTIAL_DECORATOR

#include "HFOperator.h"

/** Template class to add an extra local potential to a HF operator.
    Extra local potential should be stored in directPotential (inherited from HF operator) and
    must be set by subclasses via SetCore() and/or SetODEParameters().
    Note sign is for normal electrostatic potential: V_nucleus(r) > 0.
 */
template <class Derived>
class LocalPotentialDecorator : public HFOperatorDecorator<Derived>
{
public:
    LocalPotentialDecorator(pHFOperator wrapped_hf, pOPIntegrator integration_strategy = pOPIntegrator()):
        HFOperatorDecorator<Derived>(wrapped_hf), scale(1.)
    {
        // If integration_strategy is supplied, use it.
        // Otherwise the integration_strategy from wrapped_hf will be used.
        if(integration_strategy != NULL)
            HFOperatorDecorator<Derived>::integrator = integration_strategy;
    }

    LocalPotentialDecorator(const LocalPotentialDecorator<Derived>& other):
        HFOperatorDecorator<Derived>(other), scale(other.scale)
    {}

    void SetScale(double factor) { scale = factor; }
    double GetScale() const { return scale; }

public:
    virtual RadialFunction GetDirectPotential() const override;
    virtual void GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const override;
    virtual void GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const override;
    virtual void GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const override;

public:
    virtual SpinorFunction ApplyTo(const SpinorFunction& a) const override;

protected:
    double scale;
};

/** Import a potential from a file.
    Imported file should be a text list of radial points and corresponding potentials
        R   V(R)
    with one point per line.
 */
class ImportedPotentialDecorator : public LocalPotentialDecorator<ImportedPotentialDecorator>
{
public:
    ImportedPotentialDecorator(pHFOperator wrapped_hf, const std::string& filename, pOPIntegrator integration_strategy = pOPIntegrator());
    ImportedPotentialDecorator(const ImportedPotentialDecorator& other):
        LocalPotentialDecorator(other) {}
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
class LocalExchangeApproximation : public LocalPotentialDecorator<LocalExchangeApproximation>
{
public:
    LocalExchangeApproximation(pHFOperator wrapped_hf, double Xalpha = 1.0, pOPIntegrator integration_strategy = pOPIntegrator());
    LocalExchangeApproximation(const LocalExchangeApproximation& other):
        LocalPotentialDecorator(other), Xalpha(other.Xalpha)
    {}

    void SetXalpha(double x_alpha) { Xalpha = x_alpha; }
    double GetXalpha() const { return Xalpha; }

    virtual void SetCore(pCoreConst hf_core) override;
    virtual void Alert() override;

protected:
    double Xalpha;
};

typedef std::shared_ptr<LocalExchangeApproximation> pLocalExchangeApproximation;
typedef std::shared_ptr<const LocalExchangeApproximation> pLocalExchangeApproximationConst;

// Template functions for LocalPotentialDecorator below

template <class Derived>
RadialFunction LocalPotentialDecorator<Derived>::GetDirectPotential() const
{
    RadialFunction ret = HFOperatorDecorator<Derived>::wrapped->GetDirectPotential();
    ret +=  HFOperatorDecorator<Derived>::directPotential * scale;

    return ret;
}

template <class Derived>
void LocalPotentialDecorator<Derived>::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    HFOperatorDecorator<Derived>::wrapped->GetODEFunction(latticepoint, fg, w);

    if(latticepoint < HFOperatorDecorator<Derived>::directPotential.size())
    {   const double alpha = HFOperatorDecorator<Derived>::physicalConstant->GetAlpha() * scale;
        w[0] += alpha * HFOperatorDecorator<Derived>::directPotential.f[latticepoint] * fg.g[latticepoint];
        w[1] -= alpha * HFOperatorDecorator<Derived>::directPotential.f[latticepoint] * fg.f[latticepoint];
    }
}

template <class Derived>
void LocalPotentialDecorator<Derived>::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    HFOperatorDecorator<Derived>::wrapped->GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);

    if(latticepoint < HFOperatorDecorator<Derived>::directPotential.size())
    {   const double alpha = HFOperatorDecorator<Derived>::physicalConstant->GetAlpha() * scale;
        w_g[0] += alpha * HFOperatorDecorator<Derived>::directPotential.f[latticepoint];
        w_f[1] -= alpha * HFOperatorDecorator<Derived>::directPotential.f[latticepoint];
    }
}

template <class Derived>
void LocalPotentialDecorator<Derived>::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    HFOperatorDecorator<Derived>::wrapped->GetODEJacobian(latticepoint, fg, jacobian, dwdr);

    if(latticepoint < HFOperatorDecorator<Derived>::directPotential.size())
    {   const double alpha = HFOperatorDecorator<Derived>::physicalConstant->GetAlpha() * scale;
        jacobian[0][1] += alpha * HFOperatorDecorator<Derived>::directPotential.f[latticepoint];
        jacobian[1][0] -= alpha * HFOperatorDecorator<Derived>::directPotential.f[latticepoint];

        dwdr[0] += alpha * HFOperatorDecorator<Derived>::directPotential.dfdr[latticepoint] * fg.g[latticepoint];
        dwdr[1] -= alpha * HFOperatorDecorator<Derived>::directPotential.dfdr[latticepoint] * fg.f[latticepoint];
    }
}

template <class Derived>
SpinorFunction LocalPotentialDecorator<Derived>::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta = HFOperatorDecorator<Derived>::wrapped->ApplyTo(a);
    ta -= a * HFOperatorDecorator<Derived>::directPotential * scale;

    return ta;
}

#endif
