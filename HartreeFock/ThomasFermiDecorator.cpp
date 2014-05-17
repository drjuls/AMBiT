#include "ThomasFermiDecorator.h"
#include "Include.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/Interpolator.h"

ThomasFermiDecorator::ThomasFermiDecorator(pHFOperator decorated_object, pOPIntegrator integration_strategy):
    HFOperatorDecorator(decorated_object)
{
    // If integration_strategy is supplied, use it.
    // Otherwise the integration_strategy from wrapped_hf will be used.
    if(integration_strategy != NULL)
        integrator = integration_strategy;
}

void ThomasFermiDecorator::SetCore(pCoreConst hf_core, double hf_mixing)
{
    wrapped->SetCore(hf_core);
    core = hf_core;

    unsigned int i;
    const double* R = lattice->R();
    double C = mmax(charge, 1.);

    if(hf_mixing == 0.0)
    {
        // Thomas-Fermi potential.
        directPotential.resize(lattice->Size());
        double P = 2.235 * pow((Z/55.), 1./3.);

        for(i = 0; i < directPotential.size(); i++)
        {
            double r = R[i];
            double denom = 1. + P * r;
            denom = denom * denom * (1. + exp((r - 3.1)/0.4));
            double Zeff = (Z - C) / denom + C;

            directPotential.f[i] = Zeff/r;
        }
    }
    else
    {   RadialFunction hf_pot = wrapped->GetDirectPotential();
        directPotential = directPotential * (1. - hf_mixing) + hf_pot * hf_mixing;

        // potential should be at least C/r everywhere.
        for(i = 0; i < directPotential.size(); i++)
            directPotential.f[i] = mmax(directPotential.f[i], C/R[i]);
    }

    Interpolator I(lattice);
    I.GetDerivative(directPotential.f, directPotential.dfdr, 6);
}

RadialFunction ThomasFermiDecorator::GetDirectPotential() const
{
    return directPotential;
}

void ThomasFermiDecorator::ExtendPotential()
{
    wrapped->ExtendPotential();

    unsigned int i = size();
    const double* R = lattice->R();
    double C = mmax(charge, 1.);

    directPotential.resize(lattice->Size());

    while(i < directPotential.size())
    {   directPotential.f[i] = C/R[i];
        directPotential.dfdr[i] = -C/(R[i]*R[i]);
        i++;
    }
}

SpinorFunction ThomasFermiDecorator::GetExchange(pSingleParticleWavefunctionConst approximation) const
{
    int kappa = -1;
    if(approximation)
        kappa = approximation->Kappa();

    return SpinorFunction(kappa);
}

void ThomasFermiDecorator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    const double R = lattice->R(latticepoint);
    const double alpha = physicalConstant->GetAlpha();

    const unsigned int i = latticepoint;

    w[0] = -currentKappa/R * fg.f[i] + (2./alpha + alpha * (currentEnergy + directPotential.f[i])) * fg.g[i];
    w[1] = -alpha * (currentEnergy + directPotential.f[i]) * fg.f[i] + currentKappa/R * fg.g[i];
}

void ThomasFermiDecorator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    const double R = lattice->R(latticepoint);
    const double alpha = physicalConstant->GetAlpha();

    const unsigned int i = latticepoint;

    w_f[0] = -currentKappa/R;
    w_g[0] = 2./alpha + alpha * (currentEnergy + directPotential.f[i]);
    w_f[1] = -alpha * (currentEnergy + directPotential.f[i]);
    w_g[1] = currentKappa/R;

    w_const[0] = 0.;
    w_const[1] = 0.;
}

void ThomasFermiDecorator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    const double R = lattice->R(latticepoint);
    const double alpha = physicalConstant->GetAlpha();

    const unsigned int i = latticepoint;

    jacobian[0][0] = -currentKappa/R;
    jacobian[0][1] = 2./alpha + alpha * (currentEnergy + directPotential.f[i]);
    jacobian[1][0] = -alpha * (currentEnergy + directPotential.f[i]);
    jacobian[1][1] = currentKappa/R;

    dwdr[0] = currentKappa/(R*R) * fg.f[i] + alpha * directPotential.dfdr[i] * fg.g[i];
    dwdr[1] = -alpha * directPotential.dfdr[i] * fg.f[i] -currentKappa/(R*R) * fg.g[i];
}

SpinorFunction ThomasFermiDecorator::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta(a.Kappa(), a.size());

    const double alpha = physicalConstant->GetAlpha();
    const double alphasquared = physicalConstant->GetAlphaSquared();
    const double* R = lattice->R();

    double kappa = a.Kappa();

    std::vector<double> d2fdr2(a.size());
    std::vector<double> d2gdr2(a.size());
    Interpolator I(lattice);
    I.GetDerivative(a.dfdr, d2fdr2, 6);
    I.GetDerivative(a.dgdr, d2gdr2, 6);

    for(unsigned int i = 0; i < ta.size(); i++)
    {
        ta.f[i] = -directPotential.f[i]*a.f[i]
                  + (-a.dgdr[i] + kappa/R[i]*a.g[i])/alpha;
        ta.g[i] = (a.dfdr[i] + kappa/R[i]*a.f[i])/alpha
                  - (2./alphasquared + directPotential.f[i])*a.g[i];

        ta.dfdr[i] = - directPotential.f[i]*a.dfdr[i] - directPotential.dfdr[i]*a.f[i]
                     + (-d2gdr2[i] + kappa/R[i] * (a.dgdr[i] - a.g[i]/R[i]))/alpha;
        ta.dgdr[i] = (d2fdr2[i] + kappa/R[i] * (a.dfdr[i] - a.f[i]/R[i]))/alpha
                     - (2./alphasquared + directPotential.f[i]) * a.dgdr[i] - directPotential.dfdr[i] * a.g[i];
    }

    return ta;
}
