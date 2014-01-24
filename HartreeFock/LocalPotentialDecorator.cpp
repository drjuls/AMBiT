#include "LocalPotentialDecorator.h"
#include "Include.h"
#include "Universal/PhysicalConstant.h"
#include "CoulombOperator.h"
#include "Universal/Interpolator.h"

LocalPotentialDecorator::LocalPotentialDecorator(OneBodyOperator* wrapped_OBO, SpinorODE* wrapped_ODE, pOPIntegrator integration_strategy):
    OneBodyOperatorDecorator(wrapped_OBO, integration_strategy), SpinorODEDecorator(wrapped_ODE)
{}

void LocalPotentialDecorator::GetODEFunction(unsigned int latticepoint, const SpinorFunction& fg, double* w) const
{
    wrapped->GetODEFunction(latticepoint, fg, w);

    if(latticepoint < extraLocalPotential.Size())
    {   const double alpha = PhysicalConstant::Instance()->GetAlpha();
        w[0] += alpha * extraLocalPotential.f[latticepoint] * fg.g[latticepoint];
        w[1] -= alpha * extraLocalPotential.f[latticepoint] * fg.f[latticepoint];
    }
}

void LocalPotentialDecorator::GetODECoefficients(unsigned int latticepoint, const SpinorFunction& fg, double* w_f, double* w_g, double* w_const) const
{
    wrapped->GetODECoefficients(latticepoint, fg, w_f, w_g, w_const);

    if(latticepoint < extraLocalPotential.Size())
    {   const double alpha = PhysicalConstant::Instance()->GetAlpha();
        w_g[0] += alpha * extraLocalPotential.f[latticepoint];
        w_f[1] -= alpha * extraLocalPotential.f[latticepoint];
    }
}

void LocalPotentialDecorator::GetODEJacobian(unsigned int latticepoint, const SpinorFunction& fg, double** jacobian, double* dwdr) const
{
    wrapped->GetODEJacobian(latticepoint, fg, jacobian, dwdr);

    if(latticepoint < extraLocalPotential.Size())
    {   const double alpha = PhysicalConstant::Instance()->GetAlpha();
        jacobian[0][1] += alpha * extraLocalPotential.f[latticepoint];
        jacobian[1][0] -= alpha * extraLocalPotential.f[latticepoint];

        dwdr[0] += alpha * extraLocalPotential.dfdr[latticepoint] * fg.g[latticepoint];
        dwdr[1] -= alpha * extraLocalPotential.dfdr[latticepoint] * fg.f[latticepoint];
    }
}

SpinorFunction LocalPotentialDecorator::ApplyTo(const SpinorFunction& a) const
{
    SpinorFunction ta = component->ApplyTo(a);
    ta -= a * extraLocalPotential;
    
    return ta;
}

LocalExchangeApproximation::LocalExchangeApproximation(OneBodyOperator* wrapped_OBO, SpinorODE* wrapped_ODE):
    LocalPotentialDecorator(wrapped_OBO, wrapped_ODE)
{}

void LocalExchangeApproximation::SetCore(const Core* hf_core)
{
    SpinorODEDecorator::SetCore(hf_core);
    const double* R = lattice->R();
    
    // Get electron density function
    RadialFunction density;
    
    ConstStateIterator cs = core->GetConstStateIterator();
    while(!cs.AtEnd())
    {
        const Orbital& s = *cs.GetState();
        double number_electrons = s.Occupancy();
        
        density += s.GetDensity() * number_electrons;
        
        cs.Next();
    }

    RadialFunction y(density.Size());

    if(density.Size())
    {   CoulombOperator coulombSolver(lattice);
        coulombSolver.GetPotential(0, density, y);
    }

    unsigned int i = 0;
    const double Z = core->GetZ();
    const double Charge = core->GetCharge();

    // Get local exchange approximation
    extraLocalPotential.ReSize(density.Size());
    double C = 0.635348143228;
    for(i = 0; i < density.Size(); i++)
    {
        extraLocalPotential.f[i] = C * pow((density.f[i]/(R[i]*R[i])), 1./3.);
        extraLocalPotential.dfdr[i] = C/3. * pow((density.f[i]/(R[i]*R[i])), -2./3.) * (density.dfdr[i] - 2.*density.f[i]/R[i])/(R[i]*R[i]);

        if(extraLocalPotential.f[i] > (Z - Charge)/R[i] - y.f[i])
            break;
    }

    while(i < density.Size())
    {   extraLocalPotential.f[i] = (Z - Charge)/R[i] - y.f[i];
        extraLocalPotential.dfdr[i] = -(Z - Charge)/(R[i]*R[i]) - y.dfdr[i];
        i++;
    }
}

