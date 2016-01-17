#include "NucleusDecorator.h"
#include "Universal/MathConstant.h"

NucleusDecorator::NucleusDecorator(pHFOperator wrapped_hf, pIntegrator integration_strategy):
    BaseDecorator(wrapped_hf, integration_strategy), nuclear_radius(0.0), nuclear_thickness(0.0)
{}

void NucleusDecorator::SetFermiParameters(double radius_fm, double thickness_fm)
{
    nuclear_radius = radius_fm;
    nuclear_thickness = thickness_fm;

    directPotential.Clear();
    const double* R = lattice->R();

    if(nuclear_radius > 0.1)
    {   RadialFunction nuclear_density = CalculateNuclearDensity(nuclear_radius, nuclear_thickness);
        coulombSolver->GetPotential(nuclear_density, directPotential, Z);

        // Remove point-nucleus potential
        unsigned int i = 0;
        while(i < directPotential.size())
        {   directPotential.f[i] -= Z/R[i];
            directPotential.dfdr[i] += Z/(R[i]*R[i]);
            i++;
        }
    }
}

double NucleusDecorator::CalculateNuclearRMSRadius() const
{
    if(nuclear_radius <= 0.1)
        return 0.0;
    else if(nuclear_thickness <= 0.01)
        return std::sqrt(3./5.) * nuclear_radius;

    // nuclear_density(r) = r^2 rho(r)
    RadialFunction nuclear_density = CalculateNuclearDensity(nuclear_radius, nuclear_thickness);
    const double* R = lattice->R();
    const double* R2 = lattice->Rpower(2);

    double r2rho = integrator->Integrate(nuclear_density);

    for(int i = 0; i < nuclear_density.size(); i++)
    {
        nuclear_density.dfdr[i] = R2[i] * nuclear_density.dfdr[i] + 2. * R[i] * nuclear_density.f[i];
        nuclear_density.f[i] *= R2[i];
    }

    double r4rho = integrator->Integrate(nuclear_density);

    return std::sqrt(r4rho/r2rho) * MathConstant::Instance()->BohrRadiusInFermi();
}

RadialFunction NucleusDecorator::GetNuclearDensity() const
{
    return CalculateNuclearDensity(nuclear_radius, nuclear_thickness);
}

RadialFunction NucleusDecorator::CalculateNuclearDensity(double radius, double thickness) const
{
    RadialFunction density(lattice->size());

    const double fermi_length = 1.0/MathConstant::Instance()->BohrRadiusInFermi();
    const double pi = MathConstant::Instance()->Pi();

    radius *= fermi_length;
    thickness *= fermi_length;

    const double* R = lattice->R();

    if(radius < R[1])
        density.resize(0);
    else if(thickness > 0.01 * fermi_length)
    {
        double B = 4.*log(3.)/thickness;
        double A = 3.*Z/(radius * (radius*radius + (pi*pi)/(B*B))); // Approximate normalisation

        for(unsigned int i=0; i<lattice->size(); i++)
        {
            double X = B * (R[i] - radius);
            if(X <= -20.)
            {   density.f[i] = A * R[i] * R[i];
                density.dfdr[i] = 2.0 * A * R[i];
            }
            else if(X < 50.)
            {   double exp_X = exp(X);
                density.f[i] = A/(1. + exp(X)) * R[i] * R[i];
                density.dfdr[i] = A * 2. * R[i] * (1. + (1. - 0.5 * B * R[i]) * exp_X)
                                                 / ((1. + exp_X) * (1. + exp_X));
            }
            else
            {   density.resize(i);
                break;
            }
        }
    }
    else
    {   // zero thickness
        double A = 3.*Z/(radius * radius * radius);
        unsigned int i;
        for(i = 0; i < lattice->real_to_lattice(radius); i++)
        {   density.f[i] = A * R[i] * R[i];
            density.dfdr[i] = 2.0 * A * R[i];
        }
        density.resize(i);
    }

    // Renormalise
    double r2rho = integrator->Integrate(density);

    return density * (Z/r2rho);
}
