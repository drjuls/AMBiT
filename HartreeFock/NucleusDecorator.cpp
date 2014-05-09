#include "NucleusDecorator.h"
#include "Universal/MathConstant.h"

NucleusDecorator::NucleusDecorator(pHFOperator wrapped_hf, pOPIntegrator integration_strategy):
    LocalPotentialDecorator(wrapped_hf, integration_strategy), nuclear_radius(0.0), nuclear_thickness(0.0)
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
            directPotential.dfdr[i] += -Z/(R[i]*R[i]);
            i++;
        }
    }
}

RadialFunction NucleusDecorator::CalculateNuclearDensity(double radius, double thickness) const
{
    RadialFunction density(lattice->Size());

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
        double A = 3.*Z/(radius * (radius*radius + (pi*pi)/(B*B)));
        for(unsigned int i=0; i<lattice->Size(); i++)
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

    return density;
}