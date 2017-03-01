#include "YukawaPotential.h"

YukawaDecorator::YukawaDecorator(pHFOperator wrapped_hf, double mass, double scale, pIntegrator integration_strategy):
    BaseDecorator(wrapped_hf, integration_strategy), mass(mass)
{
    SetScale(scale);
    GenerateYukawaPotential(mass);
}

void YukawaDecorator::GenerateYukawaPotential(double new_mass)
{
    mass = new_mass;
    double r_c = physicalConstant->GetAlpha()/mass;
    directPotential.resize(lattice->size());
    auto R = lattice->R();

    for(unsigned int i = 0; i < directPotential.size(); i++)
    {
        directPotential.f[i] = exp(-R[i]/r_c)/R[i];
        directPotential.dfdr[i] = directPotential.f[i] * (-1./R[i] - 1/r_c);
    }
}

void YukawaDecorator::Alert()
{
    GenerateYukawaPotential(mass);
}

