#include "NuclearPolarisability.h"

namespace Ambit
{
NuclearPolarisability::NuclearPolarisability(pHFOperator wrapped_hf, double alphaE, double nuclear_excitation_energy_MeV, pIntegrator integration_strategy):
    BaseDecorator(wrapped_hf, integration_strategy), Ebar(nuclear_excitation_energy_MeV)
{
    scale = alphaE;
    GeneratePotential();
}

void NuclearPolarisability::GeneratePotential()
{
    directPotential.Clear();
    directPotential.resize(lattice->size());

    MathConstant* math = MathConstant::Instance();
    const double* R = lattice->R();
    const double* R3 = lattice->Rpower(3);
    const double* R4 = lattice->Rpower(4);

    double b = physicalConstant->GetAlpha() * math->Pi() * math->Pi()/std::sqrt(2);
    b = b/(19./6. + 5. * std::log(2. * Ebar * 1.e6/math->ElectronMassInEV));
    double b4 = gsl_pow_4(b);

    for(unsigned int i = 0; i < directPotential.size(); i++)
    {
        directPotential.f[i] = -0.5/(R4[i] + b4);
        directPotential.dfdr[i] = 2. * R3[i]/gsl_pow_2(R4[i] + b4);

        if(fabs(directPotential.f[i]) < 1.e-15)
        {   directPotential.resize(i);
            break;
        }
    }
}

}
