#include "Include.h"
#include "Core.h"
#include "Universal/Constant.h"
#include "Universal/CoulombIntegrator.h"

Core::Core(Lattice* lat, unsigned int atomic_number, int ion_charge):
    StateManager(lat, atomic_number, ion_charge),
    NuclearRadius(0.00001), NuclearThickness(0.000001),
    NuclearInverseMass(0.0), VolumeShiftParameter(0.0), Polarisability(0.0)
{}

void Core::Initialise()
{
    Clear();
    BuildFirstApproximation();

    // Some default values for the nuclear parameters.
    UpdateNuclearPotential();
    Update();
}

void Core::Write(FILE* fp) const
{
    fwrite(&NuclearRadius, sizeof(double), 1, fp);
    fwrite(&NuclearThickness, sizeof(double), 1, fp);

    unsigned int num_states = NumStates();

    // Output core
    fwrite(&num_states, sizeof(unsigned int), 1, fp);
    ConstStateIterator it = GetConstStateIterator();
    while(!it.AtEnd())
    {
        it.GetState()->Write(fp);
        it.Next();
    }
}

void Core::Read(FILE* fp)
{
    Clear();
    
    unsigned int num_core, i;

    fread(&NuclearRadius, sizeof(double), 1, fp);
    fread(&NuclearThickness, sizeof(double), 1, fp);

    // Read core
    fread(&num_core, sizeof(unsigned int), 1, fp);
    for(i = 0; i<num_core; i++)
    {
        DiscreteState* ds = new DiscreteState(lattice);
        ds->Read(fp);
        AddState(ds);
    }

    UpdateNuclearPotential();
    Update();
}

void Core::CalculateVolumeShiftPotential(double radius_difference)
{
    std::vector<double> new_density = CalculateNuclearDensity(NuclearRadius + radius_difference, NuclearThickness);
    std::vector<double> new_potential;
    CoulombIntegrator I(*lattice);
    I.CoulombIntegrate(new_density, new_potential, 0, Z);

    unsigned int oldsize = NuclearPotential.size(),
                 newsize = new_potential.size();
    unsigned int i = 0;
    VolumeShiftPotential.resize(mmax(newsize, oldsize));
    while(i < mmin(newsize, oldsize))
    {   VolumeShiftPotential[i] = (new_potential[i] - NuclearPotential[i]);
        i++;
    }
    while(i < newsize)
    {   VolumeShiftPotential[i] = new_potential[i];
        i++;
    }
    while(i < oldsize)
    {   VolumeShiftPotential[i] = -NuclearPotential[i];
        i++;
    }

    if(VolumeShiftParameter)
        Update();
}

std::vector<double> Core::GetHFPotential() const
{
    return HFPotential;
}

std::vector<double> Core::GetLocalExchangeApproximation() const
{
    return LocalExchangeApproximation;
}

unsigned int Core::UpdateExcitedState(State* s, const SigmaPotential* sigma, double sigma_amount) const
{
    DiscreteState* ds = dynamic_cast<DiscreteState*>(s);
    if(ds != NULL)
        return CalculateDiscreteState(ds, 1., sigma, sigma_amount);
    else
    {   ContinuumState* cs = dynamic_cast<ContinuumState*>(s);
        return CalculateContinuumState(cs);
    }
}
