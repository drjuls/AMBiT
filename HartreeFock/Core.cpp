#include "Include.h"
#include "Core.h"
#include "Universal/Constant.h"
#include "Universal/CoulombIntegrator.h"
#include "Universal/Interpolator.h"

Core::Core(Lattice* lat, unsigned int atomic_number, int ion_charge):
    StateManager(lat, atomic_number, ion_charge),
    NuclearRadius(0.0), NuclearThickness(0.0),
    NuclearInverseMass(0.0), VolumeShiftParameter(0.0), Polarisability(0.0), ClosedShellRadius(0.)
{}

Core::Core(const Core& other, Lattice* new_lattice):
    StateManager(other, new_lattice),
    NuclearRadius(other.NuclearRadius), NuclearThickness(other.NuclearThickness),
    NuclearInverseMass(other.NuclearInverseMass), VolumeShiftParameter(other.VolumeShiftParameter),
    Polarisability(other.Polarisability), ClosedShellRadius(other.ClosedShellRadius)
{
    // The StateManager part has copied/interpolated all the wavefunctions.
    UpdateNuclearPotential();
    UpdateHFPotential();

    // Copy original occupancies
    OpenShellStates = other.OpenShellStates;

    // Copy/interpolate any currently unused open shell states
    bool interpolate = !(*lattice == *other.lattice);
    Interpolator interp(other.lattice);
    unsigned int order = 6;

    const double* R_old = other.lattice->R();
    const double* R = lattice->R();
    const double* dR = lattice->dR();

    StateSet::const_iterator it = other.OpenShellStorage.begin();

    while(it != other.OpenShellStorage.end())
    {
        const StatePointer ds_old = it->second;

        // Copy kappa, pqn, etc.
        DiscreteState* ds = new DiscreteState(*ds_old.GetState());

        if(interpolate)
        {
            unsigned int new_size = lattice->real_to_lattice(R_old[ds_old->Size() - 1]);
            double dfdr, dgdr;

            ds->ReSize(new_size);
            for(unsigned int i = 0; i < new_size; i++)
            {
                interp.Interpolate(ds_old->f, R[i], ds->f[i], dfdr, order);
                interp.Interpolate(ds_old->g, R[i], ds->g[i], dgdr, order);
                ds->df[i] = dfdr * dR[i];
                ds->dg[i] = dgdr * dR[i];
            }
        }

        OpenShellStorage[it->first] = ds;
        it++;
    }
}

void Core::Initialise(std::string configuration)
{
    Clear();
    BuildFirstApproximation(configuration);

    // Some default values for the nuclear parameters.
    UpdateNuclearPotential();
    Update();
}

void Core::Clear()
{
    StateManager::Clear();

    // Delete stored states
    StateSet::iterator it = OpenShellStorage.begin();
    while(it != OpenShellStorage.end())
    {   it->second.DeleteState();
        it++;
    }
    OpenShellStorage.clear();

    OpenShellStates.clear();
}

void Core::Write(FILE* fp) const
{
    fwrite(&NuclearRadius, sizeof(double), 1, fp);
    fwrite(&NuclearThickness, sizeof(double), 1, fp);
    fwrite(&NuclearInverseMass, sizeof(double), 1, fp);
    fwrite(&Polarisability, sizeof(double), 1, fp);

    // Output core
    StateManager::Write(fp);

    // Output open shell states
    unsigned int num_open = OpenShellStorage.size();
    fwrite(&num_open, sizeof(unsigned int), 1, fp);

    StateSet::const_iterator it = OpenShellStorage.begin();
    while(it != OpenShellStorage.end())
    {
        it->second->Write(fp);
        it++;
    }

    // Output original occupancies
    num_open = OpenShellStates.size();
    fwrite(&num_open, sizeof(unsigned int), 1, fp);

    std::map<StateInfo, double>::const_iterator info_it = OpenShellStates.begin();
    while(info_it != OpenShellStates.end())
    {
        unsigned int pqn = info_it->first.PQN();
        int kappa = info_it->first.Kappa();
        double occ = info_it->second;

        fwrite(&pqn, sizeof(unsigned int), 1, fp);
        fwrite(&kappa, sizeof(int), 1, fp);
        fwrite(&occ, sizeof(double), 1, fp);

        info_it++;
    }
}

void Core::Read(FILE* fp)
{
    Clear();    
    unsigned int num_core, i;

    fread(&NuclearRadius, sizeof(double), 1, fp);
    fread(&NuclearThickness, sizeof(double), 1, fp);
    fread(&NuclearInverseMass, sizeof(double), 1, fp);
    fread(&Polarisability, sizeof(double), 1, fp);

    // Read core
    fread(&num_core, sizeof(unsigned int), 1, fp);
    for(i = 0; i<num_core; i++)
    {
        DiscreteState* ds = new DiscreteState();
        ds->Read(fp);
        AddState(ds);
    }

    // Read open shell states
    unsigned int num_open;
    fread(&num_open, sizeof(unsigned int), 1, fp);
    for(i = 0; i<num_open; i++)
    {
        DiscreteState* ds = new DiscreteState();
        ds->Read(fp);

        StateInfo info(ds);
        StatePointer sp(ds);
        OpenShellStorage.insert(StateSet::value_type(info, sp));
    }

    // Read original occupancies
    fread(&num_open, sizeof(unsigned int), 1, fp);
    for(i = 0; i<num_open; i++)
    {
        unsigned int pqn;
        int kappa;
        double occ;

        fread(&pqn, sizeof(unsigned int), 1, fp);
        fread(&kappa, sizeof(int), 1, fp);
        fread(&occ, sizeof(double), 1, fp);

        StateInfo info(pqn, kappa);
        OpenShellStates.insert(std::pair<StateInfo, double>(info, occ));
    }

    UpdateNuclearPotential();
    if(Polarisability)
        CalculateClosedShellRadius();
}

void Core::CalculateVolumeShiftPotential(double radius_difference)
{
    // Calculate new potential
    double old_radius = NuclearRadius;
    NuclearRadius = NuclearRadius + radius_difference;
    UpdateNuclearPotential();
    std::vector<double> new_potential(NuclearPotential);

    // Return to old potential
    NuclearRadius = old_radius;
    UpdateNuclearPotential();

    unsigned int oldsize = NuclearPotential.size(),
                 newsize = new_potential.size();
    unsigned int i = 0;
    VolumeShiftPotential.resize(mmax(newsize, oldsize));
    while(i < mmin(newsize, oldsize))
    {   VolumeShiftPotential[i] = (new_potential[i] - NuclearPotential[i]);
        i++;
    }
    while(i < newsize)
    {   VolumeShiftPotential[i] = new_potential[i] - Z/lattice->R(i);
        i++;
    }
    while(i < oldsize)
    {   VolumeShiftPotential[i] = Z/lattice->R(i) - NuclearPotential[i];
        i++;
    }
}

std::vector<double> Core::GetHFPotential() const
{
    return HFPotential;
}

const std::vector<double>& Core::GetConstHFPotential() const
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
    {
        // Number of iterations required. Zero shows that the state existed previously.
        unsigned int loop = 0;

        const DiscreteState* core_state = GetState(StateInfo(ds));
        StateSet::const_iterator it = OpenShellStorage.find(StateInfo(ds));
        if(core_state != NULL)
        {   // Try to find in core (probably open shells).
            *ds = *core_state;
        }
        else if(it != OpenShellStorage.end())
        {   // Try to find in unoccupied (but previously calculated) open shells.
            *ds = *(it->second.GetState());
        }
        else
        {   // Hartree-Fock loops
            bool debugHF = DebugOptions.LogHFIterations();

            double deltaE;
            DiscreteState* new_ds = new DiscreteState(*ds);
            double prop_new = 0.5;
            do
            {   loop++;
                CoupledFunction exchange;
                CalculateExchange(*new_ds, exchange, sigma, sigma_amount);
                deltaE = IterateDiscreteStateGreens(new_ds, &exchange);

                if(debugHF)
                    *logstream << "  " << std::setw(4) << ds->Name() 
                            << "  E = " << std::setprecision(12) << ds->Energy()
                            << "  deltaE = " << std::setprecision(3) << deltaE
                            << "  size: (" << new_ds->Size()
                            << ") " << lattice->R(new_ds->Size()) << std::endl;

                ds->Scale(1. - prop_new);
                new_ds->Scale(prop_new);
                ds->ReSize(mmax(ds->Size(), new_ds->Size()));
                new_ds->ReSize(mmax(ds->Size(), new_ds->Size()));

                for(unsigned int i = 0; i<ds->Size(); i++)
                {   ds->f[i] += new_ds->f[i];
                    ds->g[i] += new_ds->g[i];
                    ds->df[i] += new_ds->df[i];
                    ds->dg[i] += new_ds->dg[i];
                }

                // Renormalise core states (should be close already) and update energy.
                ds->ReNormalise(lattice);
                ds->SetEnergy((1. - prop_new) * ds->Energy() + prop_new * new_ds->Energy());
                *new_ds = *ds;
                deltaE = fabs(deltaE/new_ds->Energy());

            }while((deltaE > StateParameters::EnergyTolerance) && (loop < StateParameters::MaxHFIterations));

            delete new_ds;

            if(loop >= StateParameters::MaxHFIterations)
                *errstream << "Core: Failed to converge excited HF state " << ds->Name() << std::endl;

            if(DebugOptions.OutputHFExcited())
            {
                *logstream << std::setprecision(12);
                if(DebugOptions.HartreeEnergyUnits() || DebugOptions.InvCmEnergyUnits())
                {
                    double energy = ds->Energy();
                    if(DebugOptions.InvCmEnergyUnits())
                        energy *= Constant::HartreeEnergy_cm;
                    *logstream << ds->Name() << "  E = " << energy << "  loops: " << loop << "  size: " << ds->Size() << std::endl;
                }
                else
                    *logstream << ds->Name() << "  nu = " << ds->Nu() << "  loops: " << loop << "  size: " << ds->Size() << std::endl;
            }
        }
        return loop;
    }
    else
    {   ContinuumState* cs = dynamic_cast<ContinuumState*>(s);
        return CalculateContinuumState(cs);
    }
}
