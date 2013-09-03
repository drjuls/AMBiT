#include "Include.h"
#include "ContinuumBuilder.h"
#include "Universal/ExpLattice.h"
#include "StateIntegrator.h"
#include "Universal/Interpolator.h"
#include "Atom/Debug.h"
#include "Universal/PhysicalConstant.h"

ContinuumBuilder::ContinuumBuilder(const Core* other_core):
    lattice(NULL), core(NULL), norm_type(Cowan)
{
    CopyCore(other_core, true);
}

ContinuumBuilder::~ContinuumBuilder()
{
    if(lattice)
        delete lattice;
    if(core)
        delete core;
}

void ContinuumBuilder::CopyLattice(const Lattice* lat)
{
    if(lattice)
        delete lattice;

    const ExpLattice* explat = dynamic_cast<const ExpLattice*>(lat);
    if(explat != NULL)
        lattice = new ExpLattice(*explat);
    else
        lattice = new Lattice(*lat);
}

void ContinuumBuilder::CreateNewLattice(unsigned int numpoints, double r_min, double r_max)
{
    if(lattice)
        delete lattice;
    
    lattice = new Lattice(numpoints, r_min, r_max);
}

void ContinuumBuilder::CopyCore(const Core* other_core, bool import_lattice)
{
    if(core)
        delete core;

    if(import_lattice || !lattice)
    {   CopyLattice(other_core->GetLattice());
    }

    core = new Core(*other_core, lattice);
}

unsigned int ContinuumBuilder::CalculateContinuumWave(ContinuumWave* s, Lattice* external_lattice) const
{
    const std::vector<double>& HFPotential = core->GetConstHFPotential();

    unsigned int loop = 0;
    double final_amplitude, final_phase;
    SpinorFunction exchange(s->Kappa(), HFPotential.size());
    SpinorFunction new_exchange(s->Kappa(), HFPotential.size());
    double ds, old_phase = 0.;
    unsigned int start_sine = 0;

    s->ReSize(HFPotential.size());
    unsigned int lattice_extensions = 0;

    do
    {   loop++;
        
        StateIntegrator I(lattice);
        start_sine = I.IntegrateContinuum(*s, HFPotential, exchange, core->GetZ(), 0.01, final_amplitude, final_phase);
        if(!start_sine)
        {   // Likely reason for not reaching start_sine is that the lattice is too small. Extend it and try again.
            lattice_extensions++;
            if(lattice_extensions > 5)
            {   *errstream << "ContinuumBuilder::CalculateContinuumWave:\n"
                           << "    start_sine not reached; lattice_extensions = " << lattice_extensions << std::endl;
                return 0;
            }

            lattice->R(HFPotential.size()+1);
            core->ExtendPotential();
            s->Clear();
            s->ReSize(HFPotential.size());
            *logstream << "Resizing continuum lattice; new size = " << HFPotential.size() << std::endl;
            old_phase = 0.0;
            start_sine = 0;
            ds = 1.0;   // Do another loop!
        }
        else
        {   ds = (final_phase - old_phase)/MathConstant::Instance()->Pi();
            if(fabs(ds) > Core::StateParameters::EnergyTolerance)
            {   core->CalculateExchange(*s, new_exchange);
                exchange.ReSize(new_exchange.Size());
                for(unsigned int i = 0; i<exchange.Size(); i++)
                    exchange.f[i] = 0.5 * exchange.f[i] + 0.5 * new_exchange.f[i];
                old_phase = final_phase;
            }
        }
    }
    while((loop < Core::StateParameters::MaxHFIterations) && (fabs(ds) > Core::StateParameters::EnergyTolerance));

    // Actual amplitude of wavefunction as r->Infinity (from IntegrateContinuum),
    //      A = final_amplitude/(2E)^(1/4)
    switch(norm_type)
    {
        case Flambaum:
            // Flambaum normalization:  A = 2 * Pi^(-1/2) * E^(1/2)
            final_amplitude = sqrt(2./(MathConstant::Instance()->Pi()*pow(s->GetNu(),3.)))/final_amplitude;
            break;

        case Cowan:
            // Cowan normalization:     A = Pi^(-1/2) * (2/E)^(1/4)
            final_amplitude = sqrt(2./MathConstant::Instance()->Pi())/final_amplitude;
            break;

        case Unitary:
            // Unitary normalization:   A = 1
            final_amplitude = sqrt(sqrt(2.*s->GetEnergy()))/final_amplitude;
            break;
    }

    (*s) *= final_amplitude;

    // Interpolate back onto external lattice
    if(external_lattice && !(*external_lattice == *lattice))
    {
        // Copy current state
        ContinuumWave cs_old(*s);

        Interpolator interp(lattice);
        const double* extR = external_lattice->R();
        unsigned int order = 6;
        double dfdr, dgdr;

        // Determine size of new continuum state
        s->Clear();
        unsigned int size = external_lattice->Size();
        if(external_lattice->MaxRealDistance() > lattice->MaxRealDistance())
        {   size = external_lattice->real_to_lattice(lattice->MaxRealDistance());
        }
        s->ReSize(size);

        // Interpolate
        for(unsigned int i = 0; i < size; i++)
        {
            interp.Interpolate(cs_old.f, extR[i], s->f[i], dfdr, order);
            interp.Interpolate(cs_old.g, extR[i], s->g[i], dgdr, order);
            s->dfdr[i] = dfdr;
            s->dgdr[i] = dgdr;
        }
    }

    if(DebugOptions.LogHFContinuum())
    {
        *logstream << std::setprecision(8);

        double energy = s->GetEnergy();
        if(DebugOptions.InvCmEnergyUnits())
            energy *= MathConstant::Instance()->HartreeEnergyInInvCm();
        *logstream << s->Name() << "  E = " << energy;

        *logstream << std::setprecision(4);
        *logstream << "  loops: " << loop << "  start sine: (" << start_sine << ") " << lattice->R(start_sine) << std::endl;
    }

    return loop;
}

bool ContinuumBuilder::ReadContinuumWave(ContinuumWave* s, Lattice* external_lattice, const std::string& upper_file, const std::string& lower_file)
{
    SetNormalisationType(Unitary);

    /*  Each file is formatted as follows: first line has energy, L.
        Subsequent lines have lattice point and wavefunction value (either upper or lower), one per line.
     */
    FILE* fp_upper = fopen(upper_file.c_str(), "rt");
    FILE* fp_lower = fopen(lower_file.c_str(), "rt");

    if(fp_upper && fp_lower)
    {
        // Get Energy and check L
        double energy_upper, energy_lower;
        int l_upper, l_lower;
        fscanf(fp_upper, "%lf%d", &energy_upper, &l_upper);
        fscanf(fp_lower, "%lf%d", &energy_lower, &l_lower);
        s->SetEnergy(energy_upper);
        if(energy_upper != energy_lower)
        {   *errstream << "ReadContinuumWave: energy doesn't match" << std::endl;
            return false;
        }
        if((l_upper != l_lower) || (l_lower != s->L()))
        {   *errstream << "ReadContinuumWave: l doesn't match" << std::endl;
            return false;
        }

        // Read upper component
        double r, val;
        std::vector<double> R, f, g;
        while(fscanf(fp_upper, "%lf%lf", &r, &val) == 2)
        {   R.push_back(r);
            f.push_back(val);
        }

        // Read lower component
        unsigned int i = 0;
        double alpha = PhysicalConstant::Instance()->GetAlpha();
        while(fscanf(fp_lower, "%lf%lf", &r, &val) == 2)
        {   if(R[i] != r)
            {   *errstream << "ReadContinuumWave: lattice points don't match at "
                           << R[i] << " (i = " << i << ")" << std::endl;
                return false;
            }
            g.push_back(val/alpha);
            i++;
        }

        // Interpolate onto external lattice
        unsigned int size;
        if(external_lattice->MaxRealDistance() < R[R.size()-1])
            size = external_lattice->Size();
        else
            size = external_lattice->real_to_lattice(R[R.size()-1]);
        s->ReSize(size);

        unsigned int order = 6;
        Interpolator interp(R, order);

        const double* extR = external_lattice->R();
        const double* extdR = external_lattice->dR(); 

        double dfdr, dgdr;
        for(unsigned int i = 0; i < s->Size(); i++)
        {
            interp.Interpolate(f, extR[i], s->f[i], dfdr, order);
            interp.Interpolate(g, extR[i], s->g[i], dgdr, order);
            s->dfdr[i] = dfdr;
            s->dgdr[i] = dgdr;
        }
    }
    
    return (fp_upper && fp_lower);
}
