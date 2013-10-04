#include "HartreeFocker.h"
#include "Include.h"
#include "Universal/PhysicalConstant.h"
#include "ODESolver.h"
#include "GreensMethodODE.h"

/** Iterate all orbitals in core until self-consistency is reached. */
void HartreeFocker::SolveCore(Core* core, HFOperator* hf, SpinorODE* hf_decorated)
{
    /*
    bool debug = DebugOptions.LogHFIterations();
    
    // Copy states for next iteration.
    StateManager next_states(lattice, (unsigned int)Z, (int)Charge);
    
    StateIterator core_it(this);
    core_it.First();
    while(!core_it.AtEnd())
    {   next_states.AddState(new Orbital(*core_it.GetState()));
        core_it.Next();
    }
    
    // Iterate Hartree-Fock. At each step:
    // 1. Calculate new solutions to the potential.
    //    Keep the new solutions (states) separate from the core so that the direct and
    //    exchange potentials are consistent.
    // 2. Mix old wavefunctions with new ones.
    // 3. Update potentials.
    double deltaE, max_deltaE;
    StateIterator it(&next_states);
    double prop_new = 0.5;
    StateIntegrator SI(lattice);
    
    unsigned int loop = 0;
    
    do
    {   loop++;
        max_deltaE = 0.;
        
        if(debug)
            *logstream << "HF Iteration :" << loop << std::endl;
        
        // Calculate new states.
        it.First();
        while(!it.AtEnd())
        {
            Orbital* new_state = it.GetState();
            double old_energy = new_state->GetEnergy();
            
            SpinorFunction exchange(new_state->Kappa());
            CalculateExchange(*new_state, exchange);
            deltaE = IterateOrbitalGreens(new_state, &exchange);
            
            if(debug)
                *logstream << "  " << std::setw(4) << new_state->Name()
                << "  E = " << std::setprecision(12) << old_energy
                << "  deltaE = " << std::setprecision(4) << deltaE
                << "  size: (" << new_state->Size()
                << ") " << lattice->R(new_state->Size()) << std::endl;
            
            deltaE = fabs(deltaE/new_state->GetEnergy());
            max_deltaE = mmax(deltaE, max_deltaE);
            
            it.Next();
        }
        
        // Mix new and old states.
        core_it.First();
        while(!core_it.AtEnd())
        {
            Orbital* core_state = core_it.GetState();
            Orbital* new_state = next_states.GetState(OrbitalInfo(core_state));
            
            // Add proportion of new states to core states.
            *core_state *= (1. - prop_new);
            *new_state *= (prop_new);
            core_state->ReSize(mmax(core_state->Size(), new_state->Size()));
            new_state->ReSize(mmax(core_state->Size(), new_state->Size()));
            
            for(unsigned int i = 0; i<core_state->Size(); i++)
            {   core_state->f[i] += new_state->f[i];
                core_state->g[i] += new_state->g[i];
                core_state->dfdr[i] += new_state->dfdr[i];
                core_state->dgdr[i] += new_state->dgdr[i];
            }
            
            // Renormalise core states (should be close already) and update energy.
            core_state->ReNormalise(lattice);
            core_state->CheckSize(lattice, StateParameters::WavefunctionTolerance);
            
            double energy = (1. - prop_new) * core_state->GetEnergy() + prop_new * new_state->GetEnergy();
            core_state->SetEnergy(energy);
            *new_state = *core_state;
            
            core_it.Next();
        }
        
        // Update potential.
        UpdateHFPotential();
        
    }while((max_deltaE > StateParameters::EnergyTolerance) && (loop < StateParameters::MaxHFIterations));
    
    if(loop >= StateParameters::MaxHFIterations)
        *errstream << "Failed to converge Hartree-Fock in Core." << std::endl;
    
    if(debug)
        *logstream << "Core Orthogonality test: " << TestOrthogonality() << std::endl;
     */
}

/** Find self-consistent solution to hf operator, including exchange. */
double HartreeFocker::SolveOrbital(Orbital* orbital, SpinorODE* hf)
{
    double E = orbital->GetEnergy();
    double initial_energy = E;
    double delta_E = 0.0;
    
    SpinorFunction exchange = hf->GetExchange(orbital);
    double prop_new = 0.5;

    unsigned int loop;
    for(loop = 0; loop < MaxHFIterations; loop++)
    {
        // Get solution with current exchange
        IterateOrbital(orbital, hf, &exchange);
        delta_E = orbital->GetEnergy() - E;

        // Check number of nodes for pqn
        //*outstream << "Loop " << loop << ": " << E << " " << delta_E << std::endl;

        if(fabs(delta_E/E) < WavefunctionEnergyTolerance)
            break;

        // Adjust exchange
        SpinorFunction new_exchange = hf->GetExchange(orbital);
        exchange = exchange * (1. - prop_new) + new_exchange * prop_new;
        E = E  + delta_E * prop_new;
        orbital->SetEnergy(E);
    }

    return E - initial_energy;
}

/** Find energy eigenvalue for orbital with a given exchange potential.
 If exchange is not given, generate from hf.
 Note: this function does not iterate/update the exchange potential,
 so the final orbital is not an eigenvalue of the hf operator.
 */
double HartreeFocker::IterateOrbital(Orbital* orbital, SpinorODE* hf, SpinorFunction* exchange)
{
    Lattice* lattice = hf->GetLattice();
    ODESolver* odesolver = new AdamsSolver(lattice);
    const double alpha = PhysicalConstant::Instance()->GetAlpha();

    SpinorFunction ex(orbital->Kappa());
    if(exchange)
        ex = *exchange;
    else
        ex = hf->GetExchange(orbital);

    double delta_E = 0.0;
    double E = orbital->GetEnergy();
    double initial_energy = E;

    do
    {   // Use greens method to iterate orbital
        orbital->ReNormalise(lattice);
        orbital->CheckSize(lattice, WavefunctionTolerance);
        E = orbital->GetEnergy();

        // Get solutions to homogenous equation (no exchange)
        hf->SetODEParameters(orbital->Kappa(), E);
        hf->IncludeExchangeInODE(false);
        Orbital originregular(*orbital);
        Orbital infinityregular(*orbital);
        odesolver->IntegrateBackwards(hf, &infinityregular);
        odesolver->IntegrateForwards(hf, &originregular);

        GreensMethodODE greens(lattice);
        greens.SetHomogenousSolutions(originregular, infinityregular);

        RadialFunction G0(orbital->Size());
        RadialFunction GInf(orbital->Size());

        greens.SetSourceTerm(ex * alpha, true);
        odesolver->IntegrateForwards(&greens, &G0);
        greens.SetSourceTerm(ex * alpha, false);
        odesolver->IntegrateBackwards(&greens, &GInf);

        *orbital = originregular * GInf - infinityregular * G0;

        // Now modify energy if required
        double norm = orbital->Norm(lattice);

        // Get delta_psi using Greens operator with psi as the source term
        greens.SetSourceTerm(*orbital, true);
        odesolver->IntegrateForwards(&greens, &G0);
        greens.SetSourceTerm(*orbital, false);
        odesolver->IntegrateBackwards(&greens, &GInf);

        Orbital delta_psi = originregular * GInf - infinityregular * G0;
        delta_psi *= PhysicalConstant::Instance()->GetAlpha();

        double var = orbital->Overlap(delta_psi, lattice);

        delta_E = (1. - norm)/(2. * var);
        if(fabs(delta_E/E) > 0.5)
            delta_E *= 0.5 * fabs(E/delta_E);

        orbital->SetEnergy(E + delta_E);
    
    } while (fabs(delta_E/E) > WavefunctionEnergyTolerance);

    // Get better derivative
    hf->SetODEParameters(orbital->Kappa(), orbital->GetEnergy(), &ex);
    hf->IncludeExchangeInODE(true);
    hf->GetDerivative(*orbital);

    delete odesolver;
    return orbital->GetEnergy() - initial_energy;
}
