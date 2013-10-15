#include "HartreeFocker.h"
#include "Include.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/MathConstant.h"
#include "ODESolver.h"
#include "GreensMethodODE.h"

/** Iterate all orbitals in core until self-consistency is reached. */
void HartreeFocker::SolveCore(Core* core, SpinorODE* hf)
{
    bool debug = DebugOptions.LogHFIterations();

    // Iterate Hartree-Fock. At each step:
    // 1. Calculate new solutions to the potential.
    //    Keep the new solutions (states) separate from the core so that the direct and
    //    exchange potentials are consistent.
    // 2. Mix old wavefunctions with new ones.
    // 3. Update potentials.

    Core next_states(*core);
    StateIterator it(&next_states);
    double prop_new = 0.5;
    
    double deltaE, max_deltaE;
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
            exchange = hf->GetExchange(new_state);

            deltaE = IterateOrbital(new_state, hf, &exchange);

            if(debug)
                *logstream << "  " << std::setw(4) << new_state->Name()
                << "  E = " << std::setprecision(12) << old_energy
                << "  deltaE = " << std::setprecision(4) << deltaE
                << "  size: (" << new_state->Size()
                << ") " << core->GetLattice()->R(new_state->Size()) << std::endl;
            
            deltaE = fabs(deltaE/new_state->GetEnergy());
            max_deltaE = mmax(deltaE, max_deltaE);
            
            it.Next();
        }
        
        // Mix new and old states.
        StateIterator core_it(core);
        core_it.First();
        while(!core_it.AtEnd())
        {
            Orbital* core_state = core_it.GetState();
            Orbital* new_state = next_states.GetState(OrbitalInfo(core_state));
            
            // Add proportion of new states to core states.
            *core_state *= (1. - prop_new);
            *new_state *= (prop_new);

            *core_state += *new_state;
            
            // Renormalise core states (should be close already) and update energy.
            core_state->ReNormalise(core->GetLattice());
            core_state->CheckSize(core->GetLattice(), WavefunctionTolerance);
            
            double energy = (1. - prop_new) * core_state->GetEnergy() + prop_new * new_state->GetEnergy();
            core_state->SetEnergy(energy);
            *new_state = *core_state;
            
            core_it.Next();
        }
        
        // Update potential.
        hf->SetCore(core);
        
    }while((max_deltaE > EnergyTolerance) && (loop < MaxHFIterations));

    if(loop >= MaxHFIterations)
        *errstream << "Failed to converge Hartree-Fock in Core." << std::endl;
    
    if(debug)
        *logstream << "Core Orthogonality test: " << core->TestOrthogonality() << std::endl;
}

/** Find self-consistent solution to hf operator, including exchange. */
double HartreeFocker::SolveOrbital(Orbital* orbital, SpinorODE* hf)
{
    double E = orbital->GetEnergy();
    double initial_energy = E;
    double delta_E = 0.0;

    if(!orbital->Size())
        hf->GetCore()->ConvergeStateApproximation(orbital);

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

        if(fabs(delta_E/E) < EnergyTolerance)
            break;

        // Adjust exchange
        SpinorFunction new_exchange = hf->GetExchange(orbital);
        exchange = exchange * (1. - prop_new) + new_exchange * prop_new;
        E = E  + delta_E * prop_new;
        orbital->SetEnergy(E);
    }

    return E - initial_energy;
}

unsigned int HartreeFocker::CalculateExcitedState(Orbital* orbital, SpinorODE* hf)
{
    // Number of iterations required. Zero shows that the state existed previously.
    unsigned int loop = 0;
    const Core* core = hf->GetCore();

    const Orbital* core_state = core->GetState(OrbitalInfo(orbital));
    if(core_state != NULL)
    {   // Try to find in core (probably open shells).
        *orbital = *core_state;
    }
    else
    {   if(!core->GetCharge())
        {   *errstream << "Cannot calculate excited Hartree-Fock states: Charge = 0." << std::endl;
            exit(1);
        }

        double Charge = core->GetCharge();

        // Get a trial value for nu if state hasn't already got one.
        double trial_nu;
        if(std::isnan(orbital->GetNu()))
        {
            switch(orbital->L())
            {   case 0:
                    trial_nu = 1.5;
                    break;
                case 1:
                    trial_nu = 2.;
                    break;
                case 2:
                    trial_nu = 4.;
                    break;
                default:
                    trial_nu = (double)orbital->L() + 1.;
                    break;
            }

            // Increase nu depending on how far above the core it is
            unsigned int largest_core_pqn = 0;
            ConstStateIterator i = core->GetConstStateIterator();
            while(!i.AtEnd())
            {
                const Orbital* cs = i.GetState();
                if((cs->Kappa() == orbital->Kappa()) && (cs->GetPQN() > largest_core_pqn))
                    largest_core_pqn = cs->GetPQN();
                i.Next();
            }
            
            trial_nu = trial_nu + (double)(orbital->GetPQN() - largest_core_pqn - 1)/Charge;
            orbital->SetNu(trial_nu/Charge);
        }
        
        // Calculate wavefunction with local exchange approximation.
        // Check that the number of nodes of the wavefunction is correct, otherwise
        // adjust nu and start again.
        
        int zero_difference = 0;    // Difference between required and actual number of nodes of wavefunction
        
        do
        {   // TODO: Implement local exchange approximation
            hf->IncludeExchangeInODE(false);
            loop = IterateOrbitalTailMatching(orbital, hf);

            if(loop >= MaxHFIterations)
            {   // Evil - this should never happen given our simple model.
                *errstream << "CalculateExcitedState: Local exchange approximation failed to converge.\n"
                           << "  " << orbital->Name() << "  iterations = " << loop << std::endl;
                PAUSE
                exit(1);
            }
            
            zero_difference = orbital->NumNodes() + orbital->L() + 1 - orbital->GetPQN();
            *logstream << orbital->Name() << " zero diff: " << zero_difference << ", E = " << orbital->GetEnergy() << std::endl;

            if(zero_difference)
            {   // This moves the state one pqn at a time
                trial_nu = orbital->GetNu() - zero_difference/abs(zero_difference)/Charge;
                orbital->SetNu(trial_nu);
                orbital->Clear();
            }
        }
        while(zero_difference || (loop >= MaxHFIterations));
        
        if(!core->Empty())
        {
            // Hartree-Fock loops
            bool debugHF = DebugOptions.LogHFIterations();

            // Adjust the energy slightly, otherwise the Wronskian for the Greens method may be zero.
            orbital->SetEnergy(1.1 * orbital->GetEnergy());

            double deltaE;
            loop = 0;
            Orbital new_orbital(*orbital);
            double prop_new = 0.4;
            do
            {   loop++;
                SpinorFunction exchange = hf->GetExchange(orbital);
                do
                {   if(zero_difference)
                    {   trial_nu = new_orbital.GetNu() - zero_difference/abs(zero_difference)/Charge;
                        new_orbital.SetNu(trial_nu);
                    }
                    deltaE = IterateOrbital(&new_orbital, hf, &exchange);
                    zero_difference = new_orbital.NumNodes() + orbital->L() + 1 - orbital->GetPQN();
                } while(zero_difference);
                
                if(debugHF)
                    *logstream << "  " << std::setw(4) << orbital->Name()
                               << "  E = " << std::setprecision(12) << orbital->GetEnergy()
                               << "  deltaE = " << std::setprecision(3) << deltaE
                               << "  size: (" << new_orbital.Size()
                               << ") " << hf->GetLattice()->R(new_orbital.Size()) << std::endl;

                *orbital *= (1. - prop_new);
                *orbital += new_orbital * prop_new;

                // Renormalise core states (should be close already) and update energy.
                orbital->ReNormalise(hf->GetLattice());
                orbital->SetEnergy((1. - prop_new) * orbital->GetEnergy() + prop_new * new_orbital.GetEnergy());
                new_orbital = *orbital;
                deltaE = fabs(deltaE/new_orbital.GetEnergy());
            }while((deltaE > EnergyTolerance) && (loop < MaxHFIterations));

            if(loop >= MaxHFIterations)
                *errstream << "Core: Failed to converge excited HF state " << orbital->Name() << std::endl;
        }
    }
    
    if(DebugOptions.OutputHFExcited())
    {
        *outstream << std::setprecision(12);
        if(DebugOptions.HartreeEnergyUnits() || DebugOptions.InvCmEnergyUnits())
        {
            double energy = orbital->GetEnergy();
            if(DebugOptions.InvCmEnergyUnits())
                energy *= MathConstant::Instance()->HartreeEnergyInInvCm();
            *outstream << orbital->Name() << "  E = " << energy;
        }
        else
            *outstream << orbital->Name() << "  nu = " << orbital->GetNu();
        
        *outstream << "  loops: " << loop << "  size: (" << orbital->Size() << ") " << hf->GetLattice()->R(orbital->Size()) << std::endl;
    }
    
    return loop;
}

/** Find energy eigenvalue for orbital with a given exchange potential.
 If exchange is not given, generate from hf.
 Note: this function does not iterate/update the exchange potential,
 so the final orbital is not an eigenvalue of the hf operator.
 */
double HartreeFocker::IterateOrbital(Orbital* orbital, SpinorODE* hf, SpinorFunction* exchange)
{
    Lattice* lattice = hf->GetLattice();
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
    
    } while (fabs(delta_E/E) > EnergyTolerance);

    // Get better derivative
    hf->SetODEParameters(orbital->Kappa(), orbital->GetEnergy(), &ex);
    hf->IncludeExchangeInODE(true);
    hf->GetDerivative(*orbital);

    return orbital->GetEnergy() - initial_energy;
}

unsigned int HartreeFocker::IterateOrbitalTailMatching(Orbital* orbital, SpinorODE* hf)
{
    Lattice* lattice = hf->GetLattice();
    const double Z = hf->GetCore()->GetZ();
    const double alpha = PhysicalConstant::Instance()->GetAlpha();
    AdamsSolver adamssolver(lattice);

    if(orbital->Size())
        orbital->CheckSize(lattice, WavefunctionTolerance);
    else
    {   double nu = mmax(orbital->GetNu(), 1./(9.*Z));
        orbital->SetNu(nu);
        double r_cutoff = (2. * hf->GetCore()->GetCharge() * nu + 10.) * nu;
        orbital->ReSize(lattice->real_to_lattice(r_cutoff));
    }

//    if(lattice->Size() > HFPotential.size())
//        ExtendPotential();
    
//    std::vector<double> Potential = GetHFPotential();
//    if(include_exch)
//    {   std::vector<double> LocalExchange = GetLocalExchangeApproximation();
//        for(unsigned int i=0; i<Potential.size(); i++)
//            Potential[i] = Potential[i] + LocalExchange[i];
//    }

    hf->SetODEParameters(orbital);
    
    double delta = 0.;
    unsigned int loop = 0;
    
    do
    {   loop++;

        // Get forward and backwards meeting point
        unsigned int critical_point = lattice->real_to_lattice((1./Z + orbital->GetNu())/2.);

        // Integrate backwards until first maximum
        hf->SetODEParameters(orbital->Kappa(), orbital->GetEnergy());
        Orbital backwards(*orbital);
        unsigned int peak = adamssolver.IntegrateBackwardsUntilPeak(hf, &backwards);
        if(peak < critical_point)
            peak = critical_point;

        double f_right = backwards.f[peak];
        double g_right = backwards.g[peak];

        odesolver->IntegrateForwards(hf, orbital);
        
        // Test whether forwards and backwards integrations met and modify energy accordingly.

        // Scale the backwards integration to meet the forwards one in f[peak].
        double tail_scaling = orbital->f[peak]/f_right;
        unsigned int i;
        for(i = peak+1; i < orbital->Size(); i++)
        {   orbital->f[i] = tail_scaling * backwards.f[i];
            orbital->dfdr[i] = tail_scaling * backwards.dfdr[i];
            orbital->g[i] = tail_scaling * backwards.g[i];
            orbital->dgdr[i] = tail_scaling * backwards.dgdr[i];
        }
        
        double norm = orbital->Norm(lattice);
        delta = pow(orbital->GetNu(), 3.) * orbital->f[peak] * (tail_scaling * g_right - orbital->g[peak])/norm/alpha;

        if(fabs(delta) > 0.1)
            delta = 0.1 * delta/fabs(delta);
        
        // Set new principal quantum number
        orbital->SetNu(orbital->GetNu() - delta);
        orbital->ReNormalise(lattice);
        
        orbital->CheckSize(lattice, WavefunctionTolerance);
//        if(lattice->Size() > HFPotential.size())
//        {   // Lengthen potential if necessary. Local exchange is zero outside the core.
//            ExtendPotential();
//            unsigned int old_size = Potential.size();
//            Potential.resize(HFPotential.size());
//            for(unsigned int i=old_size; i<HFPotential.size(); i++)
//                Potential[i] = HFPotential[i];
//        }
    }
    while((loop < MaxHFIterations) && (fabs(delta) > TailMatchingEnergyTolerance));

    orbital->ReNormalise(lattice);
    hf->GetDerivative(*orbital);

    return loop;
}

