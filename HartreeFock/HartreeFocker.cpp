#include "HartreeFocker.h"
#include "Include.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/MathConstant.h"
#include "ODESolver.h"
#include "GreensMethodODE.h"
#include "LocalPotentialDecorator.h"
#include "ThomasFermiDecorator.h"

void HartreeFocker::StartCore(Core* core, pHFOperator hf)
{
    // Sanity check
    if(core->GetOccupancies().size() != core->NumStates())
    {   *errstream << "HartreeFocker::StartCore(): core sanity check failed. Occupancies do not match states.";
        exit(1);
    }

    double Z = hf->GetZ();

    // Get first approximation to state energies
    double num_states = (double)core->NumStates();
    if(num_states > 1.)
    {
        StateIterator it = core->GetStateIterator();
        double i = 1.;
        while(!it.AtEnd())
        {
            double ratio = (i-1.) / (num_states-1.);
            double trial_nu = (num_states - i)/(Z * (num_states - 1.)) + ratio;
            trial_nu = trial_nu / (-4.*ratio*ratio + 4.*ratio + 1.);

            (it.GetState())->SetNu(trial_nu);

            it.Next();
            i += 1.;
        }

        // different for some states when 56 < Z < 72
        if(((56 < Z) && (Z < 72)) || (Z == 81))
        {
            // Alter 6s state
            pOrbital s = core->GetState(OrbitalInfo(6, -1));
            if(s != NULL)
                s->SetNu(1.4);

            // Alter 4f state
            s = core->GetState(OrbitalInfo(4, 3));
            if(s != NULL)
                s->SetNu(1.);

            // Alter 4f* state
            s = core->GetState(OrbitalInfo(4, -4));
            if(s != NULL)
                s->SetNu(1.);
        }
    }
    else if(!core->Empty())
    {
        pSingleParticleWavefunction s = core->GetStateIterator().GetState();
        s->SetNu(1./Z);
    }
    else
        return;

    // Wrap hf with local exchange approximation and Thomas-Fermi potential.
    pHFOperator hf_local_exch(new LocalExchangeApproximation(hf));
    pThomasFermiDecorator thomas_fermi(new ThomasFermiDecorator(hf_local_exch));

    thomas_fermi->SetCore(core);

    // Iterate states, slowly adding in direct HF + local exchange approximation.
    bool debug = DebugOptions.LogFirstBuild();

    for(unsigned int m = 0; m < 10; m++)
    {
        if(debug)
            *logstream << "First Build Iteration: " << m+1 << std::endl;

        StateIterator it = core->GetStateIterator();
        while(!it.AtEnd())
        {
            pOrbital s = it.GetState();
            unsigned int iterations = 0;
            double nu_change_factor = 0.25;
            int zero_difference = 0;        // Difference between required and actual number of nodes of wavefunction
            int previous_zero_difference = 0;
            bool is_first_iteration = true;

            double trial_nu = s->GetNu();
            double starting_nu = trial_nu;
            double starting_nu_change_factor = nu_change_factor;
            double starting_nu_percentage = 1.00;

            do
            {   s->Clear();

                // Attempt the local exchange approximation method.
                // This is encapsulated in the UpdateHFPotential method
                // when first_build = true.
                iterations = IterateOrbitalTailMatching(s, thomas_fermi);
                if(iterations >= MaxHFIterations)
                {   // Evil - this should never happen given our simple model.
                    *errstream << "    HartreeFocker::StartCore (" << m << ") " << s->Name()
                               << ": iterations = " << iterations << std::endl;
                    if(m == 9)
                    {   *outstream << "    HartreeFocker::StartCore (" << m << ") " << s->Name()
                                   << ": iterations = " << iterations << std::endl;
                        *outstream << "       ...try modifying energy/wavefunction tolerances." << std::endl;
                        exit(1);
                    }
                }

                zero_difference = s->NumNodes() + s->L() + 1 - s->GetPQN();

                if(zero_difference)
                {
                    if(debug)
                        *logstream << "    Zero difference: " << zero_difference << "  nu = " << trial_nu << std::endl;
                    s->SetNu(trial_nu - zero_difference/abs(zero_difference) * nu_change_factor * trial_nu);
                    trial_nu = s->GetNu();
                    if(!(((zero_difference > 0) && (previous_zero_difference > 0)) || ((zero_difference < 0) && (previous_zero_difference < 0))) && !is_first_iteration)
                    {
                        nu_change_factor = nu_change_factor * 0.75;
                    }
                }
                previous_zero_difference = zero_difference;
                is_first_iteration = false;
                if(s->GetNu() < 1.0e-5)
                {
                    starting_nu_percentage -= 0.1;
                    if(starting_nu_percentage <= 0.0)
                    {
                        *errstream << "HartreeFocker::StartCore: Cannot converge Hartree Fock energy." << std::endl;
                        exit(1);
                    }
                    trial_nu = starting_nu;
                    nu_change_factor = starting_nu_change_factor * starting_nu_percentage;
                }
            }
            while(zero_difference);

            if(debug)
            {   if(DebugOptions.HartreeEnergyUnits())
                *logstream << "  " << s->Name() << "  en: " << s->GetEnergy() << std::endl;
            else
                *logstream << "  " << s->Name() << " nu:   " << s->GetNu() << std::endl;
            }

            it.Next();
        }

        thomas_fermi->SetCore(core, 0.4);
    }
}

/** Iterate all orbitals in core until self-consistency is reached. */
void HartreeFocker::SolveCore(Core* core, pHFOperator hf)
{
    bool debug = DebugOptions.LogHFIterations();

    // Iterate Hartree-Fock. At each step:
    // 1. Calculate new solutions to the potential.
    //    Keep the new solutions (states) separate from the core so that the direct and
    //    exchange potentials are consistent.
    // 2. Mix old wavefunctions with new ones.
    // 3. Update potentials.

    Core next_states(core->Copy());
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
            pOrbital new_state = it.GetState();
            double old_energy = new_state->GetEnergy();
            
            pSpinorFunction exchange(new SpinorFunction(hf->GetExchange(new_state)));

            deltaE = IterateOrbital(new_state, hf, exchange);

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
            pOrbital core_state = core_it.GetState();
            pOrbital new_state = next_states.GetState(OrbitalInfo(core_state));
            
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
double HartreeFocker::SolveOrbital(pOrbital orbital, pHFOperator hf)
{
    double E = orbital->GetEnergy();
    double initial_energy = E;
    double delta_E = 0.0;

    if(!orbital->Size())
    {   pHFOperator hf_local(new LocalPotentialDecorator(hf));
        IterateOrbitalTailMatching(orbital, hf_local);
    }

    pSpinorFunction exchange(new SpinorFunction(hf->GetExchange(orbital)));

    double prop_new = 0.5;

    unsigned int loop;
    for(loop = 0; loop < MaxHFIterations; loop++)
    {
        // Get solution with current exchange
        IterateOrbital(orbital, hf, exchange);
        delta_E = orbital->GetEnergy() - E;

        // Check number of nodes for pqn
        //*outstream << "Loop " << loop << ": " << E << " " << delta_E << std::endl;

        if(fabs(delta_E/E) < EnergyTolerance)
            break;

        // Adjust exchange
        SpinorFunction new_exchange = hf->GetExchange(orbital);
        *exchange = *exchange * (1. - prop_new) + new_exchange * prop_new;
        E = E  + delta_E * prop_new;
        orbital->SetEnergy(E);
    }

    return E - initial_energy;
}

unsigned int HartreeFocker::CalculateExcitedState(pOrbital orbital, pHFOperator hf)
{
    // Number of iterations required. Zero shows that the state existed previously.
    unsigned int loop = 0;
    const Core* core = hf->GetCore();

    pOrbitalConst core_state = core->GetState(OrbitalInfo(orbital));
    if(core_state != NULL)
    {   // Try to find in core (probably open shells).
        *orbital = *core_state;
    }
    else
    {   double Charge = hf->GetCharge();

        if(!Charge)
        {   *errstream << "Cannot calculate excited Hartree-Fock states: Charge = 0." << std::endl;
            exit(1);
        }

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
                pOrbitalConst cs = i.GetState();
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
            pHFOperator hf_localexch(new LocalExchangeApproximation(hf));
            hf_localexch->SetCore(core);
            hf_localexch->IncludeExchangeInODE(false);
            loop = IterateOrbitalTailMatching(orbital, hf_localexch);
            //loop = IterateOrbitalTailMatching(orbital, hf);

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

            double deltaE;
            loop = 0;
            pOrbital new_orbital(new Orbital(*orbital));
            double prop_new = 0.4;
            pSpinorFunction exchange(new SpinorFunction(hf->GetExchange(orbital)));

            do
            {   loop++;

                int prev_zero_difference = 0;
                do
                {   if(zero_difference)
                    {
                        if(prev_zero_difference * zero_difference >= 0)
                            trial_nu = new_orbital->GetNu() - zero_difference/abs(zero_difference)/Charge;
                        else
                            trial_nu = new_orbital->GetNu() - 0.5 * zero_difference/abs(zero_difference)/Charge;

                        new_orbital->SetNu(trial_nu);
                        prev_zero_difference = zero_difference;
                    }
                    deltaE = IterateOrbital(new_orbital, hf, exchange);
                    zero_difference = new_orbital->NumNodes() + orbital->L() + 1 - orbital->GetPQN();
                } while(zero_difference);
                
                if(debugHF)
                    *logstream << "  " << std::setw(4) << orbital->Name()
                               << "  E = " << std::setprecision(12) << orbital->GetEnergy()
                               << "  deltaE = " << std::setprecision(3) << deltaE
                               << "  size: (" << new_orbital->Size()
                               << ") " << hf->GetLattice()->R(new_orbital->Size()) << std::endl;

                *exchange *= (1. - prop_new);
                *exchange += hf->GetExchange(new_orbital) * prop_new;

                *orbital *= (1. - prop_new);
                *orbital += (*new_orbital) * prop_new;

                // Renormalise core states (should be close already) and update energy.
                orbital->ReNormalise(hf->GetLattice());
                orbital->SetEnergy((1. - prop_new) * orbital->GetEnergy() + prop_new * new_orbital->GetEnergy());

                *new_orbital = *orbital;
                deltaE = fabs(deltaE/new_orbital->GetEnergy());
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
double HartreeFocker::IterateOrbital(pOrbital orbital, pHFOperator hf, pSpinorFunction exchange)
{
    pLattice lattice = hf->GetLattice();
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
        orbital->CheckSize(lattice, WavefunctionTolerance);
        orbital->ReNormalise(lattice);

        // Make sure hf operator is big enough.
        if(hf->Size() < orbital->Size())
        {   hf->ExtendPotential();

            // If it didn't work we might be in real trouble
            if(hf->Size() < orbital->Size())
            {   *errstream << "HartreeFocker::IterateOrbital(): hf->Size() < orbital->Size() and cannot be fixed.\n"
                << "    hf->Size() = " << hf->Size() << ", orbital->Size() = " << orbital->Size() << std::endl;
                orbital->ReSize(hf->Size());
            }
        }

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

unsigned int HartreeFocker::IterateOrbitalTailMatching(pOrbital orbital, pHFOperator hf)
{
    pLattice lattice = hf->GetLattice();
    const double Z = hf->GetZ();
    const double alpha = PhysicalConstant::Instance()->GetAlpha();
    AdamsSolver adamssolver(lattice);

    if(orbital->Size())
        orbital->CheckSize(lattice, WavefunctionTolerance);
    else
    {   double nu = mmax(orbital->GetNu(), 1./(9.*Z));
        orbital->SetNu(nu);
        double r_cutoff = (2. * hf->GetCharge() * nu + 10.) * nu;
        orbital->ReSize(lattice->real_to_lattice(r_cutoff));
    }

    double delta = 0.;
    unsigned int loop = 0;
    
    do
    {   loop++;

        // Make sure hf operator is big enough.
        if(hf->Size() < orbital->Size())
        {   hf->ExtendPotential();

            // If it didn't work we might be in real trouble
            if(hf->Size() < orbital->Size())
            {   *errstream << "HartreeFocker::IterateOrbitalTailMatching(): hf->Size() < orbital->Size() and cannot be fixed.\n"
                << "    hf->Size() = " << hf->Size() << ", orbital->Size() = " << orbital->Size() << std::endl;
                orbital->ReSize(hf->Size());
            }
        }

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

        odesolver->IntegrateForwards(hf, &*orbital);
        
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
    }
    while((loop < MaxHFIterations) && (fabs(delta) > TailMatchingEnergyTolerance));

    orbital->ReNormalise(lattice);
//    hf->GetDerivative(*orbital);

    return loop;
}
