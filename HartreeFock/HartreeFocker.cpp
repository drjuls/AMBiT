#include "HartreeFocker.h"
#include "Include.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/MathConstant.h"
#include "ODESolver.h"
#include "GreensMethodODE.h"
#include "LocalPotentialDecorator.h"
#include "ThomasFermiDecorator.h"
#include "Universal/Interpolator.h"

void HartreeFocker::StartCore(pCore core, pHFOperator hf)
{
    // Sanity check
    if(core->GetOccupancies().size() != core->size())
    {   *errstream << "HartreeFocker::StartCore(): core sanity check failed. Occupancies do not match states.";
        exit(1);
    }

    double Z = hf->GetZ();

    // Get first approximation to state energies
    double num_states = (double)core->size();
    if(num_states > 1.)
    {
        auto it = core->begin();
        double i = 1.;
        while(it != core->end())
        {
            double ratio = (i-1.) / (num_states-1.);
            double trial_nu = (num_states - i)/(Z * (num_states - 1.)) + ratio;
            trial_nu = trial_nu / (-4.*ratio*ratio + 4.*ratio + 1.);

            it->second->SetNu(trial_nu);

            it++;
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
    else if(!core->empty())
    {
        pOrbital s = core->begin()->second;
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

        auto it = core->begin();
        while(it != core->end())
        {
            pOrbital s = it->second;
            unsigned int iterations = 0;
            double nu_change_factor = 0.25;
            int zero_difference = 0;        // Difference between required and actual number of nodes of wavefunction
            int previous_zero_difference = 0;
            bool is_first_iteration = true;

            double trial_nu = s->Nu();
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

                zero_difference = s->NumNodes() + s->L() + 1 - s->PQN();

                if(zero_difference)
                {
                    if(debug)
                        *logstream << "    Zero difference: " << zero_difference << "  nu = " << trial_nu << std::endl;
                    s->SetNu(trial_nu - zero_difference/abs(zero_difference) * nu_change_factor * trial_nu);
                    trial_nu = s->Nu();
                    if(!(((zero_difference > 0) && (previous_zero_difference > 0)) || ((zero_difference < 0) && (previous_zero_difference < 0))) && !is_first_iteration)
                    {
                        nu_change_factor = nu_change_factor * 0.75;
                    }
                }
                previous_zero_difference = zero_difference;
                is_first_iteration = false;
                if(s->Nu() < 1.0e-5)
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
                *logstream << "  " << s->Name() << "  en: " << s->Energy() << std::endl;
            else
                *logstream << "  " << s->Name() << " nu:   " << s->Nu() << std::endl;
            }

            it++;
        }

        thomas_fermi->SetCore(core, 0.4);
    }
}

/** Iterate all orbitals in core until self-consistency is reached. */
void HartreeFocker::SolveCore(pCore core, pHFOperator hf)
{
    bool debug = DebugOptions.LogHFIterations();

    // Iterate Hartree-Fock. At each step:
    // 1. Calculate new solutions to the potential.
    //    Keep the new solutions (states) separate from the core so that the direct and
    //    exchange potentials are consistent.
    // 2. Mix old wavefunctions with new ones.
    // 3. Update potentials.

    pCore next_states(new Core(core->Clone()));
    double prop_new = 0.5;
    
    double deltaE, max_deltaE;
    unsigned int loop = 0;

    hf->SetCore(core);

    bool include_exchange = hf->IncludeExchange();

    do
    {   loop++;
        max_deltaE = 0.;
        
        if(debug)
            *logstream << "HF Iteration :" << loop << std::endl;
        
        // Calculate new states.
        auto it = next_states->begin();
        while(it != next_states->end())
        {
            pOrbital new_state = it->second;
            double old_energy = new_state->Energy();

            if(include_exchange)
            {
                pSpinorFunction exchange(new SpinorFunction(hf->GetExchange(new_state)));
                deltaE = IterateOrbital(new_state, hf, exchange);
            }
            else
            {   IterateOrbitalTailMatching(new_state, hf);
                deltaE = new_state->Energy() - old_energy;
            }

            if(debug)
                *logstream << "  " << std::setw(4) << new_state->Name()
                << "  E = " << std::setprecision(12) << old_energy
                << "  deltaE = " << std::setprecision(4) << deltaE
                << "  size: (" << new_state->size()
                << ") " << core->GetLattice()->R(new_state->size()) << std::endl;
            
            deltaE = fabs(deltaE/new_state->Energy());
            max_deltaE = mmax(deltaE, max_deltaE);
            
            it++;
        }
        
        // Mix new and old states.
        auto core_it = core->begin();
        while(core_it != core->end())
        {
            pOrbital core_state = core_it->second;
            pOrbital new_state = next_states->GetState(OrbitalInfo(core_state));
            
            // Add proportion of new states to core states.
            *core_state *= (1. - prop_new);
            *new_state *= (prop_new);

            *core_state += *new_state;
            
            // Renormalise core states (should be close already) and update energy.
            core_state->ReNormalise(odesolver->GetIntegrator());
            core_state->CheckSize(core->GetLattice(), WavefunctionTolerance);
            
            double energy = (1. - prop_new) * core_state->Energy() + prop_new * new_state->Energy();
            core_state->SetEnergy(energy);
            *new_state = *core_state;
            
            core_it++;
        }
        
        // Update potential.
        hf->SetCore(core);
        
    }while((max_deltaE > EnergyTolerance) && (loop < MaxHFIterations));

    if(loop >= MaxHFIterations)
        *errstream << "Failed to converge Hartree-Fock in Core." << std::endl;
}

/** Find self-consistent solution to hf operator, including exchange. */
double HartreeFocker::SolveOrbital(pOrbital orbital, pHFOperator hf)
{
    double E = orbital->Energy();
    double initial_energy = E;
    double delta_E = 0.0;

    if(!orbital->size())
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
        delta_E = orbital->Energy() - E;

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
    pCoreConst core = hf->GetCore();

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
        if(std::isnan(orbital->Nu()))
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
            auto i = core->begin();
            while(i != core->end())
            {
                pOrbitalConst cs = i->second;
                if((cs->Kappa() == orbital->Kappa()) && (cs->PQN() > largest_core_pqn))
                    largest_core_pqn = cs->PQN();
                i++;
            }
            
            trial_nu = trial_nu + (double)(orbital->PQN() - largest_core_pqn - 1)/Charge;
            orbital->SetNu(trial_nu/Charge);
        }
        
        // Calculate wavefunction with local exchange approximation.
        // Check that the number of nodes of the wavefunction is correct, otherwise
        // adjust nu and start again.
        
        int zero_difference = 0;    // Difference between required and actual number of nodes of wavefunction
        
        do
        {   pHFOperator hf_localexch(new LocalExchangeApproximation(hf));
            hf_localexch->SetCore(core);
            hf_localexch->IncludeExchange(false);
            loop = IterateOrbitalTailMatching(orbital, hf_localexch);
            //loop = IterateOrbitalTailMatching(orbital, hf);

            if(loop >= MaxHFIterations)
            {   // Evil - this should never happen given our simple model.
                *errstream << "CalculateExcitedState: Local exchange approximation failed to converge.\n"
                           << "  " << orbital->Name() << "  iterations = " << loop << std::endl;
                PAUSE
                exit(1);
            }
            
            zero_difference = orbital->NumNodes() + orbital->L() + 1 - orbital->PQN();
            if(DebugOptions.LogFirstBuild())
                *logstream << orbital->Name() << " zero diff: " << zero_difference << ", E = " << orbital->Energy() << std::endl;

            if(zero_difference)
            {   // This moves the state one pqn at a time
                trial_nu = orbital->Nu() - zero_difference/abs(zero_difference)/Charge;
                orbital->SetNu(trial_nu);
                orbital->Clear();
            }
        }
        while(zero_difference || (loop >= MaxHFIterations));
        
        if(!core->empty())
        {
            // Hartree-Fock loops
            bool debugHF = DebugOptions.LogHFIterations();

            double deltaE;
            loop = 0;
            pOrbital new_orbital(new Orbital(*orbital));
            double prop_new = 0.4;
            hf->IncludeExchange(true);
            pSpinorFunction exchange(new SpinorFunction(hf->GetExchange(orbital)));

            do
            {   loop++;
                if(loop > 100)
                    prop_new = 0.1;

                int prev_zero_difference = 0;
                do
                {   if(zero_difference)
                    {
                        if(prev_zero_difference * zero_difference >= 0)
                            trial_nu = new_orbital->Nu() - zero_difference/abs(zero_difference)/Charge;
                        else
                            trial_nu = new_orbital->Nu() - 0.5 * zero_difference/abs(zero_difference)/Charge;

                        new_orbital->SetNu(trial_nu);
                        prev_zero_difference = zero_difference;
                    }

                    double start_energy = new_orbital->Energy();
                    IterateOrbitalTailMatching(new_orbital, hf, exchange);
                    deltaE = new_orbital->Energy() - start_energy;

                    zero_difference = new_orbital->NumNodes() + orbital->L() + 1 - orbital->PQN();
                } while(zero_difference);
                
                if(debugHF)
                    *logstream << "  " << std::setw(4) << orbital->Name()
                               << "  E = " << std::setprecision(12) << orbital->Energy()
                               << "  deltaE = " << std::setprecision(3) << deltaE
                               << "  size: (" << new_orbital->size()
                               << ") " << hf->GetLattice()->R(new_orbital->size()) << std::endl;

                *exchange *= (1. - prop_new);
                *exchange += hf->GetExchange(new_orbital) * prop_new;

                *orbital *= (1. - prop_new);
                *orbital += (*new_orbital) * prop_new;

                // Renormalise core states (should be close already) and update energy.
                orbital->ReNormalise(odesolver->GetIntegrator());
                orbital->SetEnergy((1. - prop_new) * orbital->Energy() + prop_new * new_orbital->Energy());

                *new_orbital = *orbital;
                deltaE = fabs(deltaE/new_orbital->Energy());
            }while((deltaE > EnergyTolerance) && (loop < MaxHFIterations));

            IterateOrbital(orbital, hf);

            if(loop >= MaxHFIterations)
                *errstream << "Core: Failed to converge excited HF state " << orbital->Name() << std::endl;
        }
    }
    
    if(DebugOptions.OutputHFExcited())
    {
        *outstream << std::setprecision(12);
        if(DebugOptions.HartreeEnergyUnits() || DebugOptions.InvCmEnergyUnits())
        {
            double energy = orbital->Energy();
            if(DebugOptions.InvCmEnergyUnits())
                energy *= MathConstant::Instance()->HartreeEnergyInInvCm();
            *outstream << orbital->Name() << "  E = " << energy;
        }
        else
            *outstream << orbital->Name() << "  nu = " << orbital->Nu();
        
        *outstream << "  loops: " << loop << "  size: (" << orbital->size() << ") " << hf->GetLattice()->R(orbital->size()) << std::endl;
    }
    
    return loop;
}

/** Find energy eigenvalue for orbital with a given exchange potential.
    If exchange is not given, generate from hf.
    Note: this function does not iterate/update the exchange potential,
    so the final orbital is not an eigenvalue of the hf operator.
 */
double HartreeFocker::IterateOrbital(pOrbital orbital, pHFOperator hf, pSpinorFunctionConst exchange)
{
    pLattice lattice = hf->GetLattice();
    const double alpha = hf->GetPhysicalConstant()->GetAlpha();
    pOPIntegrator integrator = hf->GetOPIntegrator();

    SpinorFunction ex(orbital->Kappa());
    if(exchange)
        ex = *exchange;
    else
        ex = hf->GetExchange(orbital);

    double delta_E = 0.0;
    double E = orbital->Energy();
    double initial_energy = E;

    do
    {   // Use greens method to iterate orbital
        orbital->CheckSize(lattice, WavefunctionTolerance);
        orbital->ReNormalise(integrator);

        E = orbital->Energy();

        // Get solutions to homogenous equation (no exchange)
        hf->SetODEParameters(orbital->Kappa(), E);
        hf->IncludeExchange(false);
        Orbital originregular(*orbital);
        Orbital infinityregular(*orbital);
        odesolver->IntegrateBackwards(hf, &infinityregular);
        odesolver->IntegrateForwards(hf, &originregular);

        GreensMethodODE greens(lattice);
        greens.SetHomogenousSolutions(originregular, infinityregular);

        RadialFunction G0(orbital->size());
        RadialFunction GInf(orbital->size());

        greens.SetSourceTerm(ex * alpha, true);
        odesolver->IntegrateForwards(&greens, &G0);
        greens.SetSourceTerm(ex * alpha, false);
        odesolver->IntegrateBackwards(&greens, &GInf);

        *orbital = originregular * GInf - infinityregular * G0;

        // Now modify energy if required
        double norm = orbital->Norm(integrator);

        // Get delta_psi using Greens operator with psi as the source term
        greens.SetSourceTerm(*orbital, true);
        odesolver->IntegrateForwards(&greens, &G0);
        greens.SetSourceTerm(*orbital, false);
        odesolver->IntegrateBackwards(&greens, &GInf);

        Orbital delta_psi = originregular * GInf - infinityregular * G0;
        delta_psi *= hf->GetPhysicalConstant()->GetAlpha();

        double var = integrator->GetInnerProduct(*orbital, delta_psi);

        delta_E = (1. - norm)/(2. * var);
        if(fabs(delta_E/E) > 0.5)
            delta_E *= 0.5 * fabs(E/delta_E);

        orbital->SetEnergy(E + delta_E);
    
    } while (fabs(delta_E/E) > EnergyTolerance);

    // Get better derivative
    hf->SetODEParameters(orbital->Kappa(), orbital->Energy(), &ex);
    hf->IncludeExchange(true);
    hf->GetDerivative(*orbital);

    return orbital->Energy() - initial_energy;
}

unsigned int HartreeFocker::IterateOrbitalTailMatching(pOrbital orbital, pHFOperator hf, pSpinorFunctionConst exchange)
{
    pLattice lattice = hf->GetLattice();
    const double Z = hf->GetZ();
    const double alpha = hf->GetPhysicalConstant()->GetAlpha();
    pOPIntegrator integrator = odesolver->GetIntegrator();
    AdamsSolver adamssolver(integrator);

    if(orbital->size())
        orbital->CheckSize(lattice, WavefunctionTolerance);
    else
    {   double nu = mmax(orbital->Nu(), 1./(9.*Z));
        orbital->SetNu(nu);
        double r_cutoff = (2. * hf->GetCharge() * nu + 10.) * nu;
        orbital->resize(lattice->real_to_lattice(r_cutoff));
        if(orbital->size() > lattice->size())
            lattice->resize(orbital->size());
    }

    double delta = 0.;
    unsigned int loop = 0;
    
    do
    {   loop++;

        // Get forward and backwards meeting point
        unsigned int critical_point = lattice->real_to_lattice((1./Z + orbital->Nu())/2.);

        // Integrate backwards until first maximum
        hf->SetODEParameters(orbital->Kappa(), orbital->Energy(), exchange.get());
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
        for(i = peak+1; i < orbital->size(); i++)
        {   orbital->f[i] = tail_scaling * backwards.f[i];
            orbital->dfdr[i] = tail_scaling * backwards.dfdr[i];
            orbital->g[i] = tail_scaling * backwards.g[i];
            orbital->dgdr[i] = tail_scaling * backwards.dgdr[i];
        }
        
        double norm = orbital->Norm(integrator);
        delta = pow(orbital->Nu(), 3.) * orbital->f[peak] * (tail_scaling * g_right - orbital->g[peak])/norm/alpha;

        if(fabs(delta) > 0.1)
            delta = 0.1 * delta/fabs(delta);
        
        // Set new principal quantum number
        orbital->SetNu(orbital->Nu() - delta);
        orbital->ReNormalise(integrator);
        
        orbital->CheckSize(lattice, WavefunctionTolerance);
    }
    while((loop < MaxHFIterations) && (fabs(delta) > TailMatchingEnergyTolerance));

    orbital->ReNormalise(integrator);

    return loop;
}

unsigned int HartreeFocker::CalculateContinuumWave(pContinuumWave s, pHFOperator hf)
{
    unsigned int loop = 0;
    double final_amplitude, final_phase;
    pSpinorFunction exchange(new SpinorFunction(s->Kappa()));
    double ds, old_phase = 0.;
    unsigned int start_sine = 0;
    bool all_done = false;

    pLattice lattice = hf->GetLattice();
    s->resize(lattice->size());
    double pi = MathConstant::Instance()->Pi();

    do
    {   loop++;

        start_sine = IntegrateContinuum(s, hf, exchange, final_amplitude, final_phase);
        if(!start_sine)
        {   // Likely reason for not reaching start_sine is that the lattice is too small. Extend it and try again.

            // Probably get worried if we've had, say, fifty complete oscillations and still no good.
            // k ~ sqrt(2(energy + V(r))) => 50 oscillations = 50 * 2 Pi/k
            double Vmin = hf->GetDirectPotential().f[lattice->size()-1];
            double kmin = sqrt(2. * (s->Energy() + Vmin));
            double Rmax = 100. * pi/kmin;

            if(lattice->MaxRealDistance() > mmax(Rmax, 1000))
            {   *errstream << "ContinuumBuilder::CalculateContinuumWave:\n"
                           << "    start_sine not reached; Rmax = " << hf->GetLattice()->MaxRealDistance() << std::endl;
                return 0;
            }

            lattice->resize_to_r(lattice->MaxRealDistance() * 2.);
            s->Clear();
            s->resize(lattice->size());
            old_phase = 0.0;
            start_sine = 0;
            ds = 1.0;   // Do another loop!
        }
        else
        {   ds = (final_phase - old_phase)/MathConstant::Instance()->Pi();
            if(fabs(ds) > EnergyTolerance)
            {
                *exchange *= 0.5;
                *exchange += hf->GetExchange(s) * 0.5;
                old_phase = final_phase;
            }
            else
            {   *exchange = hf->GetExchange(s);
                old_phase = final_phase;
                if(!all_done)
                {   all_done = true;
                    ds = 1.0;   // One more loop!
                }
            }
        }
    }
    while((loop < MaxHFIterations) && (fabs(ds) > EnergyTolerance));

    // Actual amplitude of wavefunction as r->Infinity (from IntegrateContinuum),
    //      A = final_amplitude/(2E)^(1/4)
    switch(continuum_normalisation_type)
    {
        case ContinuumNormalisation::LandauNu:
            // Normalization with respect to delta function in nu rather than energy:
            //      A = 2 * Pi^(-1/2) * E^(1/2)
            {   double nu = fabs(s->Nu());
                final_amplitude = sqrt(2./(pi*nu*nu*nu))/final_amplitude;
            }
            break;

        case ContinuumNormalisation::LandauEnergy:
        case ContinuumNormalisation::Cowan:
            // Cowan normalization:     A = Pi^(-1/2) * (2/E)^(1/4)
            final_amplitude = sqrt(2./pi)/final_amplitude;
            break;

        case ContinuumNormalisation::Unitary:
            // Unitary normalization:   A = 1
            final_amplitude = sqrt(sqrt(2.*s->Energy()))/final_amplitude;
            break;
    }

    (*s) *= final_amplitude;

    if(DebugOptions.LogHFContinuum())
    {
        *logstream << std::setprecision(8);

        double energy = s->Energy();
        if(DebugOptions.InvCmEnergyUnits())
            energy *= MathConstant::Instance()->HartreeEnergyInInvCm();
        *logstream << s->Name() << "  E = " << energy;

        *logstream << std::setprecision(4);
        *logstream << "  loops: " << loop << "  start sine: (" << start_sine << ") " << lattice->R(start_sine) << std::endl;
    }
    
    return loop;
}

unsigned int HartreeFocker::IntegrateContinuum(pContinuumWave s, pHFOperator hf, pSpinorFunction exchange, double& final_amplitude, double& final_phase)
{
    static double accuracy = 0.01;  // Accuracy of equality of amplitudes over one cycle

    // Start by getting wavefunction
    int kappa = s->Kappa();
    double energy = s->Energy();
    hf->SetODEParameters(kappa, energy, exchange.get());
    odesolver->IntegrateForwards(hf, s.get());

    // Now we need to check the boundary condition (sine wave) at infinity and return final amplitude and phase
    double pi = MathConstant::Instance()->Pi();
    const double* R = hf->GetLattice()->R();
    const RadialFunction hf_direct(hf->GetDirectPotential());

    // Get momentum function sqrt(E-V)
    RadialFunction P(hf_direct.size());
    double kappa_squared = kappa * (kappa + 1);
    for(unsigned int i = 0; i < P.size(); i++)
    {
        P.f[i] = sqrt(2. * fabs(energy + hf_direct.f[i] - kappa_squared/(2. * R[i] * R[i])));
        P.dfdr[i] = (hf_direct.dfdr[i] + kappa_squared/(R[i] * R[i] * R[i]))/P.f[i];
    }

    unsigned int size = P.size();
    RadialFunction S(size);

    unsigned int start_sine = 0;
    final_amplitude = 0.;
    double peak_phase = 0.;

    // Find sine
    ContinuumPhaseODE phase_ode(hf->GetLattice(), hf->GetDirectPotential(), kappa, energy);
    odesolver->IntegrateForwards(&phase_ode, &S);

    Interpolator I(hf->GetLattice());
    int i = hf->GetLattice()->real_to_lattice(1.e-4);    // Outside nucleus

    while(i < size && !start_sine)
    {
        if((hf_direct.f[i] < 2.5 * energy) && (s->dfdr[i-1]/s->dfdr[i-2] < 0.))
        {
            double x_max;
            double f_max = I.FindExtremum(s->f, i-2, x_max);
            bool maximum = (s->dfdr[i-2] > 0.);
            double p_max, dpdr;
            I.Interpolate(P.f, x_max, p_max, dpdr, 4);

            double old_amplitude = final_amplitude;
            final_amplitude = fabs(f_max * sqrt(p_max));
            if(fabs((final_amplitude - old_amplitude)/final_amplitude) < accuracy)
            {
                double old_peak_phase = peak_phase;
                I.Interpolate(S.f, x_max, peak_phase, dpdr, 4);
                if((old_peak_phase != 0.) && (fabs((peak_phase - old_peak_phase)/pi - 1.) < accuracy))
                {
                    start_sine = i;
                    if(!maximum)
                        peak_phase = peak_phase - pi;
                }
            }
            else
                peak_phase = 0.;    // Reset
        }

        i++;
    }

    if(start_sine)
    {
        final_phase = S.f[s->size() - 1] - peak_phase;
    }
    
    return start_sine;

}
