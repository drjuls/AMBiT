#include "HartreeFocker.h"
#include "Include.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/MathConstant.h"
#include "ODESolver.h"
#include "GreensMethodODE.h"
#include "LocalPotentialDecorator.h"
#include "ThomasFermiDecorator.h"
#include "Universal/Interpolator.h"

#define PRINT_HF_LOOP_ORBITALS false

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

    bool include_exchange = hf->IncludeExchange();
    thomas_fermi->IncludeExchange(false);

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

            s->Clear();
            ConvergeOrbital(s, thomas_fermi, pSpinorFunction(), &HartreeFocker::IterateOrbitalTailMatching, 1.e-6);

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

    hf->IncludeExchange(include_exchange);
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
    int zero_difference = 0;
    int abs_zero_diff = 0;
    double energy_tolerance = 1.e-6;

    hf->SetCore(core);
    bool include_exchange = hf->IncludeExchange();

    do
    {   loop++;
        max_deltaE = 0.;
        abs_zero_diff = 0;

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
                deltaE = ConvergeOrbital(new_state, hf, exchange, &HartreeFocker::IterateOrbital, energy_tolerance);
            }
            else
            {
                pSpinorFunction exchange(new SpinorFunction(new_state->Kappa()));
                deltaE = ConvergeOrbital(new_state, hf, exchange, &HartreeFocker::IterateOrbitalTailMatching, energy_tolerance);
            }

            zero_difference = new_state->NumNodes() + new_state->L() + 1 - new_state->PQN();
            abs_zero_diff += abs(zero_difference);

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

        if((energy_tolerance > EnergyTolerance) && (abs_zero_diff == 0))
            energy_tolerance = mmax(energy_tolerance * 0.1, EnergyTolerance);
        
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

unsigned int HartreeFocker::CalculateExcitedState(pOrbital orbital, pHFOperator hf)
{
    // Number of iterations required. Zero shows that the state existed previously.
    unsigned int loop = 0;
    pCoreConst core = hf->GetCore();
    bool exchange_included = hf->IncludeExchange();

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
        if(std::isnan(orbital->Nu()) || orbital->Energy() >= 0.)
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

        // Calculate wavefunction without exchange.
        if(!orbital->size())
        {
            hf->SetCore(core);
            hf->IncludeExchange(false);
            ConvergeOrbital(orbital, hf, pSpinorFunction(), &HartreeFocker::IterateOrbitalTailMatching, TailMatchingEnergyTolerance);
        }
    }

    // Hartree-Fock loops
    hf->IncludeExchange(exchange_included);
    ConvergeOrbitalAndExchange(orbital, hf, &HartreeFocker::IterateOrbitalTailMatching, TailMatchingEnergyTolerance);

    if(!core->empty() && exchange_included)
        ConvergeOrbitalAndExchange(orbital, hf, &HartreeFocker::IterateOrbital, EnergyTolerance);
    else
        ConvergeOrbitalAndExchange(orbital, hf, &HartreeFocker::IterateOrbitalTailMatching, EnergyTolerance);

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
        
        *outstream << "  size: (" << orbital->size() << ") " << hf->GetLattice()->R(orbital->size()) << std::endl;
    }
    
    return loop;
}

double HartreeFocker::ConvergeOrbital(pOrbital orbital, pHFOperator hf, pSpinorFunctionConst exchange, IteratorFunction iterator, double energy_tolerance)
{
    pLattice lattice = hf->GetLattice();
    pOPIntegrator integrator = hf->GetOPIntegrator();

    if(orbital->size())
    {   orbital->ReNormalise(integrator);
        orbital->CheckSize(lattice, WavefunctionTolerance);
    }
    else
    {   double nu = mmax(orbital->Nu(), 1./(9. * hf->GetZ()));
        orbital->SetNu(nu);
        double r_cutoff = (2. * mmax(1., hf->GetCharge()) * nu + 10.) * nu;
        orbital->resize(lattice->real_to_lattice(r_cutoff));
        if(orbital->size() > lattice->size())
            lattice->resize(orbital->size());
    }

    bool include_exchange = hf->IncludeExchange() && exchange;

    double delta_E = 0.0;
    double Eupper = 0.0;
    double Elower = -(hf->GetZ() * hf->GetZ());

    double E = orbital->Energy();
    double initial_energy = E;
    int zero_difference = 0;
    unsigned int loop = 0;

    do
    {   loop++;

        E = orbital->Energy();
        delta_E = iterator(*this, orbital, hf, exchange);

        // Check zeros and modify energy as required
        zero_difference = orbital->NumNodes() + orbital->L() + 1 - orbital->PQN();
        if(zero_difference > 0)
        {
            // Too many nodes: decrease E
            Eupper = E;
            delta_E = 0.1 * E;
        }
        else if(zero_difference < 0)
        {
            // Too few nodes: increase E
            Elower = E;
            delta_E = -0.1 * E;
        }
        else
        {   // Limit magnitude of delta_nu to 0.1
            double max_delta_E = (-2. * E);
            max_delta_E *= sqrt(max_delta_E) * 0.1;

            if(fabs(delta_E) > max_delta_E)
                delta_E = (delta_E > 0.0? max_delta_E: -max_delta_E);
        }

        if(E + delta_E >= Eupper)
            delta_E = (Eupper - E)/2.;
        else if(E + delta_E <= Elower)
            delta_E = (Elower - E)/2.;

        orbital->SetEnergy(E + delta_E);
        orbital->ReNormalise(integrator);
        orbital->CheckSize(lattice, WavefunctionTolerance);

        if(DebugOptions.LogHFInnerLoop())
            *logstream << "  " << std::setw(4) << orbital->Name()
                       << "  E = " << std::setprecision(12) << E
                       << "  deltaE = " << std::setprecision(3) << delta_E
                       << "  size: (" << orbital->size()
                       << ") " << hf->GetLattice()->R(orbital->size())
                       << "  zerodiff = " << zero_difference << std::endl;

    } while((loop < MaxHFIterations) && (fabs(delta_E/E) > energy_tolerance));

    if(loop >= MaxHFIterations)
        *errstream << "ConvergeOrbital: Failed to converge HF state " << orbital->Name() << std::endl;

    orbital->ReNormalise(integrator);
    orbital->CheckSize(lattice, WavefunctionTolerance);

    // Get better derivative
    if(include_exchange)
    {
        hf->SetODEParameters(orbital->Kappa(), orbital->Energy(), exchange.get());
        hf->IncludeExchange(true);
        hf->GetDerivative(*orbital);
    }

    return orbital->Energy() - initial_energy;
}

double HartreeFocker::ConvergeOrbitalAndExchange(pOrbital orbital, pHFOperator hf, IteratorFunction iterator, double energy_tolerance)
{
    pLattice lattice = hf->GetLattice();
    const double Z = hf->GetZ();
    pOPIntegrator integrator = odesolver->GetIntegrator();
    AdamsSolver adamssolver(integrator);

    pSpinorFunction exchange(new SpinorFunction(orbital->Kappa()));
    bool include_exchange = hf->IncludeExchange();

    if(orbital->size())
    {   orbital->ReNormalise(integrator);
        orbital->CheckSize(lattice, WavefunctionTolerance);
        if(include_exchange)
            *exchange = hf->GetExchange(orbital);
    }
    else
    {   double nu = mmax(orbital->Nu(), 1./(9.*Z));
        orbital->SetNu(nu);
        double r_cutoff = (2. * mmax(1., hf->GetCharge()) * nu + 10.) * nu;
        orbital->resize(lattice->real_to_lattice(r_cutoff));
        if(orbital->size() > lattice->size())
            lattice->resize(orbital->size());
    }

    unsigned int loop = 0;
    double E = orbital->Energy();
    double initial_energy = E;

    double delta_E = 0.0;
    double Eupper = 0.0;
    double Elower = -(hf->GetZ() * hf->GetZ());

    int zero_difference = 0;
    Orbital prev(*orbital);

    do
    {   loop++;
        E = orbital->Energy();

        delta_E = iterator(*this, orbital, hf, exchange);

        // Check zeros and modify energy as required
        zero_difference = orbital->NumNodes() + orbital->L() + 1 - orbital->PQN();
        if(zero_difference > 0)
        {
            // Too many nodes: decrease E
            Eupper = E;
            delta_E = 0.1 * E;

            #if PRINT_HF_LOOP_ORBITALS
            if(DebugOptions.LogHFInnerLoop())
            {
                std::string filename = orbital->Name() + "_" + itoa(loop) + ".orb";
                orbital->Print(filename, lattice);
            }
            #endif

            *orbital = prev;
            orbital->SetEnergy(E + delta_E);
        }
        else if(zero_difference < 0)
        {
            // Too few nodes: increase E
            Elower = E;
            delta_E = -0.1 * E;

            #if PRINT_HF_LOOP_ORBITALS
            if(DebugOptions.LogHFInnerLoop())
            {
                std::string filename = orbital->Name() + "_" + itoa(loop) + ".orb";
                orbital->Print(filename, lattice);
            }
            #endif

            *orbital = prev;
            orbital->SetEnergy(E + delta_E);
        }
        else
        {   // Limit magnitude of delta_nu to 0.1
            double original_delta_E = delta_E;
            double max_delta_E = (-2. * E);
            max_delta_E *= sqrt(max_delta_E) * 0.1;

            if(fabs(delta_E) > max_delta_E)
                delta_E = (delta_E > 0.0? max_delta_E: -max_delta_E);

            // Check upper and lower energy limits
            if(E + delta_E >= Eupper)
            {
                delta_E = (Eupper - E)/2.;
            }
            else if(E + delta_E <= Elower)
            {
                delta_E = (Elower - E)/2.;
            }

            orbital->SetEnergy(E + delta_E);

            // Mix old and new orbitals
            double prop_new = 0.1;
            if(fabs(original_delta_E) > energy_tolerance)
                prop_new *= fabs(delta_E/original_delta_E);

            (*orbital) *= prop_new;
            (*orbital) += prev * (1. - prop_new);

            orbital->CheckSize(lattice, WavefunctionTolerance);
            orbital->ReNormalise(integrator);
            prev = *orbital;

            #if PRINT_HF_LOOP_ORBITALS
            if(DebugOptions.LogHFInnerLoop())
            {
                std::string filename = orbital->Name() + "_" + itoa(loop) + ".orb";
                orbital->Print(filename, lattice);
            }
            #endif

            if(include_exchange)
                *exchange = hf->GetExchange(orbital);
        }

        if(DebugOptions.LogHFInnerLoop())
            *logstream << "  " << std::setw(4) << orbital->Name()
                       << "  E = " << std::setprecision(12) << E
                       << "  deltaE = " << std::setprecision(3) << delta_E
                       << "  size: (" << orbital->size()
                       << ") " << hf->GetLattice()->R(orbital->size())
                       << "  zerodiff = " << zero_difference << std::endl;

    } while((loop < MaxHFIterations) && (fabs(delta_E/E) > energy_tolerance));

    if(loop >= MaxHFIterations)
        *errstream << "ConvergeOrbitalAndExchange: Failed to converge HF state " << orbital->Name() << std::endl;

    // Get better derivative
    hf->SetODEParameters(orbital->Kappa(), orbital->Energy(), exchange.get());
    hf->IncludeExchange(include_exchange);
    hf->GetDerivative(*orbital);

    return orbital->Energy() - initial_energy;
}

double HartreeFocker::IterateOrbital(pOrbital orbital, pHFOperator hf, pSpinorFunctionConst exchange)
{
    pLattice lattice = hf->GetLattice();
    const double alpha = hf->GetPhysicalConstant()->GetAlpha();
    pOPIntegrator integrator = hf->GetOPIntegrator();
    const SpinorFunction& ex(*exchange);

    orbital->ReNormalise(integrator);
    orbital->CheckSize(lattice, WavefunctionTolerance);

    // Get solutions to homogenous equation (no exchange)
    hf->SetODEParameters(orbital->Kappa(), orbital->Energy());
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

    double norm = orbital->Norm(integrator);

    // Get delta_psi using Greens operator with psi as the source term
    greens.SetSourceTerm(*orbital, true);
    odesolver->IntegrateForwards(&greens, &G0);
    greens.SetSourceTerm(*orbital, false);
    odesolver->IntegrateBackwards(&greens, &GInf);

    Orbital delta_psi = originregular * GInf - infinityregular * G0;
    delta_psi *= hf->GetPhysicalConstant()->GetAlpha();

    double var = integrator->GetInnerProduct(*orbital, delta_psi);

    return (1. - norm)/(2. * var);
}

double HartreeFocker::IterateOrbitalTailMatching(pOrbital orbital, pHFOperator hf, pSpinorFunctionConst exchange)
{
    pLattice lattice = hf->GetLattice();
    const double Z = hf->GetZ();
    const double alpha = hf->GetPhysicalConstant()->GetAlpha();
    pOPIntegrator integrator = odesolver->GetIntegrator();
    AdamsSolver adamssolver(integrator);

    // Get forward and backwards meeting point
    unsigned int critical_point = lattice->real_to_lattice((1./Z + orbital->Nu())/2.);

    // Integrate backwards until first maximum
    if(exchange)
        hf->SetODEParameters(orbital->Kappa(), orbital->Energy(), exchange.get());
    else
        hf->SetODEParameters(orbital->Kappa(), orbital->Energy());

    // Get classical turning point
    RadialFunction V = hf->GetDirectPotential();
    const double* R = lattice->R();
    int classical_turning_point = V.size();
    double centrifugal = -0.5 * orbital->Kappa() * (orbital->Kappa() + 1);
    double P;
    do
    {   classical_turning_point--;
        P = orbital->Energy() + V.f[classical_turning_point]
            + centrifugal/(R[classical_turning_point] * R[classical_turning_point]);
    }while(P < 0. && classical_turning_point > 0);

    Orbital backwards(*orbital);
    unsigned int peak = adamssolver.IntegrateBackwardsUntilPeak(hf, &backwards, classical_turning_point);
    if(peak < critical_point)
        peak = critical_point;

    double f_right = backwards.f[peak];
    double g_right = backwards.g[peak];

    odesolver->IntegrateForwards(hf, orbital.get());
    
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
    double delta_nu = pow(orbital->Nu(), 3.) * orbital->f[peak] * (tail_scaling * g_right - orbital->g[peak])/norm/alpha;
    double new_nu = orbital->Nu() - delta_nu;

    return -0.5/(new_nu * new_nu) - orbital->Energy();
}

unsigned int HartreeFocker::CalculateContinuumWave(pContinuumWave s, pHFOperator hf)
{
    unsigned int loop = 0;
    double final_amplitude, final_phase;
    pContinuumWave previous_s(new ContinuumWave(s->Kappa()));
    pSpinorFunction exchange(new SpinorFunction(s->Kappa()));
    double ds, old_phase = 0.;
    unsigned int start_sine = 0;
    bool all_done = false;

    pLattice lattice = hf->GetLattice();
    s->resize(lattice->size());
    double pi = MathConstant::Instance()->Pi();

    do
    {   loop++;

        *previous_s = *s;
        start_sine = IntegrateContinuum(s, hf, exchange, final_amplitude, final_phase);
        if(!start_sine)
        {   // Likely reason for not reaching start_sine is that the lattice is too small. Extend it and try again.

            // Probably get a bit worried if we've had, say, fifty complete oscillations and still no good.
            // k ~ sqrt(2(energy + V(r))) => 50 oscillations = 50 * 2 Pi/k
            // Generate a warning and expand the lattice.
            double Vmin = hf->GetDirectPotential().f[lattice->size()-1];
            double kmin = sqrt(2. * (s->Energy() + Vmin));
            double Rmax = 100. * pi/kmin;

            if(lattice->MaxRealDistance() > mmax(Rmax, 1000))
            {   *errstream << "HartreeFocker::CalculateContinuumWave:\n"
                           << "    start_sine not reached; Rmax = " << hf->GetLattice()->MaxRealDistance()
                           << "\n    ...expanding lattice and continuing." << std::endl;
            }

            lattice->resize_to_r(lattice->MaxRealDistance() * 2.);
            s->Clear();
            s->resize(lattice->size());
            old_phase = -1.0;
            start_sine = 0;
            ds = 1.0;   // Do another loop!
        }
        else
        {   final_phase = final_phase/MathConstant::Instance()->Pi();
            final_phase = final_phase - std::floor(final_phase);
            ds = final_phase - old_phase;
            if(fabs(ds) > EnergyTolerance)
            {
                *s *= 0.5;
                *s += (*previous_s * 0.5);
                *exchange = hf->GetExchange(s);
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
    while((loop < MaxHFIterations) && (fabs(ds) > TailMatchingEnergyTolerance));

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
    ContinuumPhaseODE phase_ode(hf->GetLattice(), hf_direct, kappa, energy);
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
