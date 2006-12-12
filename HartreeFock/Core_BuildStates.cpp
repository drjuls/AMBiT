#include "Include.h"
#include <algorithm>
#include "Core.h"
#include "Universal/Constant.h"
#include "Universal/CoulombIntegrator.h"
#include "StateIntegrator.h"
#include "BoundStateIntegrator.h"
#include "GreensIntegrator.h"

#define NON_REL_SCALING true

void Core::Update()
{
    // Check that there is something to start with
    if((Z > 1.) && Empty())
        BuildFirstApproximation();
    UpdateHFPotential();

    bool debug = DebugOptions.LogHFIterations();

    // Copy states for next iteration.
    StateManager next_states(lattice, (unsigned int)Z, (int)Charge);

    StateIterator core_it(this);
    core_it.First();
    while(!core_it.AtEnd())
    {   next_states.AddState(new DiscreteState(*core_it.GetState()));
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
            DiscreteState* new_state = it.GetState();
            double old_energy = new_state->Energy();

            CoupledFunction exchange;
            CalculateExchange(*new_state, exchange);
            deltaE = IterateDiscreteStateGreens(new_state, &exchange);

            if(debug)
                *logstream << "  " << std::setw(4) << new_state->Name()
                           << "  E = " << std::setprecision(12) << old_energy
                           << "  deltaE = " << std::setprecision(4) << deltaE
                           << "  size: (" << new_state->Size()
                           << ") " << lattice->R(new_state->Size()) << std::endl;

            deltaE = fabs(deltaE/new_state->Energy());
            max_deltaE = mmax(deltaE, max_deltaE);

            it.Next();
        }

        // Mix new and old states.
        core_it.First();
        while(!core_it.AtEnd())
        {
            DiscreteState* core_state = core_it.GetState();
            DiscreteState* new_state = next_states.GetState(StateInfo(core_state));

            // Add proportion of new states to core states.
            core_state->Scale(1. - prop_new);
            new_state->Scale(prop_new);
            core_state->ReSize(mmax(core_state->Size(), new_state->Size()));
            new_state->ReSize(mmax(core_state->Size(), new_state->Size()));

            for(unsigned int i = 0; i<core_state->Size(); i++)
            {   core_state->f[i] += new_state->f[i];
                core_state->g[i] += new_state->g[i];
                core_state->df[i] += new_state->df[i];
                core_state->dg[i] += new_state->dg[i];
            }

            // Renormalise core states (should be close already) and update energy.
            core_state->ReNormalise();
            
            double energy = (1. - prop_new) * core_state->Energy() + prop_new * new_state->Energy();
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
}

void Core::UpdateHFPotential(double proportion_new, bool first_build)
{
    HFPotential.resize(lattice->Size());

    // Get electron density function
    std::vector<double> density;
    density.clear();

    ConstStateIterator cs = GetConstStateIterator();
    while(!cs.AtEnd())
    {
        const DiscreteState& s = *cs.GetState();
        const std::vector<double>& f = s.f;
        const std::vector<double>& g = s.g;
        double number_electrons = s.Occupancy();

        if(s.Size() > density.size())
            density.resize(s.Size());

        for(unsigned int j=0; j<s.Size(); j++)
        {
            density[j] = density[j] + (pow(f[j], 2.) + Constant::AlphaSquared * pow(g[j], 2.))*number_electrons;
        }

        cs.Next();
    }

    std::vector<double> y(density.size());
    if(density.size())
    {   CoulombIntegrator I(lattice);
        I.CoulombIntegrate(density, y, 0, Z-Charge);
    }

    // Get local exchange approximation
    std::vector<double> new_local_exchange(HFPotential.size());
    double C = 0.635348143228;
    unsigned int i;
    for(i=0; i<new_local_exchange.size() && i<density.size(); i++)
    {   
        new_local_exchange[i] = C * pow((density[i]/pow(lattice->R(i),2.)), 1./3.);
    }

    std::vector<double> new_potential(HFPotential.size());

    i=0;
    while((i < NuclearPotential.size()) && (i < y.size()))
    {   new_potential[i] = NuclearPotential[i] - y[i];
        i++;
    }
    while(i < NuclearPotential.size())
    {   new_potential[i] = NuclearPotential[i];
        i++;
    }
    while(i < y.size())
    {   new_potential[i] = Z/lattice->R(i) - y[i];
        i++;
    }
    while(i < HFPotential.size())
    {   new_potential[i] = Charge/lattice->R(i);
        i++;
    }

    if(VolumeShiftPotential.size() && VolumeShiftParameter)
    {   for(i=0; i<VolumeShiftPotential.size(); i++)
            new_potential[i] = new_potential[i] + VolumeShiftParameter * VolumeShiftPotential[i];
    }


    if(first_build)
    {
        double HFCharge = Charge;
        if(HFCharge == 0.)
            HFCharge = 1.;

        for(i=0; i<HFPotential.size(); i++)
        {   new_potential[i] = new_potential[i] + new_local_exchange[i];
            if(new_potential[i] < HFCharge/lattice->R(i))
                new_potential[i] = HFCharge/lattice->R(i);
            HFPotential[i] = new_potential[i] * proportion_new + HFPotential[i] * (1. - proportion_new);
        }
    }
    else
    {   for(i=0; i<HFPotential.size(); i++)
            HFPotential[i] = proportion_new * new_potential[i]
                            + (1. - proportion_new) * HFPotential[i];
        LocalExchangeApproximation.resize(new_local_exchange.size());
        for(unsigned int i=0; i<new_local_exchange.size(); i++)
        {   if(new_local_exchange[i] > (HFPotential[i] - Charge/lattice->R(i)))
                LocalExchangeApproximation[i] = HFPotential[i] - Charge/lattice->R(i);
            else
                LocalExchangeApproximation[i] = new_local_exchange[i];
        }
    }
}

void Core::ExtendPotential() const
{
    unsigned int old_size = HFPotential.size();
    HFPotential.resize(lattice->Size());

    for(unsigned int i=old_size; i<lattice->Size(); i++)
        HFPotential[i] = Charge/lattice->R(i);

    LocalExchangeApproximation.resize(HFPotential.size());
}

void Core::UpdateNuclearPotential()
{
    if(NuclearRadius < 0.1/Constant::AtomicToFermi)
    {   // Point nucleus
        NuclearPotential.clear();
    }
    else if(NuclearThickness < 0.1/Constant::AtomicToFermi)
    {   // Uniform hard sphere
        NuclearPotential.resize(lattice->real_to_lattice(NuclearRadius));
        for(unsigned int i=0; i<NuclearPotential.size(); i++)
        {   double radius_ratio = lattice->R(i)/NuclearRadius;
            NuclearPotential[i] = (1.5 - 0.5*radius_ratio*radius_ratio)*Z/NuclearRadius;
        }
    }
    else
    {   // Create nuclear density function
        std::vector<double> density = CalculateNuclearDensity(NuclearRadius, NuclearThickness);

        CoulombIntegrator I(lattice);
        I.CoulombIntegrate(density, NuclearPotential, 0, Z);
    }
}

std::vector<double> Core::CalculateNuclearDensity(double radius, double thickness) const
{
    std::vector<double> density(lattice->Size());
    double B = 4.*log(3.)/thickness;
    double A = 3.*Z/(radius * (pow(radius,2.) + pow(Constant::Pi/B, 2.)));
    for(unsigned int i=0; i<lattice->Size(); i++)
    {   
        double X = B * (lattice->R(i) - radius);
        if(X <= -20.)
            density[i] = A * pow(lattice->R(i), 2.);
        else if(X < 50)
            density[i] = A/(1. + exp(X)) * pow(lattice->R(i), 2.);
        else
        {   density.resize(i);
            break;
        }
    }

    return density;
}

void Core::CalculateClosedShellRadius()
{
    // Get electron density function
    std::vector<double> density;
    density.clear();
    unsigned int j;

    ConstStateIterator cs = GetConstStateIterator();
    while(!cs.AtEnd())
    {
        const DiscreteState& s = *cs.GetState();
        if(OpenShellStates.find(StateInfo(&s)) == OpenShellStates.end())
        {
            const std::vector<double>& f = s.f;
            const std::vector<double>& g = s.g;
            double number_electrons = s.Occupancy();

            if(s.Size() > density.size())
                density.resize(s.Size());

            for(unsigned int j=0; j<s.Size(); j++)
            {
                density[j] = density[j] + (f[j] * f[j] + Constant::AlphaSquared * g[j] * g[j])*number_electrons;
            }
        }

        cs.Next();
    }

    std::vector<double>::const_iterator it = density.end();
    it--;
    double last_max = fabs(*it);
    while(it != density.begin())
    {
        if(fabs(*it) >= last_max)
            last_max = fabs(*it);
        else break;
        it--;
    }

    last_max *= 0.5;
    it = density.end();
    it--;
    j = density.size() - 1;
    while(it != density.begin())
    {
        if(fabs(*it) >= last_max)
            break;
        it--;
        j--;
    }

    ClosedShellRadius = lattice->R(j);
}

unsigned int Core::ConvergeStateApproximation(DiscreteState* s, bool include_exch) const
{
    StateIntegrator I(lattice);
    I.SetAdamsOrder(10);

    if(s->Size() > (I.AdamsOrder()-1))
        s->CheckSize(StateParameters::WavefunctionTolerance);
    else
    {   double r_cutoff = (2. * Charge * (s->Nu()) + 10.) * (s->Nu());
        s->ReSize(lattice->real_to_lattice(r_cutoff));
    }

    if(lattice->Size() > HFPotential.size())
        ExtendPotential();

    std::vector<double> Potential = GetHFPotential();
    if(include_exch)
    {   std::vector<double> LocalExchange = GetLocalExchangeApproximation();
        for(unsigned int i=0; i<Potential.size(); i++)
            Potential[i] = Potential[i] + LocalExchange[i];
    }

    double delta=0.;
    double old_delta1 = 0.;

    unsigned int loop = 0;

    do
    {   loop++;

        // Get forward and backwards meeting point
        unsigned int critical_point = lattice->real_to_lattice((1./Z + s->Nu())/2.);

        DiscreteState no_exchange_state(*s);
        unsigned int peak = I.IntegrateBackwardsUntilPeak(no_exchange_state, Potential, critical_point);

        // assert(peak > 3)
        I.IntegrateBackwards(*s, Potential, NULL, peak-1);
        double f_right = s->f[peak];
        double g_right = s->g[peak];

        I.IntegrateForwards(*s, Potential, NULL, peak+1);

        // Test whether forwards and backwards integrations met and modify energy accordingly.

        // Scale the backwards integration to meet the forwards one in f[peak].
        double tail_scaling = (s->f[peak] - f_right)/no_exchange_state.f[peak];
        unsigned int i;
        for(i = peak+1; i<s->Size(); i++)
        {   s->f[i] = s->f[i] + tail_scaling * no_exchange_state.f[i];
            s->g[i] = s->g[i] + tail_scaling * no_exchange_state.g[i];
        }
        f_right += tail_scaling * no_exchange_state.f[peak];
        g_right += tail_scaling * no_exchange_state.g[peak];

        double norm = s->Norm();
        delta = pow(s->Nu(), 3.) * f_right * (g_right - s->g[peak])/norm;

        if(fabs(delta) > 0.1)
                delta = 0.1 * delta/fabs(delta);

        // Set new principal quantum number
        s->SetNu(s->Nu() - delta);
        s->ReNormalise();

        s->CheckSize(StateParameters::WavefunctionTolerance);
        if(lattice->Size() > HFPotential.size())
        {   // Lengthen potential if necessary. Local exchange is zero outside the core.
            ExtendPotential();
            unsigned int old_size = Potential.size();
            Potential.resize(HFPotential.size());
            for(unsigned int i=old_size; i<HFPotential.size(); i++)
                    Potential[i] = HFPotential[i];
        }
    }
    while((loop < StateParameters::MaxHFIterations) && (fabs(delta) > StateParameters::EnergyTolerance));

    s->ReNormalise();

    return loop;
}

double Core::IterateDiscreteStateGreens(DiscreteState* s, CoupledFunction* exchange) const
{
    StateIntegrator I(lattice);
    I.SetAdamsOrder(10);
    BoundStateIntegrator BI(lattice);
    BI.SetAdamsOrder(10);

    if(s->Size() > (I.AdamsOrder()-1))
        s->CheckSize(StateParameters::WavefunctionTolerance);
    else
    {   double r_cutoff = (2. * Charge * (s->Nu()) + 10.) * (s->Nu());
        s->ReSize(lattice->real_to_lattice(r_cutoff));
    }

    if(lattice->Size() > HFPotential.size())
        ExtendPotential();

    const double* R = lattice->R();
    const double* dR = lattice->dR();

    const std::vector<double>& Potential(HFPotential);

    unsigned int i;
    double delta_E = 0.;

    // Solve Dirac equation
    DiscreteState s0(*s);
    DiscreteState sInf(*s);
    // This routine can modify the size
    I.IntegrateBackwards(sInf, Potential, NULL, -1);

    s0.ReSize(sInf.Size());
    s->ReSize(sInf.Size());
    exchange->ReSize(sInf.Size());
    BI.IntegrateForwards(s0, Potential, NULL, s0.Size());

    GreensIntegrator greens(lattice);
    greens.SetSolutionOrigin(&s0);
    greens.SetSolutionInfinity(&sInf);
    greens.SetGreensIntegrand(exchange);

    std::vector<double> Ginf = greens.GetGreensInfinity();
    std::vector<double> G0 = greens.GetGreensOrigin();

    for(i=0; i<s->Size(); i++)
    {
        s->f[i] = -(s0.f[i] * Ginf[i] + sInf.f[i] * G0[i]);
        s->g[i] = -(s0.g[i] * Ginf[i] + sInf.g[i] * G0[i]);
    }

    // Calculate norm and adjust energy so that norm is closer to unity.
    // E -> E + dE, f -> f + df, g -> g + dg
    double norm = s->Norm();
    double old_energy = s->Energy();

    CoupledFunction del_s(s->Size());  // ds(r)/dE

    greens.SetGreensIntegrand(s);
    Ginf = greens.GetGreensInfinity();
    G0 = greens.GetGreensOrigin();

    for(i=0; i<s->Size(); i++)
    {
        del_s.f[i] = -(s0.f[i] * Ginf[i] + sInf.f[i] * G0[i]);
        del_s.g[i] = -(s0.g[i] * Ginf[i] + sInf.g[i] * G0[i]);
    }

    double var = 0;
    for(i=0; i<s->Size(); i++)
    {   var += (del_s.f[i] * s->f[i] + Constant::AlphaSquared * del_s.g[i] * s->g[i])*dR[i];
    }

    delta_E = (1. - norm)/(2. * var);
    // Limit change in energy to 20%. Place warning in error file.
    if(fabs(delta_E/old_energy) > 0.2)
    {   *errstream << "IterateDiscreteStateGreens: large delta_E: = " << fabs(delta_E/old_energy)*100 << "%" << std::endl;
        delta_E = fabs(old_energy*delta_E*0.2)/delta_E;
    }
    s->SetEnergy(old_energy + delta_E);

    // Solve Dirac equation again with new energy
    s->CheckSize(StateParameters::WavefunctionTolerance);

    s0 = *s;
    sInf = *s;

    I.IntegrateBackwards(sInf, Potential, NULL, -1);

    s0.ReSize(sInf.Size());
    s->ReSize(sInf.Size());
    exchange->ReSize(sInf.Size());
    BI.IntegrateForwards(s0, Potential, NULL, s0.Size());

    greens.SetSolutionOrigin(&s0);
    greens.SetSolutionInfinity(&sInf);
    greens.SetGreensIntegrand(exchange);

    Ginf = greens.GetGreensInfinity();
    G0 = greens.GetGreensOrigin();

    for(i=0; i<s->Size(); i++)
    {
        s->f[i] = -(s0.f[i] * Ginf[i] + sInf.f[i] * G0[i]);
        s->g[i] = -(s0.g[i] * Ginf[i] + sInf.g[i] * G0[i]);

        s->df[i] = (- s->Kappa()/R[i] * s->f[i] +
                    (2. + Constant::AlphaSquared*(s->Energy() + Potential[i])) * s->g[i]
                    + Constant::AlphaSquared * exchange->g[i]) * dR[i];
        s->dg[i] = (s->Kappa()/R[i] * s->g[i] - (s->Energy() + Potential[i]) * s->f[i]
                    - exchange->f[i]) *dR[i];
    }

    s->CheckSize(StateParameters::WavefunctionTolerance);
    return delta_E;
}

unsigned int Core::CalculateExcitedState(State* s) const
{
    // Forward continuum states to other method.
    DiscreteState* ds = dynamic_cast<DiscreteState*>(s);
    if(ds == NULL)
    {   ContinuumState* cs = dynamic_cast<ContinuumState*>(s);
        return CalculateContinuumState(cs);
    }

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
    {   if(!Charge)
        {   *errstream << "Cannot calculate excited Hartree-Fock states: Charge = 0." << std::endl;
            exit(1);
        }

        // Calculate
        // Get a trial value for nu
        double trial_nu;
        switch(ds->L())
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
                trial_nu = (double)ds->L() + 1.;
                break;
        }

        // Increase nu depending on how far above the core it is
        unsigned int largest_core_pqn = 1;
        ConstStateIterator i = GetConstStateIterator();
        while(!i.AtEnd())
        {   
            const DiscreteState* cs = i.GetState();
            if((cs->Kappa() == ds->Kappa()) && (cs->RequiredPQN() > largest_core_pqn))
                    largest_core_pqn = cs->RequiredPQN();
            i.Next();
        }

        trial_nu = trial_nu + (double)(ds->RequiredPQN() - largest_core_pqn - 1)/Charge;
        ds->SetNu(trial_nu/Charge);

        // Calculate wavefunction with local exchange approximation.
        // Check that the number of zeroes of the wavefunction is correct, otherwise
        // adjust nu and start again.

        int zero_difference = 0;    // Difference between required and actual number of zeroes of wavefunction

        do
        {   loop = ConvergeStateApproximation(ds);

            if(loop >= StateParameters::MaxHFIterations)
            {   // Evil - this should never happen given our simple model.
                *errstream << "CalculateExcitedState: Local exchange approximation failed to converge.\n"
                           << "  " << ds->Name() << "  iterations = " << loop << std::endl;
                PAUSE
                exit(1);
            }

            zero_difference = ds->NumZeroes() + ds->L() + 1 - ds->RequiredPQN();

            if(zero_difference)
            {   *logstream << "    " << ds->Name()
                           << std::setprecision(7) << "  E = " << ds->Energy()
                           << "  loops: " << loop  << "  size: " << ds->Size()
                           << "  zerodiff: " << zero_difference << std::endl;

                // This moves the state one pqn at a time
                trial_nu = ds->Nu() - zero_difference/abs(zero_difference)/Charge;
                ds->SetNu(trial_nu);
            }
        }
        while(zero_difference || (loop >= StateParameters::MaxHFIterations));

        // Hartree-Fock loops
        bool debugHF = DebugOptions.LogHFIterations();

        double deltaE;
        loop = 0;
        DiscreteState* new_ds = new DiscreteState(*ds);
        double prop_new = 0.4;
        do
        {   loop++;
            CoupledFunction exchange;
            CalculateExchange(*new_ds, exchange);
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
            ds->ReNormalise();
            ds->SetEnergy((1. - prop_new) * ds->Energy() + prop_new * new_ds->Energy());
            *new_ds = *ds;
            deltaE = fabs(deltaE/new_ds->Energy());

        }while((deltaE > StateParameters::EnergyTolerance) && (loop < StateParameters::MaxHFIterations));

        delete new_ds;
        if(loop >= StateParameters::MaxHFIterations)
            *errstream << "Core: Failed to converge excited HF state " << ds->Name() << std::endl;
    }

    if(DebugOptions.OutputHFExcited())
    {
        *outstream << std::setprecision(12);
        if(DebugOptions.HartreeEnergyUnits() || DebugOptions.InvCmEnergyUnits())
        {
            double energy = ds->Energy();
            if(DebugOptions.InvCmEnergyUnits())
                energy *= Constant::HartreeEnergy_cm;
            *outstream << ds->Name() << "  E = " << energy;
        }
        else
            *outstream << ds->Name() << "  nu = " << ds->Nu();

        *outstream << "  loops: " << loop << "  size: (" << ds->Size() << ") " << lattice->R(ds->Size()) << std::endl;
    }

    return loop;
}

unsigned int Core::CalculateContinuumState(ContinuumState* s) const
{
    unsigned int loop = 0;
    double final_amplitude, final_phase;
    CoupledFunction exchange(HFPotential.size());
    CoupledFunction new_exchange(HFPotential.size());
    double ds, old_phase = 0.;
    unsigned int start_sine = 0;

    s->ReSize(HFPotential.size());

    do
    {   loop++;
        
        StateIntegrator I(lattice);
        start_sine = I.IntegrateContinuum(*s, HFPotential, exchange, Z, 0.01, final_amplitude, final_phase);

        ds = (final_phase - old_phase)/Constant::Pi;
        if(fabs(ds) > StateParameters::EnergyTolerance)
        {   CalculateExchange(*s, new_exchange);
            exchange.ReSize(new_exchange.Size());
            for(unsigned int i = 0; i<exchange.Size(); i++)
                exchange.f[i] = 0.5 * exchange.f[i] + 0.5 * new_exchange.f[i];
            old_phase = final_phase;
        }
    }
    while((loop < StateParameters::MaxHFIterations) && (fabs(ds) > StateParameters::EnergyTolerance));

    final_amplitude = sqrt(2./(Constant::Pi*pow(s->Nu(),3.)))/final_amplitude;
    s->Scale(final_amplitude);

    //Graph g;
    //g.SetData(s->f, "f");
    //g.View();
    //PAUSE
    //g.CloseView();

    return loop;
}

void Core::CalculateExchange(const State& current, CoupledFunction& exchange, const SigmaPotential* sigma, double sigma_amount) const
{
    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);

    exchange.Clear();
    exchange.ReSize(current.Size());

    const DiscreteState* ds_current = dynamic_cast<const DiscreteState*>(&current);

    // Sum over all core states
    ConstStateIterator cs = GetConstStateIterator();
    while(!cs.AtEnd())
    {
        const DiscreteState& other = *(cs.GetState());
        unsigned int upper_point = mmin(other.Size(), current.Size());

        // Get overlap of wavefunctions
        std::vector<double> density(upper_point);
        for(unsigned int i=0; i<upper_point; i++)
        {
            density[i] = current.f[i] * other.f[i] + Constant::AlphaSquared * current.g[i] * other.g[i];
        }

        // Sum over all k
        for(unsigned int k = abs((int)other.L() - (int)current.L()); k <= (other.L() + current.L()); k+=2)
        {
            double coefficient = Constant::Wigner3j(k, current.J(), other.J(), 0., .5, -.5);
            coefficient = (2 * abs(other.Kappa())) * coefficient * coefficient;

            // Open shells need to be scaled
            if(IsOpenShellState(StateInfo(&other)) && (other.Occupancy() != double(2 * abs(other.Kappa()))))
            {
                double ex = 1.;
                if(NON_REL_SCALING)
                {   // Average over non-relativistic configurations
                    if(other.Kappa() == -1)
                    {
                        if((ds_current == NULL) || (StateInfo(ds_current) != StateInfo(&other)))
                            ex = other.Occupancy()/double(2 * abs(other.Kappa()));
                        else if(k)
                            ex = (other.Occupancy()-1.)/double(2 * abs(other.Kappa()) - 1);
                    }
                    else
                    {
                        int other_kappa = - other.Kappa() - 1;
                        const DiscreteState* ds = GetState(StateInfo(other.RequiredPQN(), other_kappa));

                        if((ds_current == NULL) || ((StateInfo(ds_current) != StateInfo(&other)) && (StateInfo(ds_current) != StateInfo(ds))))
                            ex = (other.Occupancy() + ds->Occupancy())/double(2 * (abs(other.Kappa()) + abs(ds->Kappa())));
                        else if(k)
                            ex = (other.Occupancy() + ds->Occupancy() - 1.)/double(2 * (abs(other.Kappa()) + abs(ds->Kappa())) - 1);
                    }
                }
                else
                {   // Average over relativistic configurations
                    if((ds_current == NULL) || (StateInfo(ds_current) != StateInfo(&other)))
                        ex = other.Occupancy()/double(2 * (abs(other.Kappa())));
                    else if(k)
                        ex = (other.Occupancy() - 1.)/double(2 * (abs(other.Kappa())) - 1);
                }

                coefficient = coefficient * ex;
            }

            // Integrate density to get (1/r)Y(ab,r)
            std::vector<double> potential;
            I.CoulombIntegrate(density, potential, k);

            for(unsigned int i=0; i<upper_point; i++)
            {
                exchange.f[i] = exchange.f[i] + coefficient * potential[i] * other.f[i];
                exchange.g[i] = exchange.g[i] + coefficient * potential[i] * other.g[i];
            }

            if(NuclearInverseMass && (k == 1))
            {
                std::vector<double> P(upper_point);
                double sms = SI.IsotopeShiftIntegral(current, other, &P);
                
                for(unsigned int i=0; i<upper_point; i++)
                {
                    exchange.f[i] = exchange.f[i] + coefficient * NuclearInverseMass * sms *  P[i];
                }
            }
        }
        cs.Next();
    }

    if(sigma && sigma_amount)
    {
        std::vector<double> potential = sigma->GetPotential(current.f);
        for(unsigned int i=0; i<mmin(current.Size(), (unsigned int)potential.size()); i++)
            exchange.f[i] = exchange.f[i] - sigma_amount * potential[i];
    }
}

const unsigned int Core::StateParameters::MaxHFIterations = 300;
double Core::StateParameters::WavefunctionTolerance = 1.E-11;
double Core::StateParameters::EnergyTolerance = 1.E-14;
