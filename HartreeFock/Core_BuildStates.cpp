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
    if(Empty())
        BuildFirstApproximation();

    bool debug = DebugOptions.LogHFIterations();

    double delta;
    unsigned int loop = 0;
    UpdateHFPotential();
    StateParameters::WavefunctionTolerance = 1.e-5;

    bool do_greens = true;
    if(AllStates.size() == 1)
        do_greens = false;

    do
    {   delta = 0.;
        loop++;
        if(debug)
            *logstream << "HF Iteration :" << loop << std::endl;
        DiscreteStateIterator i = GetDiscreteStateIterator();
        while(!i.AtEnd())
        {
            // Important to make a copy of the state so that the exchange 
            //  interaction doesn't get messed up.
            DiscreteState* s = new DiscreteState(*(i.GetState()));
            double old_nu = s->Nu();
            unsigned int l = s->L();
            unsigned int pqn = s->RequiredPQN();
            UpdateHFPotentialNoSelfInteraction(StateInfo(s));

            unsigned int iterations = CalculateDiscreteState(s);
            if(iterations < StateParameters::MaxHFIterations)
            {
                i.ReplaceState(s);
            }
            else
            {   *errstream << s->Name() << " failed to converge after " << StateParameters::MaxHFIterations << " iterations. " << std::endl;
                delete s;
                s = i.GetState();
            }
            double nu = s->Nu();

            delta = mmax(delta, fabs((nu - old_nu)/nu));
            if(debug)
            {   if(DebugOptions.HartreeEnergyUnits())
                    *logstream << "  " << s->Name() << "  en: " << -0.5/(nu*nu) << "  loops: " << iterations << "  size: " << s->Size() << std::endl;
                else
                    *logstream << "  " << s->Name() << "  nu: " << nu << "  loops: " << iterations << "  size: " << s->Size() << std::endl;
            }
            i.Next();
        }
        UpdateHFPotential();
    }
    while(!(do_greens && (loop >= 10)) && (delta > StateParameters::EnergyTolerance));

    if(do_greens)
    {   double tolerance = StateParameters::EnergyTolerance;
        bool converged = false;
        bool converged_once = false;
        while(!converged)
        {   converged = UpdateGreens();
            if(!converged)
            {   StateParameters::EnergyTolerance *= 10.;
                *outstream << "\n Tolerance = " << StateParameters::EnergyTolerance << std::endl;
            }
            else if(tolerance < StateParameters::EnergyTolerance)
            {   if(!converged_once)
                {   converged_once = true;
                    StateParameters::EnergyTolerance = StateParameters::EnergyTolerance/10.;
                    *outstream << "\n Tolerance = " << StateParameters::EnergyTolerance << std::endl;
                    converged = false;
                }
            }
        }
        StateParameters::EnergyTolerance = tolerance;
    }
}

bool Core::UpdateGreens()
{
    bool debug = DebugOptions.LogHFIterations();

    // Initial iterations for deltaE
    if(debug)
        *logstream << "\n Greens function iterations:" << std::endl;

    DiscreteStateIterator i = GetDiscreteStateIterator();
    std::map<StateInfo, double> deltaE;
    CoupledFunction exchange;

    while(!i.AtEnd())
    {
        DiscreteState* s = new DiscreteState(*(i.GetState()));

        UpdateHFPotentialNoSelfInteraction(StateInfo(s));
        CalculateExchangeNoSelfInteraction(*s, exchange);

        deltaE[StateInfo(s)] = IterateDiscreteStateGreens(s, &exchange)/s->Energy();
        if(debug)
            *logstream << "  " << s->Name() << "  energy: " << s->Energy() << "  deltaE: " << deltaE[StateInfo(s)] << "  norm: " << s->Norm()-1. << std::endl;

        i.ReplaceState(s);
        i.Next();
    }

    // Further iterations until full convergence
    StateParameters::WavefunctionTolerance = 1.e-7;
    double max_deltaE = 1.;
    unsigned int iterations = 0;

    while(max_deltaE > StateParameters::EnergyTolerance && iterations < StateParameters::MaxHFIterations)
    {
        max_deltaE = 0.;
        std::map<StateInfo, double>::iterator it = deltaE.begin();
        std::map<StateInfo, double>::iterator max_it = deltaE.begin();
        while(it != deltaE.end())
        {   if(fabs(it->second) > max_deltaE)
            {   max_deltaE = fabs(it->second);
                max_it = it;
            }
            it++;
        }

        DiscreteState* s = new DiscreteState(*GetState(max_it->first));

        UpdateHFPotentialNoSelfInteraction(max_it->first);
        CalculateExchangeNoSelfInteraction(*s, exchange);
        max_it->second = IterateDiscreteStateGreens(s, &exchange)/s->Energy();

        StateSet::iterator coreit = AllStates.find(max_it->first);
        coreit->second.DeleteState();
        coreit->second.SetState(s);

        if(debug)
            *logstream << "  " << s->Name() << "  energy: " << s->Energy() << "  deltaE: " << deltaE[StateInfo(s)] << "  norm: " << s->Norm()-1. << std::endl;

        iterations++;
    }

    if(iterations < StateParameters::MaxHFIterations)
    {   if(debug)
        {   *logstream << "\n Final values:" << std::endl;

            i.First();
            while(!i.AtEnd())
            {
                DiscreteState* s = new DiscreteState(*(i.GetState()));
                *logstream << "  " << s->Name() << "  en: " << std::setprecision(12) << s->Energy() << std::endl;
                i.Next();
            }
        }

        UpdateHFPotential();

        if(debug)
            *logstream << "Orthogonality test: " << TestOrthogonality() << std::endl;
        return true;
    }
    else
        return false;
}

void Core::UpdateHFPotential(double proportion_new, bool first_build)
{
    double HFCharge = Charge;
    if(HFCharge == 0. && first_build)
        HFCharge = 1.;
    HFPotential.resize(lattice->Size());

    // Get electron density function
    std::vector<double> density;
    density.clear();

    ConstDiscreteStateIterator cs = GetConstDiscreteStateIterator();
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
    CoulombIntegrator I(*lattice);
    I.CoulombIntegrate(density, y, 0, Z-Charge);

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
        {   if(new_local_exchange[i] > (HFPotential[i] - HFCharge/lattice->R(i)))
                LocalExchangeApproximation[i] = HFPotential[i] - HFCharge/lattice->R(i);
            else
                LocalExchangeApproximation[i] = new_local_exchange[i];
        }
    }
}

void Core::UpdateHFPotentialNoSelfInteraction(const StateInfo& current, double proportion_new)
{
    double HFCharge = Charge;
    HFPotential.resize(lattice->Size());

    // Get electron density function
    std::vector<double> density;
    density.clear();

    ConstDiscreteStateIterator cs = GetConstDiscreteStateIterator();
    while(!cs.AtEnd())
    {
        const DiscreteState& s = *cs.GetState();
        const std::vector<double>& f = s.f;
        const std::vector<double>& g = s.g;
        double number_electrons = s.Occupancy();

        // Open shells need to be scaled
        if(IsOpenShellState(StateInfo(&s)) && (number_electrons != double(2 * abs(s.Kappa()))))
        {   // Average over non-relativistic configurations
            if(s.Kappa() == -1)
            {   if(StateInfo(&s) == current)
                    number_electrons--;
            }
            else
            {   int other_kappa = - s.Kappa() - 1;
                StateInfo other_info(s.RequiredPQN(), other_kappa);

                if((current == StateInfo(&s)) || (current == other_info))
                    number_electrons = number_electrons - double(abs(s.Kappa()))/double(abs(s.Kappa())+abs(other_kappa));
            }
        }
        else if(StateInfo(&s) == current)
            number_electrons--;

        if(s.Size() > density.size())
            density.resize(s.Size());

        for(unsigned int j=0; j<s.Size(); j++)
        {
            density[j] = density[j] + (pow(f[j], 2.) + Constant::AlphaSquared * pow(g[j], 2.))*number_electrons;
        }

        cs.Next();
    }

    std::vector<double> y(density.size());
    CoulombIntegrator I(*lattice);
    I.CoulombIntegrate(density, y, 0, Z-Charge-1.);

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

    for(i=0; i<HFPotential.size(); i++)
        HFPotential[i] = proportion_new * new_potential[i]
                        + (1. - proportion_new) * HFPotential[i];

    LocalExchangeApproximation.resize(new_local_exchange.size());
    for(unsigned int i=0; i<new_local_exchange.size(); i++)
    {   if(new_local_exchange[i] > (HFPotential[i] - HFCharge/lattice->R(i)))
            LocalExchangeApproximation[i] = HFPotential[i] - HFCharge/lattice->R(i);
        else
            LocalExchangeApproximation[i] = new_local_exchange[i];
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
    // Create nuclear density function
    std::vector<double> density = CalculateNuclearDensity(NuclearRadius, NuclearThickness);

    CoulombIntegrator I(*lattice);
    I.CoulombIntegrate(density, NuclearPotential, 0, Z);
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

unsigned int Core::CalculateDiscreteState(DiscreteState* s, double exchange_amount, const SigmaPotential* sigma, double sigma_amount) const
{
    bool remove_self_interaction = false;
    if(GetState(StateInfo(s)))
        remove_self_interaction = true;

    StateIntegrator I(*lattice);
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

    CoupledFunction* exchange = new CoupledFunction();
    CoupledFunction* new_exchange = NULL;
    CoupledFunction* old_exchange = NULL;
    
    bool local_exchange = (exchange_amount < 0.);
    if(local_exchange)
    {
        std::vector<double> LocalExchange = GetLocalExchangeApproximation();
        for(unsigned int i=0; i<Potential.size(); i++)
            Potential[i] = Potential[i] + LocalExchange[i];
        exchange->Clear();
        exchange->ReSize(s->Size());
        exchange_amount = 1.0;
    }
    else if(exchange_amount > 0.)
    {   if(remove_self_interaction)
            CalculateExchangeNoSelfInteraction(*s, *exchange);
        else
            CalculateExchange(*s, *exchange, sigma, sigma_amount);

        exchange->Scale(exchange_amount);
    }
    else
    {   exchange->Clear();
        exchange->ReSize(s->Size());
    }

    double delta1=0.,
           delta2=0.,
           delta=0.;
    double old_delta1 = 0.;

    unsigned int loop = 0;

    do
    {   loop++;

        // Get forward and backwards meeting point
        unsigned int critical_point = lattice->real_to_lattice((1./Z + s->Nu())/2.);
        // assert (critical_point < lattice->NumPoints())
        DiscreteState no_exchange_state(*s);
        unsigned int peak = I.IntegrateBackwardsUntilPeak(no_exchange_state, Potential, critical_point);

        // assert(peak > 3)
        I.IntegrateBackwards(*s, Potential, *exchange, peak-1);
        double f_right = s->f[peak];
        double g_right = s->g[peak];

        I.IntegrateForwards(*s, Potential, *exchange, peak+1);

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
        delta1 = pow(s->Nu(), 3.) * f_right * (g_right - s->g[peak])/norm;
        
        new_exchange = new CoupledFunction();
        if(exchange_amount > 0.)
        {   if(remove_self_interaction)
                CalculateExchangeNoSelfInteraction(*s, *new_exchange);
            else
                CalculateExchange(*s, *new_exchange, sigma, sigma_amount);

            new_exchange->Scale(exchange_amount);
        }
        else
        {   new_exchange->Clear();
            new_exchange->ReSize(s->Size());
        }

        // Integrate overlap of wavefunction with change in exchange potential
        double overlap = 0.;
        for(i = 0; i<s->Size(); i++)
        {   overlap = overlap + s->f[i] * (exchange->f[i] - new_exchange->f[i]) * lattice->dR(i);
        }

        delta2 = -pow(s->Nu(), 3.) * overlap/norm;

        delta = delta1 + delta2;
        if(fabs(delta) > 0.1)
                delta = 0.1 * delta/fabs(delta);

        // Set new principal quantum number
        s->SetNu(s->Nu() - delta);

        double new_exchange_proportion;
        if(fabs(old_delta1 - delta1) >= 1.E-6)
        {   new_exchange_proportion = old_delta1/(old_delta1 - delta1);
            if((new_exchange_proportion <= 0.) || (new_exchange_proportion > 1.))
                new_exchange_proportion = 1.0;
        }
        else
            new_exchange_proportion = 1.0;

        old_delta1 = delta1;

        if((new_exchange_proportion == 1.0) || (old_exchange == NULL))
        {   delete exchange;
            if(old_exchange && (exchange != old_exchange))
                delete old_exchange;
            exchange = new_exchange;
        }
        else
        {   for(unsigned int i = 0; i < old_exchange->Size(); i++)
            {   old_exchange->f[i] = new_exchange_proportion * new_exchange->f[i]
                                    + (1. - new_exchange_proportion) * old_exchange->f[i];
                old_exchange->g[i] = new_exchange->g[i];
            }
            if(exchange != old_exchange)
                delete exchange;
            exchange = old_exchange;
        }
        old_exchange = new_exchange;

        s->CheckSize(StateParameters::WavefunctionTolerance);
        if(lattice->Size() > HFPotential.size())
        {   ExtendPotential();
            Potential = GetHFPotential();
        }

        exchange->ReSize(s->Size());
        old_exchange->ReSize(s->Size());

        if(local_exchange)
        {   
            Potential = GetHFPotential();
            local_exchange = false;
        }
    }
    while((loop < StateParameters::MaxHFIterations) && 
         ((fabs(delta) > StateParameters::EnergyTolerance) || (fabs(delta1) > StateParameters::EnergyTolerance) || (fabs(delta2) > StateParameters::EnergyTolerance)));

    s->ReNormalise();

    if(old_exchange != exchange)
        delete old_exchange;
    delete exchange;

    return loop;
}

double Core::IterateDiscreteStateGreens(DiscreteState* s, CoupledFunction* exchange) const
{
    StateIntegrator I(*lattice);
    I.SetAdamsOrder(10);
    BoundStateIntegrator BI(*lattice);
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

    std::vector<double> Potential = GetHFPotential();

    unsigned int i;
    double delta_E = 0.;

    // Solve Dirac equation
    DiscreteState s0(*s);
    DiscreteState sInf(*s);
    // This routine can modify the size
    I.IntegrateBackwards(sInf, Potential, NULL, -1);

    s0.ReSize(sInf.Size());
    s->ReSize(sInf.Size());
    BI.IntegrateForwards(s0, Potential, NULL, s0.Size());

    s0.ReNormalise();
    sInf.ReNormalise();

    GreensIntegrator greens(*lattice);
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
    s->SetEnergy(old_energy + delta_E);

    // Solve Dirac equation again with new energy
    s->CheckSize(StateParameters::WavefunctionTolerance);
    exchange->ReSize(s->Size());

    s0 = *s;
    sInf = *s;

    I.IntegrateBackwards(sInf, Potential, NULL, -1);

    s0.ReSize(sInf.Size());
    s->ReSize(sInf.Size());
    BI.IntegrateForwards(s0, Potential, NULL, s0.Size());

    s0.ReNormalise();
    sInf.ReNormalise();

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

    unsigned int iterations = 0;
    const DiscreteState* core_state = GetState(StateInfo(s));
    StateSet::const_iterator it = OpenShellStorage.find(StateInfo(s));
    if(core_state != NULL)
    {   // Try to find in core (probably open shells).
        *ds = *core_state;
    }
    else if(it != OpenShellStorage.end())
    {   // Try to find in unoccupied (but previously calculated) open shells.
        *ds = *dynamic_cast<const DiscreteState*>(it->second.GetState());
    }
    else
    {   // Calculate
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
        ConstDiscreteStateIterator i = GetConstDiscreteStateIterator();
        while(!i.AtEnd())
        {   
            const DiscreteState* cs = i.GetState();
            if((cs->Kappa() == ds->Kappa()) && (cs->RequiredPQN() > largest_core_pqn))
                    largest_core_pqn = cs->RequiredPQN();
            i.Next();
        }

        trial_nu = trial_nu + (double)(ds->RequiredPQN() - largest_core_pqn - 1)/Charge;
        ds->SetNu(trial_nu/Charge);

        double nu_change_factor = 0.3;
        int zero_difference = 0,        // Difference between required and actual number of zeroes of wavefunction
            old_zero_difference = 0;
        double old_trial_nu = trial_nu;

        /* There are two methods for adding exchange:
            1. Add it in a little at a time, building up until all of it is included
            2. Use a local approximation to the exchange potential and then include exchange
            properly the next time round. Faster, but relies on the local approximation
            being "close" to the real exchange.
        exchange_addition_method is the first option
        */
        bool exchange_addition_method = false;

        do
        {   ds->Clear();

            if(exchange_addition_method)
            {   // Calculate state, adding exchange a little at a time
                *logstream << "Attempting exchange addition method. Trial nu = " << ds->Nu() << std::endl;

                double exchange_iteration_amount = 0.01;
                double previous_nu = ds->Nu();
                unsigned int small_add_count = 0;
                for(double exch = 0.; exch <= 1.; exch = exch + exchange_iteration_amount)
                {   
                    iterations = CalculateDiscreteState(ds, exch);
                    if(iterations >= StateParameters::MaxHFIterations)
                    {   exch = exch - exchange_iteration_amount;
                        exchange_iteration_amount = exchange_iteration_amount/2.;
                        ds->SetNu(previous_nu);
                        small_add_count = 30;
                    }
                    else
                    {   if(exchange_iteration_amount < 0.01)
                        {   if(small_add_count)
                                small_add_count--;
                            else
                            {   exchange_iteration_amount = exchange_iteration_amount * 2;
                                if(exchange_iteration_amount < 0.01)
                                    small_add_count = 10;
                            }
                        }
                        previous_nu = ds->Nu();
                    }
                }
                iterations = CalculateDiscreteState(ds, 1.);  // Make sure everything is included
            }
            else
            {   // Attempt the local exchange approximation method
                iterations = CalculateDiscreteState(ds, -1.);
                if(iterations >= StateParameters::MaxHFIterations)
                {   // Doesn't work, try other method
                    exchange_addition_method = true;
                }
            }

            if(iterations < StateParameters::MaxHFIterations)
            {   // Check that the number of zeroes of the wavefunction is correct, otherwise
                // adjust nu and start again.
                zero_difference = ds->NumZeroes() + ds->L() + 1 - ds->RequiredPQN();

                if(zero_difference)
                {   
                    // This checks to see whether we are getting anywhere, otherwise we need to
                    // increase the rate of change in nu.
                    if(old_zero_difference && (zero_difference == old_zero_difference))
                    {   trial_nu = old_trial_nu;
                        nu_change_factor = nu_change_factor * 2;
                    }
                    else
                    {   old_trial_nu = trial_nu;
                        old_zero_difference = zero_difference;
                    }

                    ds->SetNu(trial_nu - zero_difference/abs(zero_difference) * nu_change_factor * trial_nu);
                    trial_nu = ds->Nu();
                    nu_change_factor = nu_change_factor * 0.75;
                }
            }
        }
        while(zero_difference || (iterations >= StateParameters::MaxHFIterations));
    }

    if(DebugOptions.OutputHFExcited())
    {
        double energy = ds->Energy();
        if(DebugOptions.HartreeEnergyUnits() || DebugOptions.InvCmEnergyUnits())
        {
            if(DebugOptions.InvCmEnergyUnits())
                energy *= Constant::HartreeEnergy_cm;
            *outstream << ds->Name() << "  energy: " << energy << "  loops: " << iterations << "  size: " << ds->Size() << std::endl;
        }
        else
            *outstream << ds->Name() << "  nu: " << ds->Nu() << "  loops: " << iterations << "  size: " << ds->Size() << std::endl;
    }

    return iterations;
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
        
        StateIntegrator I(*lattice);
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
    CoulombIntegrator I(*lattice);

    exchange.Clear();
    exchange.ReSize(current.Size());

    // Sum over all core states
    ConstDiscreteStateIterator cs = GetConstDiscreteStateIterator();
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
                        if(StateInfo(&current) != StateInfo(&other))
                            ex = other.Occupancy()/double(2 * abs(other.Kappa()));
                        else if(k)
                            ex = (other.Occupancy()-1.)/double(2 * abs(other.Kappa()) - 1);
                    }
                    else
                    {
                        int other_kappa = - other.Kappa() - 1;
                        const DiscreteState* ds = GetState(StateInfo(other.RequiredPQN(), other_kappa));

                        if((StateInfo(&current) != StateInfo(&other)) && (StateInfo(&current) != StateInfo(ds)))
                            ex = (other.Occupancy() + ds->Occupancy())/double(2 * (abs(other.Kappa()) + abs(ds->Kappa())));
                        else if(k)
                            ex = (other.Occupancy() + ds->Occupancy() - 1.)/double(2 * (abs(other.Kappa()) + abs(ds->Kappa())) - 1);
                    }
                }
                else
                {   // Average over relativistic configurations
                    if(StateInfo(&current) != StateInfo(&other))
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
                double sms = I.IsotopeShiftIntegral(current, other, &P);
                
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

void Core::CalculateExchangeNoSelfInteraction(const State& current, CoupledFunction& exchange) const
{
    CoulombIntegrator I(*lattice);

    exchange.Clear();
    exchange.ReSize(current.Size());

    // Sum over all core states
    ConstDiscreteStateIterator cs = GetConstDiscreteStateIterator();
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
        unsigned int k_min;
        if(StateInfo(&current) == StateInfo(&other))
            k_min = 2;
        else
            k_min = abs((int)other.L() - (int)current.L());
            
        for(unsigned int k = k_min; k <= (other.L() + current.L()); k+=2)
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
                        if(StateInfo(&current) != StateInfo(&other))
                            ex = other.Occupancy()/double(2 * abs(other.Kappa()));
                        else if(k)
                            ex = (other.Occupancy()-1.)/double(2 * abs(other.Kappa()) - 1);
                    }
                    else
                    {
                        int other_kappa = - other.Kappa() - 1;
                        const DiscreteState* ds = GetState(StateInfo(other.RequiredPQN(), other_kappa));

                        if((StateInfo(&current) != StateInfo(&other)) && (StateInfo(&current) != StateInfo(ds)))
                            ex = (other.Occupancy() + ds->Occupancy())/double(2 * (abs(other.Kappa()) + abs(ds->Kappa())));
                        else if(k)
                            ex = (other.Occupancy() + ds->Occupancy() - 1.)/double(2 * (abs(other.Kappa()) + abs(ds->Kappa())) - 1);
                    }
                }
                else
                {   // Average over relativistic configurations
                    if(StateInfo(&current) != StateInfo(&other))
                        ex = other.Occupancy()/double(2 * (abs(other.Kappa())));
                    else if(k)
                        ex = (other.Occupancy() - 1.)/double(2 * (abs(other.Kappa())) - 1);
                }

                coefficient = coefficient * ex;
            }

            if(coefficient)
            {
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
                    double sms = I.IsotopeShiftIntegral(current, other, &P);
                    
                    for(unsigned int i=0; i<upper_point; i++)
                    {
                        exchange.f[i] = exchange.f[i] + coefficient * NuclearInverseMass * sms *  P[i];
                    }
                }
            }
        }
        cs.Next();
    }
}

const unsigned int Core::StateParameters::MaxHFIterations = 200;
double Core::StateParameters::WavefunctionTolerance = 1.E-5;
double Core::StateParameters::EnergyTolerance = 1.E-9;
