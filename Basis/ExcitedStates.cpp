#include "Include.h"
#include "ExcitedStates.h"
#include "HartreeFock/Core.h"
#include "Universal/Constant.h"
#include "MBPT/CoreMBPTCalculator.h"
#include "HartreeFock/StateIntegrator.h"
#include <stdio.h>

ExcitedStates::ExcitedStates(Lattice* lattice, const Core* atom_core):
    StateManager(lattice, (unsigned int)atom_core->GetZ(), (unsigned int)atom_core->GetCharge()),
    core(atom_core)
{}

ExcitedStates::~ExcitedStates()
{
    Clear();

    SigmaMap::iterator sigma = SecondOrderSigma.begin();
    while(sigma != SecondOrderSigma.end())
    {   delete sigma->second;
        sigma++;
    }
}

void ExcitedStates::AddState(DiscreteState* s)
{
    StateManager::AddState(s);
}

DiscreteState ExcitedStates::GetStateWithSigma(const StateInfo& info) const
{
    const DiscreteState* s = GetState(info);

    if(s != NULL)
    {   DiscreteState ds(*s);

        SigmaMap::const_iterator it = SecondOrderSigma.find(ds.Kappa());
        if((it != SecondOrderSigma.end()) && GetSigmaAmount(info))
            core->UpdateExcitedState(&ds, it->second, GetSigmaAmount(info));

        return ds;
    }
    else
        return DiscreteState();
}

void ExcitedStates::ClearSigmas()
{
    SigmaMap::iterator sigma = SecondOrderSigma.begin();
    while(sigma != SecondOrderSigma.end())
    {   delete sigma->second;
        sigma++;
    }
    SecondOrderSigma.clear();

    SecondOrderAmount.clear();
}

double ExcitedStates::GetSigmaMatrixElement(const StateInfo& info) const
{
    double matrix_element = 0.;
    const DiscreteState* s = GetState(info);
    const SigmaPotential* sigma = NULL;
    SigmaMap::const_iterator it = SecondOrderSigma.find(info.Kappa());
    if(it != SecondOrderSigma.end())
        sigma = it->second;

    if(s && sigma)
    {   matrix_element = sigma->GetMatrixElement(s->f, s->f);
    }

    return matrix_element;
}

double ExcitedStates::CreateSecondOrderSigma(const StateInfo& info, const CoreMBPTCalculator& mbpt)
{
    DiscreteState* s = GetState(info);
    if(s == NULL)
    {   *errstream << "CreateSecondOrderSigma: " << info.Name() << " is not part of ExcitedStates." << std::endl;
        exit(1);
    }

    SigmaPotential* sigma;

    std::string sigma_file = identifier + "." + itoa(info.Kappa()) + ".sigma";

    SigmaMap::iterator it = SecondOrderSigma.find(s->Kappa());
    if(it == SecondOrderSigma.end())
    {   sigma = new SigmaPotential(lattice, sigma_file, s->Size(), 100);
        sigma->Reset();
        mbpt.GetSecondOrderSigma(s->Kappa(), sigma);
        SecondOrderSigma[s->Kappa()] = sigma;
        sigma->Store();
    }
    else
    {   sigma = it->second;
    }

    double energy = sigma->GetMatrixElement(s->f, s->f);
    SetSigmaAmount(info, 1.);

    return energy;
}

bool ExcitedStates::RetrieveSecondOrderSigma(const StateInfo& info)
{
    std::string sigma_file = identifier + "." + itoa(info.Kappa()) + ".sigma";
    SigmaPotential* sigma = new SigmaPotential(lattice, sigma_file);

    if(sigma->Size())
    {
        SigmaMap::iterator it = SecondOrderSigma.find(info.Kappa());
        if(it != SecondOrderSigma.end())
        {   delete SecondOrderSigma[info.Kappa()];
            SecondOrderSigma.erase(it);
        }

        SecondOrderSigma[info.Kappa()] = sigma;
        SetSigmaAmount(info, 1.);

        return true;
    }
    else
    {   delete sigma;
        return false;
    }
}

void ExcitedStates::SetEnergyViaSigma(const StateInfo& info, double energy)
{
    DiscreteState* s = GetState(info);
    if(s == NULL)
    {   s = new DiscreteState(info.PQN(), info.Kappa());
        core->CalculateExcitedState(s);
    }

    DiscreteState* ds = s;

    SigmaPotential* sigma;
    if(SecondOrderSigma.find(s->Kappa()) == SecondOrderSigma.end())
    {   *errstream << "SetEnergyViaSigma: Sigma with kappa = " << s->Kappa() << " not found." << std::endl;
        exit(1);
    }
    sigma = SecondOrderSigma[s->Kappa()];

    double old_energy = s->Energy();
    double amount = 1.0;
    DiscreteState sigma_s(*ds);
    unsigned int iterations = core->UpdateExcitedState(&sigma_s, sigma, amount);
    double current_energy = sigma_s.Energy();

    *outstream << "  Wanted energy:  " << std::setprecision(12) << energy * Constant::HartreeEnergy_cm << std::endl;
    *outstream << "  Sigma iterated: " << current_energy * Constant::HartreeEnergy_cm << "    iterations: " << iterations << std::endl;
    double full_gap = current_energy - old_energy;

    while(fabs(current_energy/energy - 1.) > core->GetEnergyTolerance())
    {
        double gap = (energy - current_energy)/full_gap;
        amount += gap;
        iterations = core->UpdateExcitedState(&sigma_s, sigma, amount);
        current_energy = sigma_s.Energy();
        *logstream << "    " << amount << "   " << current_energy * Constant::HartreeEnergy_cm
                   << "   iterations: " << iterations << std::endl;
    }
    SetSigmaAmount(info, amount);
    *outstream << "Final sigma amount: " << std::setprecision(10) << amount << " (gives "
               << std::setprecision(12) << current_energy * Constant::HartreeEnergy_cm << ")" << std::endl;
}

void ExcitedStates::SetSigmaAmount(const StateInfo& info, double amount)
{
    SecondOrderAmount[info] = amount;
}

double ExcitedStates::GetSigmaAmount(const StateInfo& info) const
{
    SigmaAmount::const_iterator it = SecondOrderAmount.find(info);
    if(it != SecondOrderAmount.end())
        return it->second;
    else
        return 0.;
}

void ExcitedStates::MultiplyByR(const DiscreteState* previous, DiscreteState* current) const
{
    current->ReSize(previous->Size());

    std::vector<double> Potential(core->GetHFPotential());
    std::vector<double> LocalExchange(core->GetLocalExchangeApproximation());

    unsigned int i;
    for(i=0; i<Potential.size(); i++)
        Potential[i] += LocalExchange[i];

    std::vector<double> dV(Potential.size());
    StateIntegrator I(lattice);
    I.GetDerivativeStart(Potential, dV, 0);
    I.GetDerivativeEnd(Potential, dV, Potential.size());
    I.GetDerivative(Potential, dV, 2, Potential.size()-2);

    const double* R = lattice->R();
    const double* dR = lattice->dR();
    double kappa = current->Kappa();

    double AlphaSquared = Constant::AlphaSquared; // Set to zero to remove effect of core potential

    for(i=0; i<previous->Size(); i++)
    {
        current->f[i] = previous->f[i] * R[i];
        current->g[i] = (R[i]*previous->df[i]/dR[i] + (1. + kappa)*previous->f[i])
                        /(2. + AlphaSquared*Potential[i]);
        current->df[i] = previous->df[i] * R[i] + previous->f[i] * dR[i];
        current->dg[i]
            = (kappa*previous->f[i]*dR[i]/R[i] + 2.*previous->df[i]
               + AlphaSquared*dV[i]*(previous->g[i]*R[i] - current->g[i]))
                  /(2. + AlphaSquared*Potential[i])
              + previous->dg[i] * R[i];
    }

    Orthogonalise(current);
    current->ReNormalise(lattice);
    current->SetEnergy(I.HamiltonianMatrixElement(*current, *current, *core));
}

void ExcitedStates::MultiplyBySinR(const DiscreteState* previous, DiscreteState* current) const
{
    current->ReSize(previous->Size());

    std::vector<double> Potential(core->GetHFPotential());
    std::vector<double> LocalExchange(core->GetLocalExchangeApproximation());

    unsigned int i;
    for(i=0; i<Potential.size(); i++)
        Potential[i] += LocalExchange[i];

    std::vector<double> dV(Potential.size());
    StateIntegrator I(lattice);
    I.GetDerivativeStart(Potential, dV, 0);
    I.GetDerivativeEnd(Potential, dV, Potential.size());
    I.GetDerivative(Potential, dV, 2, Potential.size()-2);

    const double* R = lattice->R();
    const double* dR = lattice->dR();
    double kappa = current->Kappa();
    double k = Constant::Pi/lattice->R(previous->Size());

    double AlphaSquared = Constant::AlphaSquared; // Set to zero to remove effect of core potential

    for(i=0; i<previous->Size(); i++)
    {
        double sinp = sin(k*R[i]);
        double kcosp = k * cos(k*R[i]);

        current->f[i] = previous->f[i] * sinp;
        current->g[i] = (sinp*previous->df[i]/dR[i] + (kcosp + kappa/R[i]*sinp)*previous->f[i])
                        /(2. + AlphaSquared*Potential[i]);
        current->df[i] = previous->df[i] * sinp + previous->f[i] * kcosp * dR[i];
        current->dg[i]
            = ((kappa/R[i]*kcosp - k*k*sinp) * previous->f[i] * dR[i]
                + 2.* kcosp * previous->df[i]
                + AlphaSquared*dV[i]*(previous->g[i]*sinp - current->g[i]))
                    /(2. + AlphaSquared*Potential[i])
              + previous->dg[i] * sinp;
    }

    Orthogonalise(current);
    current->ReNormalise(lattice);
    current->SetEnergy(I.HamiltonianMatrixElement(*current, *current, *core));
}

void ExcitedStates::MultiplyByRSinR(const DiscreteState* previous, DiscreteState* current) const
{
    current->ReSize(previous->Size());

    std::vector<double> Potential(core->GetHFPotential());
    std::vector<double> LocalExchange(core->GetLocalExchangeApproximation());

    unsigned int i;
    for(i=0; i<Potential.size(); i++)
        Potential[i] += LocalExchange[i];

    std::vector<double> dV(Potential.size());
    StateIntegrator I(lattice);
    I.GetDerivativeStart(Potential, dV, 0);
    I.GetDerivativeEnd(Potential, dV, Potential.size());
    I.GetDerivative(Potential, dV, 2, Potential.size()-2);

    const double* R = lattice->R();
    const double* dR = lattice->dR();
    double kappa_c = current->Kappa();
    double kappa_p = previous->Kappa();
    double k = Constant::Pi/lattice->R(previous->Size());

    double AlphaSquared = Constant::AlphaSquared; // Set to zero to remove effect of core potential

    for(i=0; i<previous->Size(); i++)
    {
        double sinp = sin(k*R[i]);
        double kcosp = k * cos(k*R[i]);

        current->f[i] = previous->f[i] * R[i] * sinp;
        current->g[i] = (R[i] * sinp * previous->df[i]/dR[i]
                        + (R[i] * kcosp + (1.+ kappa_c) * sinp)*previous->f[i])
                        /(2. + AlphaSquared*Potential[i]);
        current->df[i] = previous->df[i] * R[i] * sinp 
                        + previous->f[i] * (R[i] * kcosp + sinp) * dR[i];
        current->dg[i]
            = (((kappa_p/R[i] - R[i]*k*k) * sinp + (2.+ kappa_c) * kcosp)
                    * previous->f[i] * dR[i]
               + (2.*R[i]*kcosp + (2. + kappa_c - kappa_p) * sinp) * previous->df[i]
               + AlphaSquared*dV[i]*(previous->g[i]*R[i]*sinp - current->g[i]))
                    /(2. + AlphaSquared*Potential[i])
              + previous->dg[i] * R[i] * sinp;
    }

    Orthogonalise(current);
    current->ReNormalise(lattice);
    current->SetEnergy(I.HamiltonianMatrixElement(*current, *current, *core));
}

void ExcitedStates::Orthogonalise(DiscreteState* current) const
{
    const double* dR = lattice->dR();
    current->ReNormalise(lattice);

    // Orthogonalise to core
    ConstStateIterator it = core->GetConstStateIterator();
    while(!it.AtEnd())
    {
        const DiscreteState* other = it.GetState();
        if((other->Kappa() == current->Kappa()) && (other->RequiredPQN() != current->RequiredPQN()))
        {
            double S = 0.;
            unsigned int i;
            for(i=0; i<mmin(other->Size(), current->Size()); i++)
            {   S = S + ((other->f[i])*(current->f[i]) + Constant::AlphaSquared*(other->g[i])*(current->g[i])) * dR[i];
            }

            for(i=0; i<mmin(other->Size(), current->Size()); i++)
            {
                current->f[i] = current->f[i] - S * other->f[i];
                current->g[i] = current->g[i] - S * other->g[i];
                current->df[i] = current->df[i] - S * other->df[i];
                current->dg[i] = current->dg[i] - S * other->dg[i];
            }
            current->ReNormalise(lattice);
        }
        it.Next();
    }

    // Orthogonalise to other excited states.
    ConstStateIterator ex_it = GetConstStateIterator();
    while(!ex_it.AtEnd())
    {
        const DiscreteState* other = ex_it.GetState();
        if((other->Kappa() == current->Kappa()) && (other->RequiredPQN() != current->RequiredPQN())
           && !core->GetState(StateInfo(other)))
        {
            double S = 0.;
            unsigned int i;
            for(i=0; i<mmin(other->Size(), current->Size()); i++)
            {   S += ((other->f[i])*(current->f[i]) + Constant::AlphaSquared*(other->g[i])*(current->g[i]))
                    * dR[i];
            }

            for(i=0; i<mmin(other->Size(), current->Size()); i++)
            {
                current->f[i] = current->f[i] - S * other->f[i];
                current->g[i] = current->g[i] - S * other->g[i];
                current->df[i] = current->df[i] - S * other->df[i];
                current->dg[i] = current->dg[i] - S * other->dg[i];
            }
            current->ReNormalise(lattice);
        }
        ex_it.Next();
    }
    StateIntegrator I(lattice);
    current->SetEnergy(I.HamiltonianMatrixElement(*current, *current, *core));
}

double ExcitedStates::TestOrthogonalityIncludingCore() const
{
    double max_orth = 0.;

    ConstStateIterator it = GetConstStateIterator();
    ConstStateIterator jt = GetConstStateIterator();

    const State* max_i;
    const State* max_j;

    it.First();
    while(!it.AtEnd())
    {
        // Core
        jt = core->GetConstStateIterator();
        jt.First();
        while(!jt.AtEnd())
        {
            double orth = fabs(it.GetState()->Overlap(*jt.GetState(), lattice));
            if(orth > max_orth)
            {   max_orth = orth;
                max_i = it.GetState();
                max_j = jt.GetState();
            }

            jt.Next();
        }
        
        // Rest of excited states
        jt = it;
        jt.Next();
        while(!jt.AtEnd())
        {
            double orth = fabs(it.GetState()->Overlap(*jt.GetState(), lattice));
            if(orth > max_orth)
            {   max_orth = orth;
                max_i = it.GetState();
                max_j = jt.GetState();
            }

            jt.Next();
        }
        it.Next();
    }

    if(DebugOptions.OutputHFExcited())
    {   *outstream << "<" << max_i->Name() << " | " << max_j->Name() << "> = " << max_orth << std::endl;
    }
    return max_orth;
}

void ExcitedStates::Clear()
{
    StateManager::Clear();
    ClearSigmas();
}
