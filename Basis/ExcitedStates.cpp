#include "Include.h"
#include "ExcitedStates.h"
#include "HartreeFock/Core.h"
#include "Universal/Constant.h"
#include "MBPT/MBPTCalculator.h"

ExcitedStates::ExcitedStates(Lattice* lattice, const Core* atom_core):
    StateManager(lattice, (unsigned int)atom_core->GetZ(), (unsigned int)atom_core->GetCharge()),
    core(atom_core)
{}

ExcitedStates::~ExcitedStates(void)
{
    Clear();

    SigmaMap::iterator sigma = SecondOrderSigma.begin();
    while(sigma != SecondOrderSigma.end())
    {   delete sigma->second;
        sigma++;
    }
}

void ExcitedStates::AddState(State* s)
{
    StateManager::AddState(s);
}

DiscreteState ExcitedStates::GetStateWithSigma(const StateInfo& info) const
{
    const DiscreteState* s = dynamic_cast<const DiscreteState*>(GetState(info));

    if(s != NULL)
    {   DiscreteState ds(*s);

        SigmaMap::const_iterator it = SecondOrderSigma.find(ds.Kappa());
        if((it != SecondOrderSigma.end()) && GetSigmaAmount(info))
            core->UpdateExcitedState(&ds, it->second, GetSigmaAmount(info));

        return ds;
    }
    else
        return DiscreteState(lattice);
}

double ExcitedStates::GetSecondOrderSigma(const StateInfo& info)
{
    State* s = GetState(info);
    if(s == NULL)
    {   s = new DiscreteState(lattice, info.RequiredPQN(), info.Kappa());
        core->CalculateExcitedState(s);
    }

    SigmaPotential* sigma;

    std::string sigma_file = *identifier + "." + s->Name() + ".sigma";
    sigma = new SigmaPotential(lattice, sigma_file, s->Size(), 100);
    sigma->Reset();

    SigmaMap::iterator it = SecondOrderSigma.find(s->Kappa());
    if(it != SecondOrderSigma.end())
    {   delete SecondOrderSigma[s->Kappa()];
        SecondOrderSigma.erase(it);
    }

    MBPTCalculator mbpt(lattice, core, this);
    double energy = mbpt.GetSecondOrderSigma(s, sigma);

    SecondOrderSigma[s->Kappa()] = sigma;
    sigma->Store();
    SetSigmaAmount(info, 1.);

    return energy;
}

bool ExcitedStates::RetrieveSecondOrderSigma(const StateInfo& info)
{
    State* s = GetState(info);
    if(s == NULL)
    {   s = new DiscreteState(lattice, info.RequiredPQN(), info.Kappa());
        core->CalculateExcitedState(s);
    }

    std::string sigma_file = *identifier + "." + s->Name() + ".sigma";
    SigmaPotential* sigma = new SigmaPotential(lattice, sigma_file);

    if(sigma->Size())
    {
        SigmaMap::iterator it = SecondOrderSigma.find(s->Kappa());
        if(it != SecondOrderSigma.end())
        {   delete SecondOrderSigma[s->Kappa()];
            SecondOrderSigma.erase(it);
        }

        SecondOrderSigma[s->Kappa()] = sigma;
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
    if(!info.Discrete())
        return;

    State* s = GetState(info);
    if(s == NULL)
    {   s = new DiscreteState(lattice, info.RequiredPQN(), info.Kappa());
        core->CalculateExcitedState(s);
    }

    DiscreteState* ds = dynamic_cast<DiscreteState*>(s);

    SigmaPotential* sigma;
    if(SecondOrderSigma.find(s->Kappa()) == SecondOrderSigma.end())
        GetSecondOrderSigma(info);
    sigma = SecondOrderSigma[s->Kappa()];

    double old_energy = s->Energy();
    double amount = 1.0;
    DiscreteState sigma_s(*ds);
    unsigned int iterations = core->UpdateExcitedState(&sigma_s);
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
