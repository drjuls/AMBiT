#include "DiscreteExcitedStates.h"
#include "Include.h"

#include "HartreeFock/StateIntegrator.h"

DiscreteState* DiscreteExcitedStates::GetState(const StateInfo& info)
{
    StateSet::iterator it;

    if((it = AllStates.find(info)) != AllStates.end())
        return dynamic_cast<DiscreteState*>(it->second.GetState());
    else
        return NULL;
}

const DiscreteState* DiscreteExcitedStates::GetState(const StateInfo& info) const
{
    StateSet::const_iterator it;

    if((it = AllStates.find(info)) != AllStates.end())
        return dynamic_cast<const DiscreteState*>(it->second.GetState());
    else
        return NULL;
}

DiscreteStateIterator DiscreteExcitedStates::GetDiscreteStateIterator()
{
    DiscreteStateIterator it(this);
    return it;
}

ConstDiscreteStateIterator DiscreteExcitedStates::GetConstDiscreteStateIterator() const
{
    ConstDiscreteStateIterator it(this);
    return it;
}

void DiscreteExcitedStates::Read(FILE* fp)
{
    Clear();
    
    unsigned int num_core, i;

    // Read states
    fread(&num_core, sizeof(unsigned int), 1, fp);
    for(i = 0; i<num_core; i++)
    {
        DiscreteState* ds = new DiscreteState(lattice);
        ds->Read(fp);
        AddState(ds);
    }
}

void DiscreteExcitedStates::Write(FILE* fp) const
{
    unsigned int num_states = NumStates();

    fwrite(&num_states, sizeof(unsigned int), 1, fp);
    ConstDiscreteStateIterator it = GetConstDiscreteStateIterator();
    while(!it.AtEnd())
    {
        it.GetState()->Write(fp);
        it.Next();
    }
}

void DiscreteExcitedStates::MultiplyByR(const DiscreteState* previous, DiscreteState* current) const
{
    current->ReSize(previous->Size());

    std::vector<double> Potential(core->GetHFPotential());
    std::vector<double> LocalExchange(core->GetLocalExchangeApproximation());

    unsigned int i;
    for(i=0; i<Potential.size(); i++)
        Potential[i] += LocalExchange[i];

    std::vector<double> dV(Potential.size());
    StateIntegrator I(*lattice);
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
    current->ReNormalise();
    current->SetEnergy(I.HamiltonianMatrixElement(*current, *current, *core));
}

void DiscreteExcitedStates::MultiplyBySinR(const DiscreteState* previous, DiscreteState* current) const
{
    current->ReSize(previous->Size());

    std::vector<double> Potential(core->GetHFPotential());
    std::vector<double> LocalExchange(core->GetLocalExchangeApproximation());

    unsigned int i;
    for(i=0; i<Potential.size(); i++)
        Potential[i] += LocalExchange[i];

    std::vector<double> dV(Potential.size());
    StateIntegrator I(*lattice);
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
    current->ReNormalise();
    current->SetEnergy(I.HamiltonianMatrixElement(*current, *current, *core));
}

void DiscreteExcitedStates::MultiplyByRSinR(const DiscreteState* previous, DiscreteState* current) const
{
    current->ReSize(previous->Size());

    std::vector<double> Potential(core->GetHFPotential());
    std::vector<double> LocalExchange(core->GetLocalExchangeApproximation());

    unsigned int i;
    for(i=0; i<Potential.size(); i++)
        Potential[i] += LocalExchange[i];

    std::vector<double> dV(Potential.size());
    StateIntegrator I(*lattice);
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
    current->ReNormalise();
    current->SetEnergy(I.HamiltonianMatrixElement(*current, *current, *core));
}

void DiscreteExcitedStates::Orthogonalise(DiscreteState* current) const
{
    const double* dR = lattice->dR();
    current->ReNormalise();

    // Orthogonalise to core
    ConstDiscreteStateIterator it = core->GetConstDiscreteStateIterator();
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
//            std::cout << "  orth: " << current->Name() << " " << other->Name() << " \t" << S << std::endl;

            for(i=0; i<mmin(other->Size(), current->Size()); i++)
            {
                current->f[i] = current->f[i] - S * other->f[i];
                current->g[i] = current->g[i] - S * other->g[i];
                current->df[i] = current->df[i] - S * other->df[i];
                current->dg[i] = current->dg[i] - S * other->dg[i];
            }
            current->ReNormalise();
        }
        it.Next();
    }

    // Orthogonalise to other excited states
    ConstDiscreteStateIterator ex_it = GetConstDiscreteStateIterator();
    while(!ex_it.AtEnd())
    {
        const DiscreteState* other = dynamic_cast<const DiscreteState*>(ex_it.GetState());
        if((other->Kappa() == current->Kappa()) && (other->RequiredPQN() != current->RequiredPQN())
           && !core->GetState(StateInfo(other)))
        {
            double S = 0.;
            unsigned int i;
            for(i=0; i<mmin(other->Size(), current->Size()); i++)
            {   S += ((other->f[i])*(current->f[i]) + Constant::AlphaSquared*(other->g[i])*(current->g[i]))
                    * dR[i];
            }
//            std::cout << "  orth: " << current->Name() << " " << other->Name() << " \t" << S << std::endl;

            for(i=0; i<mmin(other->Size(), current->Size()); i++)
            {
                current->f[i] = current->f[i] - S * other->f[i];
                current->g[i] = current->g[i] - S * other->g[i];
                current->df[i] = current->df[i] - S * other->df[i];
                current->dg[i] = current->dg[i] - S * other->dg[i];
            }
            current->ReNormalise();
        }
        ex_it.Next();
    }
}
