#include "SubsetBasis.h"
#include "Include.h"

SubsetBasis::SubsetBasis(pLattice lattice, pExcitedStates superset_basis):
    ExcitedStates(lattice, superset_basis->GetCore())
{
    SetSuperset(superset_basis);
}

SubsetBasis::~SubsetBasis()
{
    // After this, the ExcitedState destructor will be called. It in turn calls
    // StateManager::Clear(), which deletes all states. We do not own these states,
    // therefore we delete the StateSet here.
    AllStates.clear();
}

void SubsetBasis::CreateExcitedStates(const std::vector<unsigned int>& num_states_per_l)
{
    // Need to treat all kappas
    std::vector<unsigned int> num_states_so_far(2 * num_states_per_l.size() - 1, 0);
    
    StateIterator it = superset->GetStateIterator();
    it.First();
    while(!it.AtEnd())
    {
        OrbitalInfo si(it.GetOrbitalInfo());
        if(si.L() < num_states_per_l.size())
        {
            int kappa = si.Kappa();

            unsigned int index;
            if(kappa > 0)
                index = 2*kappa - 1;
            else
                index = 2*abs(kappa) - 2;

            if(num_states_so_far[index] < num_states_per_l[si.L()])
            {
                AddState(it.GetState());
                num_states_so_far[index]++;
            }
        }

        it.Next();
    }
}

void SubsetBasis::Clear()
{
    ClearSigmas();
    AllStates.clear();
}

pExcitedStatesConst SubsetBasis::GetSuperset() const
{
    return superset;
}

void SubsetBasis::SetSuperset(pExcitedStates superset_basis)
{
    superset = superset_basis;
}
