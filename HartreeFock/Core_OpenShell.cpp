#include "Include.h"
#include "Core.h"

void Core::Ionise(StateInfo removed_electron)
{
    DiscreteState* from_state = GetState(removed_electron);
    if(from_state != NULL)
    {
        if(OpenShellStates.find(removed_electron) == OpenShellStates.end())
        {   OpenShellStates[removed_electron] = from_state->Occupancy();
        }
        if(from_state->Occupancy() > 1.)
            from_state->SetOccupancy(from_state->Occupancy() - 1.);
        else
        {   OpenShellStorage[removed_electron] = from_state;
            AllStates.erase(AllStates.find(removed_electron));
        }
        
        Charge++;
        UpdateHFPotential();
    }
}

void Core::ToggleClosedShellCore()
{
    std::map<StateInfo, double>::iterator it = OpenShellStates.begin();
    while(it != OpenShellStates.end())
    {
        DiscreteState* ds = GetState(it->first);
        if(ds)
        {   StateSet::iterator kill_it = AllStates.find(it->first);
            if(kill_it != AllStates.end())
            {   Charge = Charge + ds->Occupancy();
                AllStates.erase(kill_it);
            }
            OpenShellStorage[it->first] = ds;
        }
        it++;
    }
    UpdateHFPotential();
}

void Core::ToggleOpenShellCore()
{
    std::map<StateInfo, double>::iterator it = OpenShellStates.begin();
    while(it != OpenShellStates.end())
    {
        double new_occupancy = it->second;
        DiscreteState* ds;

        if(AllStates.find(it->first) == AllStates.end())
        {   Charge = Charge - new_occupancy;
            ds = OpenShellStorage[it->first].GetState();
            ds->SetOccupancy(new_occupancy);
            AddState(ds);
        }
        else
        {   ds = AllStates.find(it->first)->second.GetState();
            if(ds->Occupancy() != new_occupancy)
            {   Charge = Charge + ds->Occupancy() - new_occupancy;
                ds->SetOccupancy(new_occupancy);
            }
        }
        it++;
    }
    UpdateHFPotential();
    OpenShellStorage.clear();
}

bool Core::IsOpenShellState(const StateInfo& info) const
{
    return (OpenShellStates.find(info) != OpenShellStates.end());
}
