#include "Include.h"
#include "StateManager.h"
#include "StateIterator.h"

StateManager::StateManager(Lattice* lat, unsigned int atomic_number, int ion_charge):
    lattice(lat), Z(atomic_number), Charge(ion_charge)
{}

StateManager::~StateManager(void)
{
    Clear();
}

const DiscreteState* StateManager::GetState(const StateInfo& info) const
{
    StateSet::const_iterator it;

    if((it = AllStates.find(info)) != AllStates.end())
        return it->second.GetState();
    else
        return NULL;
}

DiscreteState* StateManager::GetState(const StateInfo& info)
{
    StateSet::iterator it;

    if((it = AllStates.find(info)) != AllStates.end())
        return it->second.GetState();
    else
        return NULL;
}

void StateManager::AddState(DiscreteState* s)
{
    StateInfo info(s);
    StatePointer sp(s);

    AllStates.insert(StateSet::value_type(info, sp));
}

void StateManager::Clear()
{
    // Delete current states
    StateSet::iterator it = AllStates.begin();
    while(it != AllStates.end())
    {   it->second.DeleteState();
        it++;
    }
    AllStates.clear();
}

double StateManager::TestOrthogonality() const
{
    double max_orth = 0.;

    ConstStateIterator it = GetConstStateIterator();
    ConstStateIterator jt = GetConstStateIterator();

    it.First();
    while(!it.AtEnd())
    {
        jt = it;
        jt.Next();
        while(!jt.AtEnd())
        {
            double orth = fabs(it.GetState()->Overlap(*jt.GetState()));
            if(orth > max_orth)
                max_orth = orth;

            jt.Next();
        }
        it.Next();
    }

    return max_orth;
}

StateIterator StateManager::GetStateIterator()
{
    StateIterator it(this);
    return it;
}

ConstStateIterator StateManager::GetConstStateIterator() const
{
    ConstStateIterator it(this);
    return it;
}

void StateManager::Write(FILE* fp) const
{
    unsigned int num_states = NumStates();

    fwrite(&num_states, sizeof(unsigned int), 1, fp);
    ConstStateIterator it = GetConstStateIterator();
    while(!it.AtEnd())
    {
        it.GetState()->Write(fp);
        it.Next();
    }
}

void StateManager::Read(FILE* fp)
{
    unsigned int num_core, i;
    DiscreteState ds(lattice);

    // Read states
    fread(&num_core, sizeof(unsigned int), 1, fp);
    for(i = 0; i<num_core; i++)
    {
        ds.Read(fp);
        DiscreteState* current = GetState(StateInfo(&ds));
        if(current)
        {   *current = ds;
        }
    }
}
