#include "StateIterator.h"

StateIterator::StateIterator(StateManager* state_manager):
    manager(state_manager)
{
    First();
}

void StateIterator::First()
{
    it = manager->AllStates.begin();
}

void StateIterator::Next()
{
    if(!AtEnd())
        it++;
}

bool StateIterator::AtEnd()
{
    return (it == manager->AllStates.end());
}

Orbital* StateIterator::GetState()
{
    return it->second.GetState();
}

StateInfo StateIterator::GetStateInfo()
{
    return it->first;
}

void StateIterator::ReplaceState(Orbital* s)
{
    StateInfo new_state(s);
    if(it->first != new_state)
        return;

    it->second.DeleteState();
    it->second.SetState(s);
}

StateIterator& StateIterator::operator=(const StateIterator& other)
{
    manager = other.manager;
    it = other.it;

    return *this;
}

ConstStateIterator::ConstStateIterator(const StateManager* state_manager):
    manager(state_manager)
{
    First();
}

void ConstStateIterator::First()
{
    it = manager->AllStates.begin();
}

void ConstStateIterator::Next()
{
    if(!AtEnd())
        it++;
}

bool ConstStateIterator::AtEnd()
{
    return (it == manager->AllStates.end());
}

const Orbital* ConstStateIterator::GetState()
{
    return it->second.GetState();
}

StateInfo ConstStateIterator::GetStateInfo()
{
    return it->first;
}

ConstStateIterator& ConstStateIterator::operator=(const ConstStateIterator& other)
{
    manager = other.manager;
    it = other.it;

    return *this;
}

bool ConstStateIterator::operator==(const ConstStateIterator& other) const
{
    if(manager == other.manager)
        if(it == other.it)
            return true;

    return false;
}

bool ConstStateIterator::operator!=(const ConstStateIterator& other) const
{
    if(manager == other.manager)
        if(it == other.it)
            return false;

    return true;
}
