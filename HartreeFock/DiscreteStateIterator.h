#ifndef DISCRETE_STATE_ITERATOR_H
#define DISCRETE_STATE_ITERATOR_H

#include "StateIterator.h"

class DiscreteStateIterator : public StateIterator
{
public:
    DiscreteStateIterator(StateManager* state_manager): StateIterator(state_manager) {}
    virtual ~DiscreteStateIterator(void) {}

    virtual DiscreteState* GetState();
};

class ConstDiscreteStateIterator : public ConstStateIterator
{
public:
    ConstDiscreteStateIterator(const StateManager* state_manager): ConstStateIterator(state_manager) {}
    virtual ~ConstDiscreteStateIterator(void) {}

    virtual const DiscreteState* GetState();
};

inline DiscreteState* DiscreteStateIterator::GetState()
{
    return dynamic_cast<DiscreteState*>(StateIterator::GetState());
}

inline const DiscreteState* ConstDiscreteStateIterator::GetState()
{
    return dynamic_cast<const DiscreteState*>(ConstStateIterator::GetState());
}

#endif
