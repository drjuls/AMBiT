#ifndef STATE_ITERATOR_H
#define STATE_ITERATOR_H

#include "StateManager.h"

class StateIterator
{
public:
    StateIterator(StateManager* state_manager);
    virtual ~StateIterator() {}

    virtual void First();
    virtual void Next();
    virtual bool AtEnd();

    virtual Orbital* GetState();
    virtual StateInfo GetStateInfo();
    virtual void ReplaceState(Orbital* s);

    virtual double Weight() { return 1.; }
    virtual StateIterator& operator=(const StateIterator& other);

protected:
    StateManager* manager;
    StateSet::iterator it;
};

class ConstStateIterator
{
public:
    ConstStateIterator(const StateManager* state_manager);
    virtual ~ConstStateIterator() {}

    virtual void First();
    virtual void Next();
    virtual bool AtEnd();

    virtual const Orbital* GetState();
    virtual StateInfo GetStateInfo();

    virtual double Weight() { return 1.; }
    virtual ConstStateIterator& operator=(const ConstStateIterator& other);

    virtual bool operator==(const ConstStateIterator& other) const;
    virtual bool operator!=(const ConstStateIterator& other) const;

protected:
    const StateManager* manager;
    StateSet::const_iterator it;
};

#endif
