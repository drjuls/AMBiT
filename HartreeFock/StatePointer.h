#ifndef STATE_POINTER_H
#define STATE_POINTER_H

#include "State.h"
#include "StateInfo.h"

class StatePointer
{
public:
    StatePointer(Orbital* s = NULL): p(s) {}

    State& operator*() { return *p; }
    const State& operator*() const { return *p; }

    State* operator->() { return p; }
    const State* operator->() const { return p; }

    bool IsNull() const { return p == NULL; }
    Orbital* GetState() { return p; }
    const Orbital* GetState() const { return p; }
    void SetState(Orbital* s) { p = s; }

    void DeleteState()
    {   delete p;
        p = NULL;
    }

    inline bool operator<(const StatePointer& other) const;
private:
    Orbital* p;
};

inline bool StatePointer::operator<(const StatePointer& other) const
{
    StateInfo first_info(this->p);
    StateInfo second_info(other.p);

    return (first_info < second_info);
}

#endif
