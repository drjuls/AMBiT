#ifndef STATE_POINTER_H
#define STATE_POINTER_H

#include "State.h"
#include "StateInfo.h"

class StatePointer
{
public:
    StatePointer(State* s = NULL): p(s) {}

    State& operator*() { return *p; }
    const State& operator*() const { return *p; }

    State* operator->() { return p; }
    const State* operator->() const { return p; }

    bool IsNull() const { return p == NULL; }
    State* GetState() { return p; }
    const State* GetState() const { return p; }
    void SetState(State* s) { p = s; }

    void DeleteState()
    {   delete p;
        p = NULL;
    }

    inline bool operator<(const StatePointer& other) const;
private:
    State* p;
};

inline bool StatePointer::operator<(const StatePointer& other) const
{
    StateInfo first_info(this->p);
    StateInfo second_info(other.p);

    return (first_info < second_info);
}

#endif