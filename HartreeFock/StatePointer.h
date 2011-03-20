#ifndef STATE_POINTER_H
#define STATE_POINTER_H

#include "SingleParticleWavefunction.h"
#include "OrbitalInfo.h"

class StatePointer
{
public:
    StatePointer(Orbital* s = NULL): p(s) {}

    SingleParticleWavefunction& operator*() { return *p; }
    const SingleParticleWavefunction& operator*() const { return *p; }

    SingleParticleWavefunction* operator->() { return p; }
    const SingleParticleWavefunction* operator->() const { return p; }

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
    OrbitalInfo first_info(this->p);
    OrbitalInfo second_info(other.p);

    return (first_info < second_info);
}

#endif
