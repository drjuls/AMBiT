#ifndef STATE_INFO_H
#define STATE_INFO_H

#include "Orbital.h"
#include <string>

class StateInfo
{
    /** Single particle state information.
        Stores pqn and kappa (L). Has an inbuilt ordering.
     */
public:
    StateInfo(unsigned int principal_qn, int kappa);
    StateInfo(const Orbital* s);
    StateInfo(const StateInfo& other);
    virtual ~StateInfo(void) {}

    virtual bool operator<(const StateInfo& other) const;
    virtual bool operator==(const StateInfo& other) const;
    virtual bool operator!=(const StateInfo& other) const;

    inline unsigned int PQN() const { return pqn; }
    inline int Kappa() const { return kappa; }
    inline unsigned int L() const { return l; }
    inline double J() const;
    inline unsigned int TwoJ() const;

    inline unsigned int MaxNumElectrons() const { return 2*abs(kappa); }
    virtual std::string Name() const;

protected:
    unsigned int pqn;
    int kappa;
    unsigned int l;
};

inline StateInfo::StateInfo(const Orbital* s)
{
    pqn = s->RequiredPQN();
    kappa = s->Kappa();
    if(kappa > 0)
        l = (unsigned int)kappa;
    else
        l = (unsigned int)(-kappa-1);
}

inline StateInfo::StateInfo(unsigned int principal_qn, int kap):
    pqn(principal_qn), kappa(kap)
{
    if(kappa > 0)
        l = (unsigned int)kappa;
    else
        l = (unsigned int)(-kappa-1);
}

inline StateInfo::StateInfo(const StateInfo& other):
    pqn(other.pqn), kappa(other.kappa), l(other.l)
{}

inline double StateInfo::J() const
{
    return fabs(double(kappa)) - 0.5;
}

inline unsigned int StateInfo::TwoJ() const
{
    return (2*abs(kappa) - 1);
}

#endif
