#ifndef STATE_INFO_H
#define STATE_INFO_H

#include "State.h"
#include "DiscreteState.h"

class StateInfo
{
    /** StateInfo is an ordered indexing class.
        It contains all the information needed to identify a state.
        There is a zero (null) StateInfo.
     */
public:
    inline StateInfo();
    inline StateInfo(int null_value) { StateInfo(); };
    inline StateInfo(const State* s);
    inline StateInfo(double principal_qn, int kappa, bool discrete = true);
    ~StateInfo() {}

    inline bool operator<(const StateInfo& other) const;
    inline bool operator==(const StateInfo& other) const;
    inline bool operator!=(const StateInfo& other) const;

    inline bool Discrete() const { return discrete; }
    inline double RequiredPQN() const { return requiredPQN; }
    inline int Kappa() const { return kappa; }

    inline bool IsNull() const { return (kappa == 0); }

private:
    bool discrete;
    double requiredPQN;
    int kappa;
    unsigned int l;
};

inline StateInfo::StateInfo()
{
    discrete = true;    // Good for ordering
    requiredPQN = 0.;
    kappa = 0;          // This is a good marker for a null state
    l = 0;
}

inline StateInfo::StateInfo(const State* s)
{
    const DiscreteState* ds = dynamic_cast<const DiscreteState*>(s);
    if(ds != NULL)
    {
        discrete = true;
        requiredPQN = ds->RequiredPQN();
        kappa = ds->Kappa();
        l = ds->L();
    }
    else
    {
        discrete = false;
        requiredPQN = s->Nu();
        kappa = s->Kappa();
        l = s->L();
    }
}

inline StateInfo::StateInfo(double principal_qn, int kap, bool is_discrete):
    discrete(is_discrete), requiredPQN(principal_qn), kappa(kap)
{
    if(kappa > 0)
        l = (unsigned int)kappa;
    else
        l = (unsigned int)(-kappa-1);
}

inline bool StateInfo::operator<(const StateInfo& other) const
{
    // Discrete states < continuum states
    if(this->discrete && !other.discrete)
        return true;
    else if(!this->discrete && other.discrete)
        return false;
    
    if(this->requiredPQN < other.requiredPQN)
       return true;
    else if(this->requiredPQN > other.requiredPQN)
        return false;

    if(this->l < other.l)
        return true;
    else if(this->l > other.l)
        return false;
    else return (this->kappa > other.kappa);
}

inline bool StateInfo::operator==(const StateInfo& other) const
{
    if((this->discrete == other.discrete) &&
       (this->requiredPQN == other.requiredPQN) &&
       (this->kappa == other.kappa))
        return true;
    else
        return false;
}

inline bool StateInfo::operator!=(const StateInfo& other) const
{
    return !(*this == other);
}

#endif