#ifndef NON_REL_INFO_H
#define NON_REL_INFO_H

#include "StateInfo.h"

class NonRelInfo : public StateInfo
{
    /** Stores a non-relativistic single particle state info. */
public:
    NonRelInfo(unsigned int principal_qn, unsigned int l):
        StateInfo(principal_qn, -(int(l) + 1))
    {}
    NonRelInfo(const StateInfo& other):
        StateInfo(other.PQN(), -(int(other.L()) + 1))
    {}
    virtual ~NonRelInfo(void) {}

    /** Get relativistic state info corresponding to kappa = -(L+1). */
    StateInfo GetFirstRelativisticInfo() const;
    /** Get relativistic state info corresponding to kappa = L. */
    StateInfo GetSecondRelativisticInfo() const;

    virtual std::string Name() const;
};

inline StateInfo NonRelInfo::GetFirstRelativisticInfo() const
{
    return StateInfo(pqn, kappa);
}

inline StateInfo NonRelInfo::GetSecondRelativisticInfo() const
{
    if(l == 0)
        return StateInfo(pqn, -1);
    else
        return StateInfo(pqn, l);
}

#endif
