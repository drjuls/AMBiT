#ifndef NON_REL_INFO_H
#define NON_REL_INFO_H

#include <set>
#include "OrbitalInfo.h"

class NonRelInfo : public OrbitalInfo
{
    /** Stores a non-relativistic single particle state info. */
public:
    NonRelInfo(unsigned int principal_qn, unsigned int ll):
        OrbitalInfo(principal_qn, -(int(ll) + 1))
    {}
    NonRelInfo(const NonRelInfo& other):
        OrbitalInfo(other.pqn, other.kappa)
    {}
    NonRelInfo(const OrbitalInfo& other):
        OrbitalInfo(other.PQN(), -(int(other.L()) + 1))
    {}
    virtual ~NonRelInfo(void) {}

    /** Get relativistic state info corresponding to kappa = -(L+1). */
    OrbitalInfo GetFirstRelativisticInfo() const;
    /** Get relativistic state info corresponding to kappa = L. */
    OrbitalInfo GetSecondRelativisticInfo() const;

    virtual std::string Name() const;
};

inline OrbitalInfo NonRelInfo::GetFirstRelativisticInfo() const
{
    return OrbitalInfo(pqn, kappa);
}

inline OrbitalInfo NonRelInfo::GetSecondRelativisticInfo() const
{
    if(l == 0)
        return OrbitalInfo(pqn, -1);
    else
        return OrbitalInfo(pqn, l);
}


class NonRelInfoSet : public std::set<NonRelInfo>
{
public:
    void AddConfigs(const char* basis_def);
    void EraseConfigs(const char* basis_def);
};

#endif
