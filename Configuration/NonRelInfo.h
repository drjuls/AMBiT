#ifndef NON_REL_INFO_H
#define NON_REL_INFO_H

#include "SingleParticleInfo.h"
#include "RelativisticInfo.h"
#include <string>
#include "Universal/Constant.h"

class NonRelInfo : public SingleParticleInfo
{
    /** Stores a non-relativistic single particle state info. */
public:
    NonRelInfo(unsigned int principal_qn, unsigned int l):
        SingleParticleInfo(principal_qn, -(int(l) + 1))
    {}
    NonRelInfo(const SingleParticleInfo& other):
        SingleParticleInfo(other.PQN(), -(int(other.L()) + 1))
    {}
    virtual ~NonRelInfo(void) {}

    /** Get relativistic state info corresponding to kappa = -(L+1). */
    RelativisticInfo GetFirstRelativisticInfo() const;
    /** Get relativistic state info corresponding to kappa = L. */
    RelativisticInfo GetSecondRelativisticInfo() const;

    inline std::string Name() const;
};

inline RelativisticInfo NonRelInfo::GetFirstRelativisticInfo() const
{
    return RelativisticInfo(pqn, kappa);
}

inline RelativisticInfo NonRelInfo::GetSecondRelativisticInfo() const
{
    if(l == 0)
        return RelativisticInfo(pqn, -1);
    else
        return RelativisticInfo(pqn, l);
}

inline std::string NonRelInfo::Name() const
{
    char buffer[20];
    sprintf(buffer, "%d", pqn);
    std::string ret(buffer);

    ret.append(1, Constant::SpectroscopicNotation[L()]);
    return ret;
}

#endif
