#ifndef NON_REL_INFO_H
#define NON_REL_INFO_H

#include "OrbitalInfo.h"
#include "Universal/MathConstant.h"

namespace Ambit
{
class NonRelInfo : public OrbitalInfo
{
    /** Stores a non-relativistic single particle state info. */
public:
    NonRelInfo(int principal_qn, int kappa_or_l):
        OrbitalInfo(principal_qn, kappa_or_l)
    {   if(kappa_or_l >= 0)
            kappa = -(kappa_or_l + 1);
    }
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

    virtual int MaxNumElectrons() const override { return -4 * kappa - 2; } //!< Max occupancy of shell = 2(2L+1)
    virtual std::string Name() const override;
};

inline OrbitalInfo NonRelInfo::GetFirstRelativisticInfo() const
{
    return OrbitalInfo(pqn, kappa);
}

inline OrbitalInfo NonRelInfo::GetSecondRelativisticInfo() const
{
    if(kappa == -1)
        return OrbitalInfo(pqn, kappa);
    else
        return OrbitalInfo(pqn, -kappa-1);
}

inline std::string NonRelInfo::Name() const
{
    char buffer[20];
    sprintf(buffer, "%d", pqn);
    std::string ret(buffer);

    ret.append(1, MathConstant::Instance()->GetSpectroscopicNotation(L()));
    return ret;
}

}
#endif
