#include "Include.h"
#include "OrbitalInfo.h"
#include "Universal/MathConstant.h"

const OrbitalInfo& OrbitalInfo::operator=(const OrbitalInfo& other)
{
    pqn = other.pqn;
    kappa = other.kappa;

    return *this;
}

bool OrbitalInfo::operator<(const OrbitalInfo& other) const
{
    // Sort on abs(kappa) first:
    //  |-1| < |1| < |-2| < |2| < |-3| ...
    if(abs(this->kappa) != abs(other.kappa))
        return (abs(this->kappa) < abs(other.kappa));
    else if(this->kappa != other.kappa)
        return (this->kappa < other.kappa);
    else
        return (this->pqn < other.pqn);
}

bool OrbitalInfo::operator==(const OrbitalInfo& other) const
{
    if((this->pqn == other.pqn) &&
       (this->kappa == other.kappa))
        return true;
    else
        return false;
}

bool OrbitalInfo::operator!=(const OrbitalInfo& other) const
{
    return !(*this == other);
}

std::string OrbitalInfo::Name() const
{
    char buffer[20];
    sprintf(buffer, "%d", pqn);
    std::string ret(buffer);

    ret.append(1, MathConstant::Instance()->GetSpectroscopicNotation(L()));

#ifdef USE_ALT_STATE_NOTATION
    if(kappa > 0)
        ret.append(1, '-');
    else
        ret.append(1, ' ');
#else
    if(kappa < -1)
        ret.append(1, '+');
#endif

    return ret;
}
