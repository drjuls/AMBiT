#include "Include.h"
#include "OrbitalInfo.h"
#include "Universal/MathConstant.h"

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

    if(kappa < -1)
        ret.append(1, '+');

    return ret;
}
