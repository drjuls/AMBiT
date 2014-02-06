#include "Include.h"
#include "OrbitalInfo.h"
#include "Universal/MathConstant.h"

bool OrbitalInfo::operator<(const OrbitalInfo& other) const
{
    // Sort on pqn
    if(this->pqn < other.pqn)
       return true;
    else if(this->pqn > other.pqn)
        return false;

    // Sort on abs(kappa):
    //  |-1| < |1| < |-2| < |2| < |-3| ...
    if(abs(this->kappa) < abs(other.kappa))
        return true;
    else if(abs(this->kappa) > abs(other.kappa))
        return false;
    // And finally sort on kappa itself
    else return (this->kappa < other.kappa);
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
