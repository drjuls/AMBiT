#include "Include.h"
#include "ElectronInfo.h"

bool ElectronInfo::operator<(const ElectronInfo& other) const
{
    // We redo the OrbitalInfo rather than calling OrbitalInfo::operator<()
    // to save a bit of time.

    // Sort on abs(kappa):
    //  |-1| < |1| < |-2| < |2| < |-3| ...
    if(abs(this->kappa) < abs(other.kappa))
        return true;
    else if(abs(this->kappa) > abs(other.kappa))
        return false;
    // Sort on kappa itself
    else if(this->kappa < other.kappa)
        return true;
    else if(this->kappa > other.kappa)
        return false;

    // Sort on pqn
    if(this->pqn < other.pqn)
        return true;
    else if(this->pqn > other.pqn)
        return false;

    // And sort on m such that AngularData generates projections sorted correctly
    // (so they don't need to be resorted).
    // Sorted Projections are required for GetProjectionDifferences()
    return (this->two_m > other.two_m);
}

bool ElectronInfo::operator==(const ElectronInfo& other) const
{
    if((this->pqn == other.pqn) &&
       (this->kappa == other.kappa) &&
       (this->two_m == other.two_m))
        return true;
    else
        return false;
}

bool ElectronInfo::operator !=(const ElectronInfo& other) const
{
    return !(*this == other);
}

std::string ElectronInfo::Name() const
{
    std::string name;
    if(is_hole)
        name.push_back('-');
    name.append(OrbitalInfo::Name());

    char buffer[20];
    sprintf(buffer, "(%d)", two_m);
    name.append(buffer);

    return name;
}
