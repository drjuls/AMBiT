#include "Include.h"
#include "ElectronInfo.h"

const ElectronInfo& ElectronInfo::operator=(const ElectronInfo& other)
{
    OrbitalInfo::operator=(other);
    two_m = other.two_m;
    is_hole = other.is_hole;

    return *this;
}

bool ElectronInfo::operator<(const ElectronInfo& other) const
{
    // We redo the OrbitalInfo rather than calling OrbitalInfo::operator<()
    // to save a bit of time.

    // Sort on abs(kappa):
    //  |-1| < |1| < |-2| < |2| < |-3| ...
    if(abs(this->kappa) != abs(other.kappa))
        return (abs(this->kappa) < abs(other.kappa));
    else if(this->kappa != other.kappa)
        return (this->kappa < other.kappa);
    else if(this->pqn != other.pqn)
        return (this->pqn < other.pqn);

    // And sort on m such that AngularData generates projections sorted correctly
    // (so they don't need to be resorted).
    // Sorted Projections are required for GetProjectionDifferences()
    else if(is_hole)
        return (this->two_m < other.two_m);
    else
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
