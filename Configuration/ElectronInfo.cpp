#include "Include.h"
#include "ElectronInfo.h"

bool ElectronInfo::operator<(const ElectronInfo& other) const
{
    if(this->pqn < other.pqn)
       return true;
    else if(this->pqn > other.pqn)
        return false;

    if(this->l < other.l)
        return true;
    else if(this->l > other.l)
        return false;

    if(this->kappa < other.kappa)
        return true;
    else if(this->kappa > other.kappa)
        return false;
    
    return (this->two_m < other.two_m);
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
    std::string name = StateInfo::Name();

    char buffer[20];
    sprintf(buffer, "(%d)", two_m);
    name.append(buffer);

    return name;
}
