#include "Include.h"
#include "StateInfo.h"
#include "Universal/Constant.h"

bool StateInfo::operator<(const StateInfo& other) const
{
    if(this->pqn < other.pqn)
       return true;
    else if(this->pqn > other.pqn)
        return false;

    if(this->l < other.l)
        return true;
    else if(this->l > other.l)
        return false;
    else return (this->kappa > other.kappa);
}

bool StateInfo::operator==(const StateInfo& other) const
{
    if((this->pqn == other.pqn) &&
       (this->kappa == other.kappa))
        return true;
    else
        return false;
}

bool StateInfo::operator!=(const StateInfo& other) const
{
    return !(*this == other);
}

std::string StateInfo::Name() const
{
    char buffer[20];
    sprintf(buffer, "%d", pqn);
    std::string ret(buffer);

    ret.append(1, Constant::SpectroscopicNotation[L()]);
    if(kappa < -1)
        ret.append(1, '+');

    return ret;
}
