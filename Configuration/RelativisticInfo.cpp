#include "Include.h"
#include "RelativisticInfo.h"
#include "NonRelInfo.h"
#include "Universal/Constant.h"

NonRelInfo RelativisticInfo::GetNonRelInfo() const
{
    return NonRelInfo(pqn, l);
}

StateInfo RelativisticInfo::GetStateInfo() const
{
    return StateInfo(pqn, kappa);
}

std::string RelativisticInfo::Name() const
{
    char buffer[20];
    sprintf(buffer, "%d", pqn);
    std::string ret(buffer);

    ret.append(1, Constant::SpectroscopicNotation[L()]);
    if(kappa < -1)
        ret.append(1, '+');

    return ret;
}
