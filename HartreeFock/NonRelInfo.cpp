#include "Include.h"
#include "NonRelInfo.h"
#include "Universal/Constant.h"

std::string NonRelInfo::Name() const
{
    char buffer[20];
    sprintf(buffer, "%d", pqn);
    std::string ret(buffer);

    ret.append(1, Constant::SpectroscopicNotation[L()]);
    return ret;
}
