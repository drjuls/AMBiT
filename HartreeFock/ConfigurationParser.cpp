#include "ConfigurationParser.h"
#include "Include.h"
#include "Configuration/Configuration.h"
#include "NonRelInfo.h"

OccupationMap ConfigurationParser::ParseFractionalConfiguration(const std::string& configuration)
{
    OccupationMap ret;
    Configuration config(configuration);

    // Split non-relativistic configuration into relativistic states
    config.First();
    while(!config.AtEnd())
    {
        pOrbital s1, s2;

        NonRelInfo info(config.GetInfo());
        double occupancy = config.GetOccupancy();
        int L = info.L();

        if(L == 0)
        {   ret[info.GetFirstRelativisticInfo()] = occupancy;
        }
        else
        {   // Split electrons between subshells
            double first_fraction = double(L + 1)/double(2 * L + 1);
            ret[info.GetFirstRelativisticInfo()] = first_fraction * occupancy;
            ret[info.GetSecondRelativisticInfo()] = (1. - first_fraction) * occupancy;
        }

        config.Next();
    }

    return ret;
}