#include "ConfigurationParser.h"
#include "Include.h"
#include "Configuration/NonRelConfiguration.h"
#include "NonRelInfo.h"
#include "Universal/MathConstant.h"

OccupationMap ConfigurationParser::ParseFractionalConfiguration(const std::string& configuration)
{
    OccupationMap ret;
    NonRelConfiguration config(configuration);

    // Split non-relativistic configuration into relativistic states
    NonRelConfiguration::iterator it = config.begin();
    while(it != config.end())
    {
        pOrbital s1, s2;

        NonRelInfo info(it->first);
        double occupancy = it->second;
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

        it++;
    }

    return ret;
}

std::vector<int> ConfigurationParser::ParseBasisSize(const std::string& basis)
{
    unsigned int p = 0;
    unsigned int L = 0;

    const char* basis_def = basis.c_str();
    std::vector<int> num_states;

    // Set num_states[L] to highest pqn
    while(basis_def[p])
    {
        int pqn = atoi(basis_def + p);
        while(basis_def[p] && isdigit(basis_def[p]))
            p++;
        while(basis_def[p] && !isdigit(basis_def[p]))
        {
            // Get L
            L = MathConstant::Instance()->GetL(tolower(basis_def[p]));
            if(L >= num_states.size())
                num_states.resize(L+1);
            num_states[L] = pqn;
            p++;
        }
    }

    return num_states;
}
