#include "ConfigurationParser.h"
#include "Include.h"
#include "Configuration/Configuration.h"
#include "NonRelInfo.h"
#include "Universal/MathConstant.h"

Configuration ConfigurationParser::ParseConfiguration(const std::string& configuration)
{
    Configuration ret;

    // Set num_states[L] to highest pqn
    unsigned int p = 0;
    const char* all_states = configuration.c_str();
    unsigned int pqn, L, occupancy;

    // Skip whitespace
    while(all_states[p] && isspace(all_states[p]))
        p++;

    while(all_states[p])
    {
        // Get pqn
        pqn = atoi(all_states + p);
        while(all_states[p] && isdigit(all_states[p]))
            p++;
        // Get L
        if(all_states[p])
        {   // Get L
            L = MathConstant::Instance()->GetL(all_states[p]);
            p++;
        }
        // Get occupancy
        occupancy = atoi(all_states + p);
        while(all_states[p] && isdigit(all_states[p]))
            p++;

        // Check everything is okay and skip any whitespace
        if((L >= 10) || (pqn < L+1) || (occupancy > 4*L + 2))
        {   *errstream << "Configuration() initialised with problematic string " << configuration << std::endl;
            exit(1);
        }
        ret[NonRelInfo(pqn, L)] = occupancy;
        while(all_states[p] && isspace(all_states[p]))
            p++;
    }

    return ret;
}

OccupationMap ConfigurationParser::ParseFractionalConfiguration(const std::string& configuration)
{
    OccupationMap ret;
    Configuration config(configuration);

    // Split non-relativistic configuration into relativistic states
    Configuration::iterator it = config.begin();
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
