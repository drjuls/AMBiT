#include "ConfigurationParser.h"
#include "Include.h"
#include "Configuration/NonRelConfiguration.h"
#include "NonRelInfo.h"
#include "Universal/MathConstant.h"

namespace Ambit
{
OccupationMap ConfigurationParser::ParseFractionalConfiguration(const std::string& configuration)
{
    OccupationMap ret;
    auto config = ParseConfiguration<OrbitalInfo, double>(configuration);

    // Split non-relativistic configuration into relativistic states
    auto it = config.begin();
    while(it != config.end())
    {
        double occupancy = it->second;

        // Split non-rel info with too many electrons
        if((it->first.Kappa() > 0) && (occupancy > it->first.MaxNumElectrons()))
        {
            OrbitalInfo second_orbital(it->first.PQN(), -(it->first.L()+1));
            if(ret.GetOccupancy(second_orbital))
            {   *errstream << "ConfigurationParser::ParseFractionalConfiguration() initialised with problematic string " << configuration << std::endl;
                exit(1);
            }
            ret.insert(std::make_pair(second_orbital, occupancy - it->first.MaxNumElectrons()));
            occupancy = it->first.MaxNumElectrons();
        }

        // Check that orbital doesn't already exist and add
        if(ret.GetOccupancy(it->first))
        {   *errstream << "ConfigurationParser::ParseFractionalConfiguration() initialised with problematic string " << configuration << std::endl;
            exit(1);
        }

        ret.insert(std::make_pair(it->first, occupancy));

        it++;
    }

    // Sanity check
    for(const auto& orbital: ret)
    {   if(orbital.second > orbital.first.MaxNumElectrons())
        {   *errstream << "ConfigurationParser::ParseFractionalConfiguration() initialised with problematic string " << configuration << std::endl;
            exit(1);
        }
    }

    return ret;
}

NonRelInfo ConfigurationParser::ParseOrbital(const std::string& orbital)
{
    unsigned int p = 0;
    int L = -1;

    const char* basis_def = orbital.c_str();

    int pqn = atoi(basis_def + p);
    while(basis_def[p] && (isdigit(basis_def[p]) || isblank(basis_def[p])))
        p++;
    if(basis_def[p] && !isdigit(basis_def[p]) && !isblank(basis_def[p]))
    {
        L = MathConstant::Instance()->GetL(tolower(basis_def[p]));
    }
    if(L < 0)
    {   *errstream << "ConfigurationParser::ParseOrbital() bad string: " << orbital << std::endl;
        exit(1);
    }

    return NonRelInfo(pqn, L);
}

std::vector<int> ConfigurationParser::ParseBasisSize(const std::string& basis)
{
    unsigned int p = 0;
    int L = 0;

    const char* basis_def = basis.c_str();
    std::vector<int> num_states;

    // Set num_states[L] to highest pqn
    while(basis_def[p])
    {
        int pqn = atoi(basis_def + p);
        while(basis_def[p] && (isdigit(basis_def[p]) || isblank(basis_def[p])))
            p++;
        while(basis_def[p] && !isdigit(basis_def[p]) && !isblank(basis_def[p]))
        {
            // Get L
            L = MathConstant::Instance()->GetL(tolower(basis_def[p]));
            if(L < 0)
            {   *errstream << "ConfigurationParser::ParseBasisSize() bad string: " << basis << std::endl;
                exit(1);
            }
            else if(L >= num_states.size())
                num_states.resize(L+1);
            num_states[L] = pqn;
            p++;
        }
    }

    return num_states;
}
}
