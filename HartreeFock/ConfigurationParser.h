#ifndef CONFIGURATION_PARSER_H
#define CONFIGURATION_PARSER_H

#include "Include.h"
#include "Core.h"
#include "Configuration.h"
#include "Universal/MathConstant.h"
#include <sstream>

class ConfigurationParser
{
    /** Class that holds static functions to parse strings to configurations. */
public:
    /** Convert orbital name (e.g. "5s" or "4p+") to OrbitalType. */
    template <class OrbitalType>
    static OrbitalType ParseOrbital(const std::string& orbital);

    /** Get Configuration<OrbitalType, OccupancyType> from string. */
    template <class OrbitalType, class OccupancyType>
    static Configuration<OrbitalType, OccupancyType> ParseConfiguration(const std::string& configuration);

    /** Split non-relativistic configurations from string according to occupancy weighting. */
    static OccupationMap ParseFractionalConfiguration(const std::string& configuration);

    /** Set num_states[L] to highest pqn. 0 indicates that the L wasn't specified in the string.
        e.g. "6sp5d4f" -> [6, 6, 5, 4]
     */
    static std::vector<int> ParseBasisSize(const std::string& basis);
};

template <class OrbitalType>
OrbitalType ConfigurationParser::ParseOrbital(const std::string& orbital)
{
    unsigned int p = 0;
    int L = -1;
    int kappa = 0;

    const char* orbital_def = orbital.c_str();
    // Get pqn
    int pqn = atoi(orbital_def);
    while(orbital_def[p] && (isdigit(orbital_def[p]) || isblank(orbital_def[p])))
        p++;

    // Get L
    if(orbital_def[p])
    {   L = MathConstant::Instance()->GetL(tolower(orbital_def[p]));
        p++;
    }
    kappa = (L==0? -1: L);

    // Check for "+"
    if(orbital_def[p] == '+')
    {
        kappa = -(L+1);
    }

    if((L == -1) || (pqn < L+1))
    {   *errstream << "ConfigurationParser::ParseOrbital() initialised with problematic string " << orbital << std::endl;
        exit(1);
    }

    return OrbitalType(pqn, kappa);
}

template <class OrbitalType, class OccupancyType>
Configuration<OrbitalType, OccupancyType> ConfigurationParser::ParseConfiguration(const std::string& configuration)
{
    Configuration<OrbitalType, OccupancyType> ret;
    std::stringstream sstream(configuration);

    char temp_char;
    int pqn, L, kappa;
    OccupancyType occupancy;

    try{
        // Skip whitespace
        while(sstream.good() && isspace(sstream.peek()))
            sstream.get();

        // Special case: 0 = vacuum.
        if(sstream.good() && sstream.peek() == '0')
            return ret;

        while(sstream.good())
        {
            // Get pqn
            sstream >> pqn;

            // Get L
            sstream.get(temp_char);
            L = MathConstant::Instance()->GetL(temp_char);
            kappa = L;

            // Peek extra char; check if +
            char next_char = sstream.peek();
            if(!isdigit(next_char))
            {
                if(next_char == '+')
                {   sstream.get(temp_char);
                    kappa = -(L+1);
                }
            }
            // s-wave can come without + (implied)
            if(kappa == 0)
            {   kappa = -1;
            }

            // Get occupancy
            sstream >> occupancy;

            // Check everything is okay and skip any whitespace
            OrbitalType info(pqn, kappa);
            if((L == -1) || (pqn < L+1) || (occupancy > info.MaxNumElectrons()))
            {   throw (1);
            }

            ret[info] = occupancy;
            while(sstream.good() && isspace(sstream.peek()))
                sstream.get();
        }
    }
    catch (...)
    {   *errstream << "ConfigurationParser::ParseConfiguration() initialised with problematic string " << configuration << std::endl;
        exit(1);
    }
    
    return ret;
}

#endif
