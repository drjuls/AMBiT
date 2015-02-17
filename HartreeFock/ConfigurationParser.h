#ifndef CONFIGURATION_PARSER_H
#define CONFIGURATION_PARSER_H

#include "Include.h"
#include "Core.h"
#include "Configuration.h"
#include "Universal/MathConstant.h"
#include <sstream>

class ConfigurationParser
{
public:
    template <class OrbitalType, class OccupancyType>
    static Configuration<OrbitalType, OccupancyType> ParseConfiguration(const std::string& configuration);

    static OccupationMap ParseFractionalConfiguration(const std::string& configuration);

    // Set num_states[L] to highest pqn. 0 indicates that the L wasn't specified in the string.
    static std::vector<int> ParseBasisSize(const std::string& basis);
};


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

            // Peek extra char: + or -
            char next_char = sstream.peek();
            if(!isdigit(next_char))
            {
                if(next_char == '+')
                {   sstream.get(temp_char);
                    kappa = -(L+1);
                }
                else if(next_char != '-')
                {   sstream.get(temp_char);
                }
            }

            // Get occupancy
            sstream >> occupancy;

            // Check everything is okay and skip any whitespace
            if((L >= 10) || (pqn < L+1) || (occupancy > OccupancyType(4*L + 2)))
            {   throw (1);
            }

            ret[OrbitalType(pqn, kappa)] = occupancy;
            while(sstream.good() && isspace(sstream.peek()))
                sstream.get();
        }
    }
    catch (...)
    {   *errstream << "Configuration() initialised with problematic string " << configuration << std::endl;
        exit(1);
    }
    
    return ret;
}

#endif
