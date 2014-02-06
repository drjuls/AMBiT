#ifndef CONFIGURATION_PARSER_H
#define CONFIGURATION_PARSER_H

#include "OrbitalInfo.h"
#include <map>
#include <iostream>

//typedef std::map<OrbitalInfo, unsigned int> Configuration;
typedef std::map<OrbitalInfo, double> OccupationMap;

class ConfigurationParser
{
public:
//    Configuration ParseConfiguration(const std::string& config);
    static OccupationMap ParseFractionalConfiguration(const std::string& configuration);

    // Set num_states[L] to highest pqn
    static std::vector<int> ParseBasisSize(const std::string& basis);
};

#endif
