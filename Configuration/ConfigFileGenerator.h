#ifndef CONFIG_FILE_GENERATOR_H
#define CONFIG_FILE_GENERATOR_H

#include "ConfigGenerator.h"

//typedef std::map<Configuration, double> ConfigPercentage;

class ConfigFileGenerator : public ConfigGenerator
{
public:
    ConfigFileGenerator(const ExcitedStates* manager):
        ConfigGenerator(manager)
    {}
    virtual ~ConfigFileGenerator(void) {}
    
    void SetInputFile(const std::string& filename)
    {   input_filename = filename;
    }

    void SetOutputFile(const std::string& filename)
    {   output_filename = filename;
    }
    
    /** Read configurations from file with percentages larger than cutoff. */
    void ReadConfigs(Parity parity, double cutoff = 0.0);
    
    /** Write configurations with percentages larger than cutoff to file. */
    void WriteConfigs(double cutoff = 0.0);

    /** Map non-rel configurations to percentages. */
    void AddPercentages(const std::map<Configuration, double> percentages);

protected:
    std::string input_filename;
    std::string output_filename;

    // Map of final configurations and their percentage contributions to eigenstates
    std::map<Configuration, double> Contributions;
};

#endif
