#ifndef CONFIG_FILE_GENERATOR_H
#define CONFIG_FILE_GENERATOR_H

#include "ConfigGenerator.h"

//typedef std::map<Configuration, double> ConfigPercentage;

class ConfigFileGenerator : public ConfigGenerator
{
public:
    ConfigFileGenerator(const ExcitedStates* manager, const std::string& atom_identifier, Symmetry sym):
        ConfigGenerator(manager, atom_identifier, sym)
    {}
    virtual ~ConfigFileGenerator(void) {}

    /** Clear the non-relativistic and relativistic config lists and free memory.
        Also clear InputContributions and Contributions.
        Use Read() (and ReadConfigs() if wanted) to restore.
     */
    virtual void Clear();

    void SetInputFile(const std::string& filename)
    {   input_filename = filename;
    }

    void SetOutputFile(const std::string& filename)
    {   output_filename = filename;
    }
    
    /** Read configurations from file with percentages larger than cutoff. */
    void ReadConfigs(double cutoff = 0.0);
    
    /** Write configurations with percentages larger than cutoff to file. */
    void WriteConfigs(double cutoff = 0.0);

    /** Map non-rel configurations to percentages. */
    void AddPercentages(const std::map<Configuration, double> percentages);

    /** Generate all non-relativistic configurations possible by exciting num_excitations
        electrons from existing configurations with percentages larger then the cutoff, as
        well as the leading configurations.
        Append the new configurations to the list.
     */
    void GenerateMultipleExcitationsFromImportantConfigs(unsigned int num_excitations, double cutoff);

protected:
    std::string input_filename;
    std::string output_filename;

    // Map of input configurations and their percentage contributions to eigenstates
    std::map<Configuration, double> InputContributions;

    // Map of final configurations and their percentage contributions to eigenstates
    std::map<Configuration, double> Contributions;
};

#endif
