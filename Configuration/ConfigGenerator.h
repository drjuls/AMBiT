#ifndef CONFIG_GENERATOR_H
#define CONFIG_GENERATOR_H

#include "Basis/ExcitedStates.h"

#include "Configuration.h"
#include "RelativisticConfiguration.h"
#include "Projection.h"
#include "HartreeFock/NonRelInfo.h"

class ConfigGenerator
{
    /** ConfigGenerator makes the set of projections for use in CI method.
        This includes a bunch of routines to create a set of
        non-relativistic configurations (eg: 3d2 4s, 4s2, ...),
        then generate the relativistic configurations and their projections.
        Additionally it stores the sets of configurations and the leading configurations
        for use by HamiltonianMatrix class.
     */
public:
    ConfigGenerator(const ExcitedStates* manager)
    {   SetExcitedStates(manager);
    }
    virtual ~ConfigGenerator(void) {}

    /** Reset the list of excited states. */
    void SetExcitedStates(const ExcitedStates* manager);

    /** Clear the non-relativistic and relativistic config lists. */
    void ClearConfigLists();

    /** Get list of leading configurations. */
    std::set<Configuration>* GetLeadingConfigs();
    const std::set<Configuration>* GetLeadingConfigs() const;

    /** Get list of non-relativistic configurations. */
    ConfigList* GetNonRelConfigs();
    const ConfigList* GetNonRelConfigs() const;

    /** Get list of relativistic configurations. */
    RelativisticConfigList* GetRelConfigs();
    const RelativisticConfigList* GetRelConfigs() const;

public:
    /** Add to the set of leading configurations. */
    void AddLeadingConfiguration(const Configuration& config);
    void AddLeadingConfigurations(const std::set<Configuration> config_set);
    
    /** Generate all non-relativistic configurations possible by exciting num_excitations
        electrons. Append the new configurations to the list.
        All new configurations must have the same parity as the first.
        Returns a sorted, unique list and appends that list to local nrlist.
     */
    void GenerateMultipleExcitations(ConfigList& configlist, unsigned int num_excitations, Parity parity);

    /** Call GenerateMultipleExcitations() with leading configurations as input. */
    void GenerateMultipleExcitationsFromLeadingConfigs(unsigned int num_electrons, Parity parity);

    /** Divide electrons between partial waves to create all possible relativistic configurations
        from the set of non-relativistic ones.
        PRE: nrlist should be unique.
        POST: rlist is sorted and unique.
      */
    void GenerateRelativisticConfigs();

    /** Make all projections of the rlist that have a projection M = two_m/2.
        Remove configurations that cannot have this projection.
        PRE: rlist should be unique.
     */
    void GenerateProjections(int two_m);

protected:
    /** Generate all non-relativistic configurations possible by exciting one electron
        of the original list. Append the new configurations to the list.
     */
    void GenerateExcitations(ConfigList& configlist);

    /** Split the current NonRelInfo of config. Recursively split the rest.
        When config.AtEnd(), add it to rlist.
      */
    void SplitNonRelInfo(Configuration config, RelativisticConfigList& rlist);

protected:
    const ExcitedStates* states;

    std::set<NonRelInfo> NonRelSet;
    std::set<StateInfo> RelativisticSet;
    std::map<ElectronInfo, unsigned int> ElectronSet;
    
    std::set<Configuration> leading_configs;
    ConfigList nrlist;
    RelativisticConfigList rlist;
};

#endif
