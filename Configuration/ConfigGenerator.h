#ifndef CONFIG_GENERATOR_H
#define CONFIG_GENERATOR_H

#include "Basis/ExcitedStates.h"

#include "Configuration.h"
#include "RelativisticConfiguration.h"
#include "Projection.h"
#include "NonRelInfo.h"

class ConfigGenerator
{
public:
    ConfigGenerator(const ExcitedStates* manager)
    {   SetExcitedStates(manager);
    }
    ~ConfigGenerator(void) {}

    void SetExcitedStates(const ExcitedStates* manager);

    /** Generate all non-relativistic configurations possible by exciting num_excitations
        electrons. Append the new configurations to the list.
        All new configurations must have the same parity as the first.
        Returns a sorted, unique list.
     */
    void GenerateMultipleExcitations(ConfigList& configlist, unsigned int num_excitations);

    /** Divide electrons between partial waves to create all possible relativistic configurations
        from the given set of non-relativistic ones.
        PRE: nrlist should be unique.
        POST: rlist is sorted and unique.
      */
    void GenerateRelativisticConfigs(const ConfigList& nrlist, RelativisticConfigList& rlist);

    /** Make all projections of the given RelativisticConfigList that have a projection M = two_m/2.
        Remove configurations that cannot have this projection.
        PRE: nrlist should be unique.
     */
    void GenerateProjections(RelativisticConfigList& rlist, int two_m);

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
    std::set<RelativisticInfo> RelativisticSet;
    std::map<ElectronInfo, unsigned int> ElectronSet;
};

#endif
