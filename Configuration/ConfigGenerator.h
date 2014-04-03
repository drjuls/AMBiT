#ifndef CONFIG_GENERATOR_H
#define CONFIG_GENERATOR_H

#include "Atom/MultirunOptions.h"
#include "HartreeFock/StateManager.h"
#include "NonRelConfiguration.h"
#include "RelativisticConfiguration.h"
#include "Projection.h"
#include "HartreeFock/NonRelInfo.h"
#include "Symmetry.h"
#include <set>

/** ConfigGenerator makes the set of projections for use in CI method.
    This includes a bunch of routines to create a set of
    non-relativistic configurations (eg: 3d2 4s, 4s2, ...),
    then generate the relativistic configurations and their projections.
    Additionally it stores the sets of configurations and the leading configurations
    for use by HamiltonianMatrix class.
 */
class ConfigGenerator
{
public:
    ConfigGenerator(pStateManagerConst coreStates, pStateManagerConst excitedStates, MultirunOptions& userInput);
    virtual ~ConfigGenerator();

    /** Clear the non-relativistic and relativistic config lists and free memory.
        Use Read() to restore.
     */
    virtual void Clear();

    /** Get list of leading configurations. */
    pConfigList GetLeadingConfigs();
    pConfigListConst GetLeadingConfigs() const;

    /** Generate non-relativistic and then relativistic configurations based on the input file.
        Returned list has CSFs with J = M.
     */
    pRelativisticConfigList GenerateRelativisticConfigurations(Parity parity, int two_j);

protected:
    /** Generate all non-relativistic configurations possible by exciting num_excitations
        electrons. Append the new configurations to the list, keeping it sorted and unique.
        All new configurations must have the same parity as the first.
     */
    void GenerateMultipleExcitations(ConfigList& configlist, unsigned int num_excitations, const NonRelInfoSet* AllowedStateSet = NULL) const;

    /** Call GenerateMultipleExcitations() with leading configurations as input. */
    ConfigList GenerateMultipleExcitationsFromLeadingConfigs(unsigned int num_excitations, const NonRelInfoSet* AllowedStateSet = NULL) const;

    /** Divide electrons between partial waves to create all possible relativistic configurations
        from the set of non-relativistic ones.
        PRE: nrlist should be unique.
        POST: rlist is sorted and unique.
      */
    pRelativisticConfigList GenerateRelativisticConfigs(const ConfigList& nrlist) const;

    /** Make all projections of the rlist that have a projection M = two_m/2.
        Remove configurations that cannot have this projection.
        PRE: rlist should be unique.
     */
    void GenerateProjections(pRelativisticConfigList rlist, int two_m) const;

    /** Store the leading_configs and rlist (RelativisticConfigList).
        Filename is "identifier.twoJ.P.configs".
     */
    DEPRECATED void Write() const;

    /** Read leading_configs and rlist. Return success. */
    DEPRECATED bool Read();

protected:
    /** Generate all non-relativistic configurations possible by exciting one electron
        of the original list. Append the new configurations to the list.
     */
    void GenerateExcitations(ConfigList& configlist, const NonRelInfoSet* states_to_be_excited_to = NULL) const;

    /** Split the current NonRelInfo of config. Recursively split the rest.
        When config.AtEnd(), add it to rlist.
      */
    void SplitNonRelInfo(const NonRelConfiguration& config, NonRelConfiguration::const_iterator current_orbital, RelativisticConfiguration& relconfig, pRelativisticConfigList& rlist) const;

    /** Restore nrlist from rlist.
        Usually we no longer want the nrlist after a read, but if it is wanted,
        nrlist is generated from rlist on the fly.
     */
    void RestoreNonRelConfigs();

protected:
    // Inputs
    MultirunOptions& user_input;
    pStateManagerConst core;
    pStateManagerConst excited;

    // Set of all valence states
    NonRelInfoSet NonRelSet;

    pConfigList leading_configs;
};

#endif
