#ifndef CONFIG_GENERATOR_H
#define CONFIG_GENERATOR_H

#include "Atom/MultirunOptions.h"
#include "HartreeFock/StateManager.h"
#include "Configuration.h"
#include "RelativisticConfiguration.h"
#include "Projection.h"
#include "HartreeFock/NonRelInfo.h"
#include "Symmetry.h"

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

    void SetTwoM(int two_m) { TwoM = two_m; }
    int GetTwoM() const { return TwoM; }
    void SetParity(Parity P) { parity = P; }
    Parity GetParity() const { return parity; }

    /** Clear the non-relativistic and relativistic config lists and free memory.
        Use Read() to restore.
     */
    virtual void Clear();

    /** Get list of leading configurations. */
    std::set<Configuration>* GetLeadingConfigs();
    const std::set<Configuration>* GetLeadingConfigs() const;

    /** Get list of non-relativistic configurations. */
    virtual ConfigList* GetNonRelConfigs();

    /** Get list of relativistic configurations. */
    virtual pRelativisticConfigList GetRelConfigs() {  return rlist; }
    virtual pRelativisticConfigListConst GetRelConfigs() const { return rlist; }

public:
    /** Generate non-relativistic and then relativistic configurations based on the input file. */
    pRelativisticConfigList GenerateRelativisticConfigurations();

    /** Add to the set of leading configurations. */
    DEPRECATED void AddLeadingConfiguration(const Configuration& config);
    DEPRECATED void AddLeadingConfigurations(const std::set<Configuration> config_set);
    DEPRECATED void AddNonRelConfiguration(const Configuration& config);

    /** Generate all non-relativistic configurations possible by exciting num_excitations
        electrons. Append the new configurations to the list.
        All new configurations must have the same parity as the first.
        Returns a sorted, unique list and appends that list to local nrlist.
     */
    void GenerateMultipleExcitations(ConfigList& configlist, unsigned int num_excitations, const NonRelInfoSet* AllowedStateSet = NULL);

    /** Call GenerateMultipleExcitations() with leading configurations as input. */
    void GenerateMultipleExcitationsFromLeadingConfigs(unsigned int num_excitations, const NonRelInfoSet* AllowedStateSet = NULL);

    /** Divide electrons between partial waves to create all possible relativistic configurations
        from the set of non-relativistic ones.
        PRE: nrlist should be unique.
        POST: rlist is sorted and unique.
      */
    void GenerateRelativisticConfigs();

    /** Make all projections of the rlist that have a projection M = two_m/2.
        Remove configurations that cannot have this projection. If two_m is not given,
        the default is M = J.
        PRE: rlist should be unique.
     */
//    void GenerateProjections(int two_m);
//    void GenerateProjections();

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
    void GenerateExcitations(ConfigList& configlist, const NonRelInfoSet* states_to_be_excited_to = NULL);

    /** Split the current NonRelInfo of config. Recursively split the rest.
        When config.AtEnd(), add it to rlist.
      */
    void SplitNonRelInfo(const Configuration& config, Configuration::const_iterator current_orbital, RelativisticConfiguration& relconfig, RelativisticConfigList& rlist);

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
    Parity parity;
    int TwoM;

    // Set of all valence states
    NonRelInfoSet NonRelSet;

    std::set<Configuration> leading_configs;
    ConfigList nrlist;
    pRelativisticConfigList rlist;
};

#endif
