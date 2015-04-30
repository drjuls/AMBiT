#ifndef CONFIG_GENERATOR_H
#define CONFIG_GENERATOR_H

#include "Atom/MultirunOptions.h"
#include "Basis/OrbitalManager.h"
#include "NonRelConfiguration.h"
#include "RelativisticConfigList.h"
#include "Projection.h"
#include "HartreeFock/NonRelInfo.h"
#include "Symmetry.h"
#include <set>

typedef std::set<NonRelInfo> NonRelInfoSet;

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
    ConfigGenerator(pOrbitalManagerConst orbitals, MultirunOptions& userInput);
    virtual ~ConfigGenerator() {}

    /** Get list of leading configurations. */
    pConfigList GetLeadingConfigs();
    pConfigListConst GetLeadingConfigs() const;

    /** Generate non-relativistic configs based on the input file.
        This may include limiting the non-relativistic configurations based on configuration average energies,
        provided one_body and two_body operators are provided.
     */
    pConfigList GenerateNonRelConfigurations(pHFOperator one_body = nullptr, pHartreeY two_body = nullptr);

    /** Generate non-relativistic and then relativistic configurations based on the input file.
        If(generate_projections == true), returned list includes projections and CSFs with J = M.
     */
    pRelativisticConfigList GenerateRelativisticConfigurations(const Symmetry& sym, bool generate_projections = true);

    /** Make all projections of the rlist that have a projection M = two_m/2.
        Remove configurations that cannot have this projection.
        Generate CSFs with angular momentum two_j/2.
        PRE: rlist should be unique.
     */
    void GenerateProjections(pRelativisticConfigList rlist, int two_m, int two_j) const;

    /** Generate configuration state functions (CSFs) with J = M = two_m/2.
        Equivalent to GenerateProjections(rlist, two_m, two_m).
     */
    void GenerateProjections(pRelativisticConfigList rlist, int two_m) const
    {   GenerateProjections(rlist, two_m, two_m);
    }

    /** For each non-relativistic configuration in nrlist, generate relconfiglist
        with projections and CSFs with correct symmetry.
        Remove all relativistic configs that do not have a CSF with the correct symmetry.
        For non-relativistic configurations that have no CSFs satisfying the symmetry,
        relconfiglist is nullptr.
     */
    void GenerateProjections(pConfigList nrlist, const Symmetry& sym, int two_m) const;

protected:
    /** Generate all non-relativistic configurations possible by exciting one electron
        of the original list. Append the new configurations to the list.
     */
    void GenerateExcitations(pConfigList configlist, const NonRelInfoSet& electron_valence, const NonRelInfoSet& hole_valence) const;

    /** Divide electrons between partial waves to create all possible relativistic configurations
        from the set of non-relativistic ones.
        PRE: nrlist should be unique.
        POST: rlist is sorted and unique.
     */
    pRelativisticConfigList GenerateRelativisticConfigs(pConfigList nrlist) const;

protected:
    // Inputs
    MultirunOptions& user_input;
    pOrbitalManagerConst orbitals;

    // Set of all valence states
    NonRelInfoSet NonRelSet;

    pConfigList leading_configs;
};

#endif
