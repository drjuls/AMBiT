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

    /** Divide electrons between partial waves to create all possible relativistic configurations
        from the set of non-relativistic ones.
        PRE: nrlist should be unique.
        POST: rlist is sorted and unique.
     */
    pRelativisticConfigList GenerateRelativisticConfigurations(pConfigList nrlist) const;

    /** Generate relativistic configurations from nrlist. Remove configs of wrong parity.
        If angular_library is supplied, returned list includes projections and CSFs with M = J.
     */
    pRelativisticConfigList GenerateRelativisticConfigurations(pConfigList nrlist, const Symmetry& sym, pAngularDataLibrary angular_library = nullptr) const;

    /** Generate non-relativistic and then relativistic configurations based on the input file.
        If angular_library is supplied, returned list includes projections and CSFs with M = J.
     */
    pRelativisticConfigList GenerateRelativisticConfigurations(const Symmetry& sym, pAngularDataLibrary angular_library = nullptr);

    /** Make all projections of the rlist that have a projection M = two_m/2.
        Remove configurations that cannot have this projection.
        Generate CSFs with angular momentum two_j/2.
        PRE: rlist should be unique.
     */
    void GenerateProjections(pRelativisticConfigList rlist, const Symmetry& sym, int two_m, pAngularDataLibrary angular_library) const;

    /** Generate configuration state functions (CSFs) with J = M = two_m/2. */
    void GenerateProjections(pRelativisticConfigList rlist, int two_m, pAngularDataLibrary angular_library) const
    {   Parity p = rlist->front().GetParity();
        GenerateProjections(rlist, Symmetry(two_m, p), two_m, angular_library);
    }

protected:
    /** Generate non-relativistic configs based on the current input file section.
        This may include limiting the non-relativistic configurations based on configuration average energies,
        provided one_body and two_body operators are provided.
     */
    pConfigList ParseAndGenerateNonRelConfigurations(pHFOperator one_body = nullptr, pHartreeY two_body = nullptr);

    /** Generate all non-relativistic configurations possible by exciting one electron
        of the original list. Append the new configurations to the list.
     */
    void GenerateExcitations(pConfigList configlist, const NonRelInfoSet& electron_valence, const NonRelInfoSet& hole_valence) const;

protected:
    // Inputs
    MultirunOptions& user_input;
    pOrbitalManagerConst orbitals;

    // Set of all valence states
    NonRelInfoSet NonRelSet;

    pConfigList leading_configs;
};

#endif
