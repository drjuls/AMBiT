#ifndef NON_REL_CONFIGURATION_H
#define NON_REL_CONFIGURATION_H

#include "HartreeFock/Configuration.h"
#include "HartreeFock/OrbitalInfo.h"
#include "HartreeFock/NonRelInfo.h"
#include "RelativisticConfigList.h"
#include "HartreeFock/OrbitalMap.h"
#include "HartreeFock/HFOperator.h"
#include "HartreeFock/HartreeY.h"
#include "SortedList.h"
#include <set>

/** Non-relativistic configuration does not know whether an orbital is above or below
    the Fermi surface (hole or particle state) so the only requirement is that
    abs(occupation) doesn't exceed allowed occupation 2.abs(kappa).
 */
class NonRelConfiguration: public Configuration<NonRelInfo, int>
{
public:
    NonRelConfiguration() = default;
    NonRelConfiguration(const BaseConfiguration& other): BaseConfiguration(other) {}
    NonRelConfiguration(BaseConfiguration&& other): BaseConfiguration(other) {}
    NonRelConfiguration(const RelativisticConfiguration& other);
    NonRelConfiguration(const std::string& name);
    virtual ~NonRelConfiguration() = default;

public:
    /** Return success. */
    bool RemoveSingleParticle(const NonRelInfo& info);
    /** Return success. */
    bool AddSingleParticle(const NonRelInfo& info);

    virtual std::string Name(bool aSpaceFirst = true) const;
    virtual std::string NameNoSpaces() const;

    /** Generate RelativisticConfigurations by distributing electrons among relativistic orbitals.
        Returns public member relconfiglist.
     */
    pRelativisticConfigList GenerateRelativisticConfigs();

    /** Calculate configuration average energy using relconfiglist. */
    double CalculateConfigurationAverageEnergy(pOrbitalMapConst orbitals, pHFOperator one_body, pHartreeY two_body);

    /** Calculate the largest projection possible for this configuration. */
    int GetTwiceMaxProjection() const;

    /** Calculate number of levels (all symmetries and relativistic configurations) corresponding to this configuration.
     */
    int GetNumberOfLevels();

    void Read(FILE* fp);        //!< Read configuration
    void Write(FILE* fp) const; //!< Write configuration only

public:
    pRelativisticConfigList relconfiglist;

protected:
    /** Split the current NonRelInfo of config. Recursively split the rest.
        When at end(), add it to relconfiglist.
     */
    void SplitNonRelInfo(NonRelConfiguration::const_iterator current_orbital, RelativisticConfiguration& relconfig);
};

/** ConfigList has a list and an unsigned int to indicate a small cut-off (Nsmall). */
typedef std::pair<std::vector<NonRelConfiguration>, unsigned int> ConfigList;
std::ostream& operator<<(std::ostream& stream, const ConfigList& config_list);

typedef std::shared_ptr<ConfigList> pConfigList;
typedef std::shared_ptr<const ConfigList> pConfigListConst;

/** Sort config_list and remove duplicates.
    Preserve separation between (0 -> Nsmall) and (Nsmall -> config_list->size()).
    Nsmall may change if elements are removed.
 */
void SortAndUnique(pConfigList& config_list);

#endif
