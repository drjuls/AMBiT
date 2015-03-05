#ifndef NON_REL_CONFIGURATION_H
#define NON_REL_CONFIGURATION_H

#include "HartreeFock/Configuration.h"
#include "HartreeFock/OrbitalInfo.h"
#include "HartreeFock/NonRelInfo.h"
#include "RelativisticConfigList.h"
#include "SortedList.h"
#include <set>

/** Non-relativistic configuration does not know whether an orbital is above or below
    the Fermi surface (hole or particle state) so the only requirement is that
    abs(occupation) doesn't exceed allowed occupation 2.abs(kappa).
 */
class NonRelConfiguration: public Configuration<NonRelInfo, int>
{
public:
    NonRelConfiguration() {}
    NonRelConfiguration(const BaseConfiguration& other): BaseConfiguration(other) {}
    NonRelConfiguration(BaseConfiguration&& other): BaseConfiguration(other) {}
    NonRelConfiguration(const NonRelConfiguration& other): BaseConfiguration(other) {}
    NonRelConfiguration(NonRelConfiguration&& other): BaseConfiguration(other) {}
    NonRelConfiguration(const RelativisticConfiguration& other);
    NonRelConfiguration(const std::string& name);
    virtual ~NonRelConfiguration() {}

    NonRelConfiguration& operator=(const NonRelConfiguration& other) { BaseConfiguration::operator=(other); return *this; }
    NonRelConfiguration& operator=(NonRelConfiguration&& other) { BaseConfiguration::operator=(other); return *this; }

public:
    /** Return success. */
    bool RemoveSingleParticle(const NonRelInfo& info);
    /** Return success. */
    bool AddSingleParticle(const NonRelInfo& info);

    virtual std::string Name(bool aSpaceFirst = true) const;
    virtual std::string ShortName() const;
};

class ConfigList : public SortedList<NonRelConfiguration>
{
public:
    ConfigList() {}
    ConfigList(const ConfigList& other): BaseSortedList(other) {}
    ConfigList(ConfigList&& other): BaseSortedList(other) {}
    ConfigList(const NonRelConfiguration& val): BaseSortedList(val) {}
    ConfigList(NonRelConfiguration&& val): BaseSortedList(val) {}

    /** Create unique non-relativistic list from rlist. */
    ConfigList(const RelativisticConfigList& rlist);

    ConfigList& operator=(const ConfigList& other)
    {   BaseSortedList::operator=(other);
        return *this;
    }
    ConfigList& operator=(ConfigList&& other)
    {   BaseSortedList::operator=(other);
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& stream, const ConfigList& config_list);
};

typedef std::shared_ptr<ConfigList> pConfigList;
typedef std::shared_ptr<const ConfigList> pConfigListConst;

class ConfigurationPair : public std::pair<NonRelConfiguration, double>
{
public:
    ConfigurationPair(const NonRelConfiguration& aConfiguration, const double& aDouble);
};

class ConfigurationPairCompare 
{
public:
    bool operator() (const ConfigurationPair& lhs, const ConfigurationPair& rhs) const
    {   return lhs.second > rhs.second;
    }
};

// Stores a set of configurations and their composition numbers sorted by descending order
class ConfigurationSet : public std::set<ConfigurationPair, ConfigurationPairCompare>
{
public:
    ConfigurationSet::iterator GetLargestConfigurationPair() { return begin(); }
    NonRelConfiguration GetLargestConfiguration() { return begin()->first; }
    void Print();
};

#endif
