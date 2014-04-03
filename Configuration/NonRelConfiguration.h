#ifndef NON_REL_CONFIGURATION_H
#define NON_REL_CONFIGURATION_H

#include "HartreeFock/Configuration.h"
#include "HartreeFock/OrbitalInfo.h"
#include "HartreeFock/NonRelInfo.h"
#include "RelativisticConfiguration.h"
#include "SortedList.h"
#include <set>

class NonRelConfiguration: public Configuration<NonRelInfo, int>
{
    /** Configuration is a map between a set of orbitals and their (integer) occupancy.
        It can be used for relativistic and non-relativistic states.
     */
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
    bool RemoveSingleParticle(const NonRelInfo& info);
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

typedef boost::shared_ptr<ConfigList> pConfigList;
typedef boost::shared_ptr<const ConfigList> pConfigListConst;

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
