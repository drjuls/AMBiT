#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "HartreeFock/OrbitalInfo.h"
#include "HartreeFock/NonRelInfo.h"
#include "RelativisticConfiguration.h"
#include "SortedList.h"
#include <map>
#include <set>
#include <string>

class Configuration
{
    /** Configuration is a map between a set of orbitals and their (integer) occupancy.
        It can be used for relativistic and non-relativistic states.
     */
public:
    Configuration() {}
    Configuration(const Configuration& other): config(other.config) {}
    Configuration(const RelativisticConfiguration& other);
    Configuration(const std::string& name);
    virtual ~Configuration() {}

public:
    typedef std::map<NonRelInfo, int>::iterator iterator;
    typedef std::map<NonRelInfo, int>::const_iterator const_iterator;

public:
    /** Get occupancy of a particular single particle state (zero if absent). */
    int GetOccupancy(const NonRelInfo& info) const;
    int& operator[](const NonRelInfo& info);

    iterator begin() { return config.begin(); }
    const_iterator begin() const { return config.begin(); }
    iterator end() { return config.end(); }
    const_iterator end() const { return config.end(); }
    iterator find(const NonRelInfo& info);
    const_iterator find(const NonRelInfo& info) const;

    void clear() { config.clear(); }
    bool empty() const { return config.empty(); }
    int size() const { return config.size(); }

    iterator erase(const_iterator position);
    int erase(const NonRelInfo& info);
    iterator erase(const_iterator first, const_iterator last);

    bool RemoveSingleParticle(const NonRelInfo& info);
    bool AddSingleParticle(const NonRelInfo& info);

    /** Number of particles = number of electrons - number of holes. */
    virtual int ParticleNumber() const;

    /** Excitation number = number of electrons + number of holes. */
    virtual int ExcitationNumber() const;
    Parity GetParity() const;

    bool operator<(const Configuration& other) const;
    bool operator==(const Configuration& other) const;
    virtual std::string Name(bool aSpaceFirst = true) const;
    virtual std::string ShortName() const;

    // File storage (binary)
    virtual void Write(FILE* fp) const;
    virtual void Read(FILE* fp);

protected:
    /** Map single particle state to occupancy. */
    std::map<NonRelInfo, int> config;
};

class ConfigList : public SortedList<Configuration>
{
public:
    ConfigList() {}
    ConfigList(const ConfigList& other): BaseSortedList(other) {}
    ConfigList(ConfigList&& other): BaseSortedList(other) {}
    ConfigList(const Configuration& val): BaseSortedList(val) {}
    ConfigList(Configuration&& val): BaseSortedList(val) {}

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

class ConfigurationPair : public std::pair<Configuration, double>
{
public:
    ConfigurationPair(const Configuration& aConfiguration, const double& aDouble);
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
    Configuration GetLargestConfiguration() { return begin()->first; }
    void Print();
};

#endif
