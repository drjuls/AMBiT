#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "HartreeFock/OrbitalInfo.h"
#include <map>
#include <set>
#include <list>
#include <string>

#ifndef PARITY_ENUM
#define PARITY_ENUM
    enum Parity { even, odd };
#endif

class Configuration
{
    /** Configuration is a map between a set of single particle states and their occupancy.
        It can be used for relativistic and non-relativistic states.
     */
public:
    Configuration() { First(); }
    Configuration(const Configuration& other);
    Configuration(const std::string& name);
    virtual ~Configuration(void) {}

    // Controllers for built-in state iterator
    inline void First() const;
    inline void Next() const;
    inline bool AtEnd() const;
    void SetIterator(const OrbitalInfo& info) const;

    OrbitalInfo GetInfo() const;
    unsigned int GetOccupancy() const;
    void SetOccupancy(unsigned int occupancy);

    // Get occupancy of a particular single particle state.
    // (zero if absent)
    unsigned int GetOccupancy(const OrbitalInfo& info) const;

    // These return the success of the operation.
    virtual bool RemoveSingleParticle(const OrbitalInfo& info);
    virtual bool AddSingleParticle(const OrbitalInfo& info);
    virtual bool SetOccupancy(const OrbitalInfo& info, unsigned int occupancy);

    void Clear();
    bool Empty() const;

    virtual unsigned int NumParticles() const;

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
    std::map<OrbitalInfo, unsigned int> Config;
    mutable std::map<OrbitalInfo, unsigned int>::const_iterator it;
};

class ConfigList : public std::list<Configuration>
{
public:
    void Print();
};

inline void Configuration::First() const
{
    it = Config.begin();
}

inline void Configuration::Next() const
{   
    if(!AtEnd())
        it++;
}

inline bool Configuration::AtEnd() const
{   
    return (it == Config.end());
}

inline void Configuration::Clear()
{
    Config.clear();
}

inline bool Configuration::Empty() const
{
    return Config.empty();
}

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

class ConfigurationSet : public std::set<ConfigurationPair, ConfigurationPairCompare>
{
public:
    ConfigurationSet::iterator GetLargestConfiguration() { return begin(); }
    void Print();
};

#endif
