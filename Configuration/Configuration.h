#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "SingleParticleInfo.h"
#include <map>
#include <list>
#include <string>

#ifndef PARITY_ENUM
#define PARITY_ENUM
    enum Parity { even, odd };
#endif

class Configuration
{
    /** Configuration is really a non-relativistic list, but could be used for
        relativistic states as well.
     */
public:
    Configuration() { First(); }
    Configuration(const Configuration& other);
    virtual ~Configuration(void) {}

    inline void First() const;
    inline void Next() const;
    inline bool AtEnd() const;
    void SetIterator(const SingleParticleInfo& info) const;

    SingleParticleInfo GetInfo() const;
    unsigned int GetOccupancy() const;
    void SetOccupancy(unsigned int occupancy);

    // These return the success of the operation.
    virtual bool RemoveSingleParticle(const SingleParticleInfo& info);
    virtual bool AddSingleParticle(const SingleParticleInfo& info);
    virtual bool SetOccupancy(const SingleParticleInfo& info, unsigned int occupancy);

    virtual unsigned int NumParticles() const;

    Parity GetParity() const;
    bool operator<(const Configuration& other) const;
    bool operator==(const Configuration& other) const;
    virtual std::string Name() const;

protected:
    /** Map single particle state to occupancy. */
    std::map<SingleParticleInfo, unsigned int> Config;
    mutable std::map<SingleParticleInfo, unsigned int>::const_iterator it;
};

typedef std::list<Configuration> ConfigList;

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

#endif