#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "HartreeFock/StateInfo.h"
#include <map>
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
    void SetIterator(const StateInfo& info) const;

    StateInfo GetInfo() const;
    unsigned int GetOccupancy() const;
    void SetOccupancy(unsigned int occupancy);

    // Get occupancy of a particular single particle state.
    // (zero if absent)
    unsigned int GetOccupancy(const StateInfo& info) const;

    // These return the success of the operation.
    virtual bool RemoveSingleParticle(const StateInfo& info);
    virtual bool AddSingleParticle(const StateInfo& info);
    virtual bool SetOccupancy(const StateInfo& info, unsigned int occupancy);

    void Clear();
    bool Empty() const;

    virtual unsigned int NumParticles() const;

    Parity GetParity() const;
    bool operator<(const Configuration& other) const;
    bool operator==(const Configuration& other) const;
    virtual std::string Name() const;

    // File storage (binary)
    virtual void Write(FILE* fp) const;
    virtual void Read(FILE* fp);

protected:
    /** Map single particle state to occupancy. */
    std::map<StateInfo, unsigned int> Config;
    mutable std::map<StateInfo, unsigned int>::const_iterator it;
};

class ConfigList : public std::list<Configuration> 
{
public:
    void print();
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

#endif
