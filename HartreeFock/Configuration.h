#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "HartreeFock/OrbitalInfo.h"
#include "HartreeFock/NonRelInfo.h"
#include "Universal/Enums.h"
#include <map>
#include <string>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

/** Configuration is similar to std::map<OrbitalType, OccupancyType>,
    with a simplified, and sometimes intentionally modified, interface.
    OrbitalType should have functions L() and Name().
 */
template <class OrbitalType, class OccupancyType>
class Configuration
{
protected:
    typedef Configuration<OrbitalType, OccupancyType> BaseConfiguration;
    typedef std::map<OrbitalType, OccupancyType> ConfigurationMapType;

public:
    Configuration() {}
    Configuration(const Configuration<OrbitalType, OccupancyType>& other): m_config(other.m_config) {}
    Configuration(Configuration<OrbitalType, OccupancyType>&& other): m_config(other.m_config) {}
    virtual ~Configuration() {}

    typedef typename std::map<OrbitalType, OccupancyType>::iterator iterator;
    typedef typename std::map<OrbitalType, OccupancyType>::const_iterator const_iterator;

    const BaseConfiguration& operator=(const BaseConfiguration& other) { m_config = other.m_config; return *this; }
    BaseConfiguration& operator=(BaseConfiguration&& other) { m_config.swap(other.m_config); return *this; }

    bool operator<(const BaseConfiguration& other) const;
    bool operator==(const BaseConfiguration& other) const;

    /** Get occupancy of a particular single particle state (zero if absent). */
    OccupancyType GetOccupancy(const OrbitalType& info) const
    {   auto r_it = m_config.find(info);
        if(r_it != m_config.end())
            return r_it->second;
        else
            return OccupancyType(0);
    }

    OccupancyType& operator[](const OrbitalType& info)
    {   auto pair = m_config.insert(std::make_pair(info, OccupancyType(0)));
        return pair.first->second;
    }

    iterator begin() { return m_config.begin(); }
    const_iterator begin() const { return m_config.begin(); }
    iterator end() { return m_config.end(); }
    const_iterator end() const { return m_config.end(); }
    iterator find(const OrbitalType& info) { return m_config.find(info); }
    const_iterator find(const OrbitalType& info) const { return m_config.find(info); }

    void clear() { m_config.clear(); }
    int count(const OrbitalType& info) { return m_config.count(info); }
    bool empty() const { return m_config.empty(); }
    unsigned int size() const { return m_config.size(); }

    iterator erase(const_iterator position)
    {   if(position == m_config.end())
            return m_config.end();
        else
            return m_config.erase(position);
    }
    int erase(const OrbitalType& info) { return m_config.erase(info); }
    iterator erase(const_iterator first, const_iterator last) { return m_config.erase(first, last); }

    /** Different to STL map insert, in that it always changes the occupancy associated
        with OrbitalType, regardless of whether it existed previously.
        Returns iterator to new value, or end() if occupancy = 0.
     */
    iterator insert(const std::pair<OrbitalType, OccupancyType>& val);

    /** Electron number = number of electrons - number of holes. */
    OccupancyType ElectronNumber() const
    {   OccupancyType ret = OccupancyType(0);
        for(auto& pair: m_config)
            ret += pair.second;
        return ret;
    }

    /** Particle number = number of electrons + number of holes. */
    template<typename U = OccupancyType>
    typename boost::enable_if<boost::is_integral<U>, OccupancyType>::type
    ParticleNumber() const
    {   OccupancyType ret = OccupancyType(0);
        for(auto& pair: m_config)
            ret += abs(pair.second);
        return ret;
    }
    /** Particle number = number of electrons + number of holes. */
    template<typename U = OccupancyType>
    typename boost::lazy_enable_if<boost::is_floating_point<U>, OccupancyType>::type
    ParticleNumber() const
    {   OccupancyType ret = OccupancyType(0);
        for(auto& pair: m_config)
            ret += fabs(pair.second);
        return ret;
    }

    /** GetParity() is only defined for integral OccupancyTypes. */
    template<typename U = OccupancyType>
    typename boost::enable_if<boost::is_integral<U>, Parity>::type
    GetParity() const
    {   OccupancyType sum = 0;
        for(auto& pair: m_config)
            sum += pair.first.L() * abs(pair.second);
        if(sum%2 == 0)
            return Parity::even;
        else
            return Parity::odd;
    }

    virtual std::string Name() const;
    friend std::ostream& operator<<(std::ostream& stream, const BaseConfiguration& config) { return stream << config.Name(); }

protected:
    std::map<OrbitalType, OccupancyType> m_config;
};

template <class OrbitalType, class OccupancyType>
bool Configuration<OrbitalType, OccupancyType>::operator<(const Configuration<OrbitalType, OccupancyType>& other) const
{   auto first = m_config.begin();
    auto second = other.m_config.begin();

    while((first != m_config.end()) && (second != other.m_config.end()))
    {
        // Order by single particle info
        if(first->first < second->first)
            return true;
        else if(second->first < first->first)
            return false;

        // Then by number of electrons
        if(first->second < second->second)
            return true;
        else if(second->second < first->second)
            return false;

        first++;
        second++;
    }

    if((first == m_config.end()) && (second != other.m_config.end()))
        return true;
    else return false;
}

template <class OrbitalType, class OccupancyType>
bool Configuration<OrbitalType, OccupancyType>::operator==(const Configuration<OrbitalType, OccupancyType>& other) const
{
    if(m_config != other.m_config)
        return false;

    if(m_config.size() != other.m_config.size())
        return false;

    auto first = m_config.begin();
    auto second = other.m_config.begin();

    while((first != m_config.end()) && (second != other.m_config.end()))
    {
        if(first->first != second->first)
            return false;

        // Then by number of electrons
        if(first->second != second->second)
            return false;

        first++;
        second++;
    }

    if((first != m_config.end()) || (second != other.m_config.end()))
        return false;
    else return true;
}

template <class OrbitalType, class OccupancyType>
typename Configuration<OrbitalType, OccupancyType>::iterator Configuration<OrbitalType, OccupancyType>::insert(const std::pair<OrbitalType, OccupancyType>& val)
{
    auto it = find(val.first);
    if(it != end())
    {
        if(val.second)
        {   it->second = val.second;
            return it;
        }
        else
        {   m_config.erase(it);
            return end();
        }
    }
    // Not found already
    if(val.second)
        return m_config.insert(val).first;
    else
        return end();
}

template <class OrbitalType, class OccupancyType>
std::string Configuration<OrbitalType, OccupancyType>::Name() const
{
    std::stringstream buffer;
    auto it = m_config.begin();
    if(m_config.size())
    {   buffer << it->first.Name() << it->second;
        it++;
    }
    while(it != m_config.end())
    {   buffer << " " << it->first.Name() << it->second;
        it++;
    }
    return buffer.str();
}

#endif
