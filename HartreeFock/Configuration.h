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
    Configuration() = default;
    virtual ~Configuration() = default;

    typedef typename std::map<OrbitalType, OccupancyType>::iterator iterator;
    typedef typename std::map<OrbitalType, OccupancyType>::const_iterator const_iterator;

    /** Add occupancies of other configuration. */
    template<typename U>
    typename boost::enable_if<boost::is_convertible<U, OccupancyType>, BaseConfiguration>::type&
        operator+=(const Configuration<OrbitalType, U>& other);

    /** Add occupancies of two configurations. */
    template<typename U>
    typename boost::enable_if<boost::is_convertible<U, OccupancyType>, BaseConfiguration>::type
        operator+(const Configuration<OrbitalType, U>& other) const;

    /** Subtract occupancies of other configuration. */
    template<typename U>
    typename boost::enable_if<boost::is_convertible<U, OccupancyType>, BaseConfiguration>::type&
        operator-=(const Configuration<OrbitalType, U>& other);

    /** Subtract occupancies of two configurations. */
    template<typename U>
    typename boost::enable_if<boost::is_convertible<U, OccupancyType>, BaseConfiguration>::type
        operator-(const Configuration<OrbitalType, U>& other) const;

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

    /** Return success. */
    bool RemoveSingleParticle(const OrbitalType& info);
    /** Return success. */
    bool AddSingleParticle(const OrbitalType& info);

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

    /** GetConfigDifferencesCount() is only defined for integral OccupancyTypes. */
    template<typename U = OccupancyType>
    inline typename boost::enable_if<boost::is_integral<U>, unsigned int>::type
        GetConfigDifferencesCount(const BaseConfiguration& other) const
    {
        unsigned int diff_count = abs(ElectronNumber() - other.ElectronNumber());
        auto it1 = begin();
        auto it2 = other.begin();

        // Note: assumes first and second are sorted
        while(it1 != end() && it2 != other.end())
        {
            if(it1->first < it2->first)
            {
                diff_count += abs(it1->second);
                it1++;
            }
            else if(it2->first < it1->first)
            {
                diff_count += abs(it2->second);
                it2++;
            }
            else
            {   diff_count += abs(it1->second - it2->second);
                it1++; it2++;
            }
        }

        while(it1 != end())
        {   diff_count += abs(it1->second);
            it1++;
        }

        while(it2 != other.end())
        {   diff_count += abs(it2->second);
            it2++;
        }

        // Each promotion of electrons or holes makes two changes, so divide by two.
        return diff_count/2;
    }

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
template <typename U>
typename boost::enable_if<boost::is_convertible<U, OccupancyType>, Configuration<OrbitalType, OccupancyType>>::type&
Configuration<OrbitalType, OccupancyType>::operator+=(const Configuration<OrbitalType, U>& other)
{
    for(auto& pair: other)
    {
        OccupancyType& occupancy = m_config[pair.first];
        occupancy += pair.second;
        if(occupancy == 0)
            erase(pair.first);
    }

    return *this;
}

template <class OrbitalType, class OccupancyType>
template <typename U>
typename boost::enable_if<boost::is_convertible<U, OccupancyType>, Configuration<OrbitalType, OccupancyType>>::type
Configuration<OrbitalType, OccupancyType>::operator+(const Configuration<OrbitalType, U>& other) const
{
    BaseConfiguration ret(*this);
    ret += other;
    return ret;
}

template <class OrbitalType, class OccupancyType>
template <typename U>
typename boost::enable_if<boost::is_convertible<U, OccupancyType>, Configuration<OrbitalType, OccupancyType>>::type&
Configuration<OrbitalType, OccupancyType>::operator-=(const Configuration<OrbitalType, U>& other)
{
    for(auto& pair: other)
    {
        OccupancyType& occupancy = m_config[pair.first];
        occupancy -= pair.second;
        if(occupancy == 0)
            erase(pair.first);
    }

    return *this;
}

template <class OrbitalType, class OccupancyType>
template <typename U>
typename boost::enable_if<boost::is_convertible<U, OccupancyType>, Configuration<OrbitalType, OccupancyType>>::type
Configuration<OrbitalType, OccupancyType>::operator-(const Configuration<OrbitalType, U>& other) const
{
    BaseConfiguration ret(*this);
    ret -= other;
    return ret;
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
bool Configuration<OrbitalType, OccupancyType>::RemoveSingleParticle(const OrbitalType& info)
{
    auto r_it = m_config.find(info);
    if(r_it != m_config.end())
    {
        if(r_it->second <= -info.MaxNumElectrons())
            return false;
        else if(r_it->second == 1)
            m_config.erase(r_it);
        else
            r_it->second--;
    }
    else
        m_config[info] = -1;

    return true;
}

template <class OrbitalType, class OccupancyType>
bool Configuration<OrbitalType, OccupancyType>::AddSingleParticle(const OrbitalType& info)
{
    iterator a_it = m_config.find(info);
    if(a_it != m_config.end())
    {
        if(a_it->second >= info.MaxNumElectrons())
            return false;
        else if(a_it->second == -1)
            m_config.erase(a_it);
        else
            a_it->second++;
    }
    else
        m_config[info] = 1;

    return true;
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
    else // Vacuum
        buffer << "0";

    while(it != m_config.end())
    {   buffer << " " << it->first.Name() << it->second;
        it++;
    }
    return buffer.str();
}

#endif
