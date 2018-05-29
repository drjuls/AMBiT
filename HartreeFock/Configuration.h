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

    It also sorts with hole states (negative occupancy) first, which then
    propagates to the sorting in Projections. This is required
    by ManyBodyOperator::GetProjectionDifferences().
 */
template <class OrbitalType, class OccupancyType>
class Configuration
{
protected:
    typedef Configuration<OrbitalType, OccupancyType> BaseConfiguration;
    typedef std::pair<OrbitalType, OccupancyType> PairType;

public:
    Configuration() = default;
    virtual ~Configuration() = default;

    typedef typename std::vector<std::pair<OrbitalType, OccupancyType>>::iterator iterator;
    typedef typename std::vector<std::pair<OrbitalType, OccupancyType>>::const_iterator const_iterator;

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

    /** Scale all occupancies by scale_factor. */
    template<typename U>
    typename boost::enable_if<boost::is_convertible<U, OccupancyType>, BaseConfiguration>::type
        operator*(U other) const;

    bool operator<(const BaseConfiguration& other) const;
    bool operator==(const BaseConfiguration& other) const
    {
        return m_config == other.m_config;
    }

    /** Get occupancy of a particular single particle state (zero if absent). */
    OccupancyType GetOccupancy(const OrbitalType& info) const;

    inline iterator begin() { return m_config.begin(); }
    inline const_iterator begin() const { return m_config.begin(); }
    inline iterator end() { return m_config.end(); }
    inline const_iterator end() const { return m_config.end(); }

    void clear() { m_config.clear(); num_hole_configs = 0; }
    bool empty() const { return m_config.empty(); }
    unsigned int size() const { return m_config.size(); }

    void erase(const OrbitalType& info)
    {
        insert(std::make_pair(info, 0));
    }

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
    {
        OccupancyType ret = OccupancyType(0);
        for(auto pair: m_config)
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
    {
        OccupancyType sum = 0;
        for(auto& pair: m_config)
            sum += pair.first.L() * abs(pair.second);
        if(sum%2 == 0)
            return Parity::even;
        else
            return Parity::odd;
    }

    virtual std::string Name(char sep = ' ') const;
    virtual std::string NameNoSpaces() const
    {   return BaseConfiguration::Name('_');
    }

    friend std::ostream& operator<<(std::ostream& stream, const BaseConfiguration& config) { return stream << config.Name(); }

    /** GetConfigDifferencesCount() is only defined for integral OccupancyTypes. */
    template<typename U = OccupancyType>
    typename boost::enable_if<boost::is_integral<U>, unsigned int>::type
        GetConfigDifferencesCount(const BaseConfiguration& other) const;

protected:
    std::vector<PairType> m_config;
    int num_hole_configs = {0};

    iterator hole_end() { return m_config.begin() + num_hole_configs; }
    const_iterator hole_end() const { return m_config.begin() + num_hole_configs; }

    iterator LocateInHoles(const OrbitalType& val)
    {
        return std::lower_bound(m_config.begin(), hole_end(), val, [](const PairType& pair, const OrbitalType& val)
        {
            return (pair.first < val);
        });
    }

    const_iterator LocateInHoles(const OrbitalType& val) const
    {
        return std::lower_bound(m_config.begin(), hole_end(), val, [](const PairType& pair, const OrbitalType& val)
        {
            return (pair.first < val);
        });
    }

    iterator LocateInElectrons(const OrbitalType& val)
    {
        return std::lower_bound(hole_end(), m_config.end(), val, [](const PairType& pair, const OrbitalType& val)
        {
            return (pair.first < val);
        });
    }

    const_iterator LocateInElectrons(const OrbitalType& val) const
    {
        return std::lower_bound(hole_end(), m_config.end(), val, [](const PairType& pair, const OrbitalType& val)
        {
            return (pair.first < val);
        });
    }
};

template <class OrbitalType, class OccupancyType>
bool Configuration<OrbitalType, OccupancyType>::operator<(const Configuration<OrbitalType, OccupancyType>& other) const
{
    auto first = m_config.begin();
    auto second = other.m_config.begin();

    while((first != m_config.end()) && (second != other.m_config.end()))
    {
        // Order by is_hole state first
        if(first->second * second->second < 0)  // If one is hole & one is particle
            return (first->second < 0);

        // Then order by single particle info
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
template <typename U>
typename boost::enable_if<boost::is_convertible<U, OccupancyType>, Configuration<OrbitalType, OccupancyType>>::type&
Configuration<OrbitalType, OccupancyType>::operator+=(const Configuration<OrbitalType, U>& other)
{
    for(auto& pair: other)
    {
        // Look for hole states and electron states separately
        auto hit = LocateInHoles(pair.first);

        // Found in hole states
        if(hit != hole_end() && hit->first == pair.first)
        {
            hit->second += pair.second;
            auto occ = hit->second;
            if(occ == 0)
            {   m_config.erase(hit);
                num_hole_configs--;
            }
            else if(occ > 0) // Move to particle states
            {
                m_config.erase(hit);
                num_hole_configs--;

                auto eit = LocateInElectrons(pair.first);
                m_config.insert(eit, std::make_pair(pair.first, occ));
            }

            continue;
        }

        // Found in electron states
        auto eit = LocateInElectrons(pair.first);
        if(eit != end() && eit->first == pair.first)
        {
            eit->second += pair.second;
            auto occ = eit->second;
            if(occ == 0)
                m_config.erase(eit);
            else if(occ < 0) // Move to hole states
            {
                m_config.erase(eit);
                auto hit = LocateInHoles(pair.first);
                m_config.insert(hit, std::make_pair(pair.first, occ));
                num_hole_configs++;
            }

            continue;
        }

        // Not found at all, so insert
        if(pair.second < 0)
        {
            m_config.insert(hit, pair);
            num_hole_configs++;
        }
        else if(pair.second > 0)
        {
            m_config.insert(eit, pair);
        }
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
    return ((*this) += other * OccupancyType(-1));
}

template <class OrbitalType, class OccupancyType>
template <typename U>
typename boost::enable_if<boost::is_convertible<U, OccupancyType>, Configuration<OrbitalType, OccupancyType>>::type
Configuration<OrbitalType, OccupancyType>::operator-(const Configuration<OrbitalType, U>& other) const
{
    BaseConfiguration ret = other * OccupancyType(-1);
    ret += (*this);
    return ret;
}

template <class OrbitalType, class OccupancyType>
template <typename U>
typename boost::enable_if<boost::is_convertible<U, OccupancyType>, Configuration<OrbitalType, OccupancyType>>::type
Configuration<OrbitalType, OccupancyType>::operator*(U scale_factor) const
{
    if(scale_factor == 0)
        return BaseConfiguration();

    BaseConfiguration ret(*this);

    for(auto& pair: ret.m_config)
    {
        pair.second *= scale_factor;
    }

    if(scale_factor < 0)
    {
        // Need to swap electron and hole states
        std::vector<PairType> new_config;
        new_config.insert(new_config.begin(), ret.m_config.begin()+num_hole_configs, ret.end());
        new_config.insert(new_config.end(), ret.m_config.begin(), ret.m_config.begin()+num_hole_configs);

        ret.m_config = new_config;
        ret.num_hole_configs = size() - num_hole_configs;
    }

    return ret;
}

template <class OrbitalType, class OccupancyType>
OccupancyType Configuration<OrbitalType, OccupancyType>::GetOccupancy(const OrbitalType& info) const
{
    auto it = LocateInHoles(info);
    if(it != hole_end() && it->first == info)
        return it->second;

    it = LocateInElectrons(info);
    if(it != end() && it->first == info)
        return it->second;

    return 0;
}

template <class OrbitalType, class OccupancyType>
typename Configuration<OrbitalType, OccupancyType>::iterator Configuration<OrbitalType, OccupancyType>::insert(const std::pair<OrbitalType, OccupancyType>& val)
{
    iterator hit = LocateInHoles(val.first);
    if(val.second < 0)
    {
        if(hit != hole_end() && hit->first == val.first)
            hit->second = val.second;
        else
        {   hit = m_config.insert(hit, val);
            num_hole_configs++;
        }
    }
    else if(hit != hole_end() && hit->first == val.first)
    {   m_config.erase(hit);
        num_hole_configs--;
        hit = m_config.end();
    }

    iterator eit = LocateInElectrons(val.first);
    if(val.second > 0)
    {
        if(eit != end() && eit->first == val.first)
            eit->second = val.second;
        else
            eit = m_config.insert(eit, val);
    }
    else if(eit != end() && eit->first == val.first)
    {
        m_config.erase(eit);
        eit = m_config.end();
    }

    if(val.second > 0)
        return eit;
    else
        return hit;
}

template <class OrbitalType, class OccupancyType>
bool Configuration<OrbitalType, OccupancyType>::RemoveSingleParticle(const OrbitalType& info)
{
    iterator it = LocateInElectrons(info);
    if(it == end() || it->first != info)
        it = LocateInHoles(info);

    if(it != end() && it->first == info)
    {
        if(it->second <= -info.MaxNumElectrons())
            return false;
        else if(it->second == 1)
            m_config.erase(it);
        else
            it->second--;
    }
    else // Not found, create hole
    {
        m_config.insert(it, std::make_pair(info, -1));
        num_hole_configs++;
    }

    return true;
}

template <class OrbitalType, class OccupancyType>
bool Configuration<OrbitalType, OccupancyType>::AddSingleParticle(const OrbitalType& info)
{
    iterator it = LocateInHoles(info);
    if(it == hole_end() || it->first != info)
        it = LocateInElectrons(info);

    if(it != end() && it->first == info)
    {
        if(it->second >= info.MaxNumElectrons())
            return false;
        else if(it->second == -1)
        {   m_config.erase(it);
            num_hole_configs--;
        }
        else
            it->second++;
    }
    else // Not found, create electron
    {
        m_config.insert(it, std::make_pair(info, 1));
    }

    return true;
}

template <class OrbitalType, class OccupancyType>
std::string Configuration<OrbitalType, OccupancyType>::Name(char sep) const
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
    {   buffer << sep << it->first.Name() << it->second;
        it++;
    }
    return buffer.str();
}

/** GetConfigDifferencesCount() is only defined for integral OccupancyTypes. */
template <class OrbitalType, class OccupancyType>
template <typename U>
typename boost::enable_if<boost::is_integral<U>, unsigned int>::type
Configuration<OrbitalType, OccupancyType>::GetConfigDifferencesCount(const BaseConfiguration& other) const
{
    unsigned int diff_count = abs(ElectronNumber() - other.ElectronNumber());
    auto it1 = begin();
    auto it2 = other.begin();

    auto pairless = [](const PairType& left, const PairType& right)
    {
        if(left.second * right.second < 0)
            return (left.second < 0);
        else
            return (left.first < right.first);
    };

    // Note: assumes first and second are sorted
    while(it1 != end() && it2 != other.end())
    {
        if(pairless(*it1, *it2))
        {
            diff_count += abs(it1->second);
            it1++;
        }
        else if(pairless(*it2, *it1))
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

#endif
