#include "Include.h"
#include "Configuration.h"
#include "NonRelInfo.h"

Configuration::Configuration(const Configuration& other):
    Config(other.Config)
{
    if(other.AtEnd())
        it = Config.end();
    else
        it = Config.find(other.GetInfo());
}

void Configuration::SetIterator(const SingleParticleInfo& info) const
{
    it = Config.find(info);
}

SingleParticleInfo Configuration::GetInfo() const
{
    return it->first;
}

unsigned int Configuration::GetOccupancy() const
{
    return it->second;
}

void Configuration::SetOccupancy(unsigned int occupancy)
{
    if(!AtEnd())
    {   SingleParticleInfo info = GetInfo();
        if(occupancy)
            Config[info] = occupancy;
        else
        {   Next();
            std::map<SingleParticleInfo, unsigned int>::iterator kill_it = Config.find(info);
            Config.erase(kill_it);
        }
    }
}

bool Configuration::RemoveSingleParticle(const SingleParticleInfo& info)
{
    std::map<SingleParticleInfo, unsigned int>::iterator r_it = Config.find(info);
    if(r_it != Config.end())
    {
        if(r_it->second <= 1)
            Config.erase(r_it);
        else
            r_it->second = r_it->second - 1;

        return true;
    }
    else
        return false;
}

bool Configuration::AddSingleParticle(const SingleParticleInfo& info)
{
    std::map<SingleParticleInfo, unsigned int>::iterator a_it = Config.find(info);
    if(a_it != Config.end())
    {
        if(a_it->second >= 4*(a_it->first.L()) + 2) // maximum number of electrons
        {   a_it->second = 4*(a_it->first.L()) + 2;
            return false;
        }
        else
            a_it->second = a_it->second + 1;
    }
    else
        Config[info] = 1;

    return true;
}

bool Configuration::SetOccupancy(const SingleParticleInfo& info, unsigned int occupancy)
{
    if(occupancy)
    {   if(occupancy <= 4 * info.L() + 2)
            Config[info] = occupancy;
        else
            return false;
    }
    else
    {   std::map<SingleParticleInfo, unsigned int>::iterator kill_it = Config.find(info);
        if(kill_it != Config.end())
            Config.erase(kill_it);
    }

    return true;
}

unsigned int Configuration::NumParticles() const
{
    unsigned int num = 0;
    std::map<SingleParticleInfo, unsigned int>::const_iterator m_it = Config.begin();
    while(m_it != Config.end())
    {   num += m_it->second;
        m_it++;
    }
    return num;
}

Parity Configuration::GetParity() const
{
    std::map<SingleParticleInfo, unsigned int>::const_iterator m_it = Config.begin();
    unsigned int sum = 0;
    while(m_it != Config.end())
    {
        sum += (m_it->first.L() * m_it->second);
        m_it++;
    }

    if(sum%2 == 0)
        return even;
    else
        return odd;
}

std::string Configuration::Name() const
{
    std::map<SingleParticleInfo, unsigned int>::const_iterator m_it = Config.begin();
    std::string name;
    while(m_it != Config.end())
    {
        name.append(" " + NonRelInfo(m_it->first).Name());

        char buffer[20];
        sprintf(buffer, "%d", m_it->second);
        name.append(buffer);

        m_it++;
    }
    return name;
}

bool Configuration::operator<(const Configuration& other) const
{
    std::map<SingleParticleInfo, unsigned int>::const_iterator first = Config.begin();
    std::map<SingleParticleInfo, unsigned int>::const_iterator second = other.Config.begin();

    while((first != Config.end()) && (second != other.Config.end()))
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

    if((first == Config.end()) && (second != other.Config.end()))
        return true;
    else return false;
}

bool Configuration::operator==(const Configuration& other) const
{
    std::map<SingleParticleInfo, unsigned int>::const_iterator first = Config.begin();
    std::map<SingleParticleInfo, unsigned int>::const_iterator second = other.Config.begin();

    while((first != Config.end()) && (second != other.Config.end()))
    {
        // Order by single particle info
        if(first->first != second->first)
            return false;
        
        // Then by number of electrons
        if(first->second != second->second)
            return false;

        first++;
        second++;
    }

    if((first != Config.end()) || (second != other.Config.end()))
        return false;
    else return true;
}
