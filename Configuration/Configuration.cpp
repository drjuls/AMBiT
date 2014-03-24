#include "Include.h"
#include "Configuration.h"
#include "HartreeFock/NonRelInfo.h"
#include "Universal/MathConstant.h"
#include "HartreeFock/ConfigurationParser.h"

Configuration::Configuration(const std::string& name)
{
    *this = ConfigurationParser::ParseConfiguration(name);
}

int Configuration::GetOccupancy(const NonRelInfo& info) const
{
    auto r_it = config.find(info);

     if(r_it != config.end())
        return r_it->second;
    else
        return 0;
}

int& Configuration::operator[](const NonRelInfo& info)
{
    if(find(info) == end())
        config[info] = 0;

    return config[info];
}

Configuration::iterator Configuration::find(const NonRelInfo& info)
{
    return config.find(info);
}

Configuration::const_iterator Configuration::find(const NonRelInfo& info) const
{
    return config.find(info);
}

Configuration::iterator Configuration::erase(const_iterator position)
{
    if(position == config.end())
        return config.end();
    else
        return config.erase(position);
}

int Configuration::erase(const NonRelInfo& info)
{
    return config.erase(info);
}

Configuration::iterator Configuration::erase(const_iterator first, const_iterator last)
{
    return config.erase(first, last);
}

bool Configuration::RemoveSingleParticle(const NonRelInfo& info)
{
    auto r_it = config.find(info);
    if(r_it != config.end())
    {
        if(r_it->second <= 1)
            config.erase(r_it);
        else
            r_it->second = r_it->second - 1;

        return true;
    }
    else
        return false;
}

bool Configuration::AddSingleParticle(const NonRelInfo& info)
{
    iterator a_it = config.find(info);
    if(a_it != config.end())
    {
        if(a_it->second >= 4*(a_it->first.L()) + 2) // maximum number of electrons
        {   a_it->second = 4*(a_it->first.L()) + 2;
            return false;
        }
        else
            a_it->second = a_it->second + 1;
    }
    else
        config[info] = 1;

    return true;
}

int Configuration::ParticleNumber() const
{
    int num = 0;
    for(auto& value : config)
        num += value.second;

    return num;
}

int Configuration::ExcitationNumber() const
{
    int num = 0;
    for(auto& value : config)
        num += value.second;

    return num;
}

Parity Configuration::GetParity() const
{
    int sum = 0;
    for(auto value : config)
        sum += (value.first.L() * abs(value.second));

    if(sum%2 == 0)
        return even;
    else
        return odd;
}

std::string Configuration::Name(bool aSpaceFirst) const
{
    std::string name;
    if(aSpaceFirst)
    {
        auto m_it = config.begin();

        while(m_it != config.end())
        {
            name.append(" " + NonRelInfo(m_it->first).Name());
    
            char buffer[20];
            sprintf(buffer, "%d", m_it->second);
            name.append(buffer);
    
            m_it++;
        }
    }
    else
    {
        auto m_it = config.begin();

        bool IsFirstTerm = true;
        while(m_it != config.end())
        {
            if(IsFirstTerm)
            {
                name.append(NonRelInfo(m_it->first).Name());
                IsFirstTerm = false;
            }
            else
            {
                name.append(" " + NonRelInfo(m_it->first).Name());
            }

    
            char buffer[20];
            sprintf(buffer, "%d", m_it->second);
            name.append(buffer);
    
            m_it++;
        }
    }
    return name;
}

std::string Configuration::ShortName() const
{
    auto m_it = config.begin();
    std::string name;
    while(m_it != config.end())
    {
        name.append(NonRelInfo(m_it->first).Name());

        char buffer[20];
        sprintf(buffer, "%d", m_it->second);
        name.append(buffer);

        m_it++;
    }
    return name;
}

bool Configuration::operator<(const Configuration& other) const
{
    auto first = config.begin();
    auto second = other.config.begin();

    while((first != config.end()) && (second != other.config.end()))
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

    if((first == config.end()) && (second != other.config.end()))
        return true;
    else return false;
}

bool Configuration::operator==(const Configuration& other) const
{
    auto first = config.begin();
    auto second = other.config.begin();

    while((first != config.end()) && (second != other.config.end()))
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

    if((first != config.end()) || (second != other.config.end()))
        return false;
    else return true;
}

void Configuration::Write(FILE* fp) const
{
    // Write config
    unsigned int size = config.size();
    fwrite(&size, sizeof(unsigned int), 1, fp);

    auto cit = config.begin();

    while(cit != config.end())
    {
        unsigned int pqn = cit->first.PQN();
        int kappa = cit->first.Kappa();
        unsigned int occupancy = cit->second;
        
        // Write PQN, Kappa, occupancy
        fwrite(&pqn, sizeof(unsigned int), 1, fp);    
        fwrite(&kappa, sizeof(int), 1, fp);    
        fwrite(&occupancy, sizeof(unsigned int), 1, fp);    

        cit++;
    }
}

void Configuration::Read(FILE* fp)
{
    // Clear the current configuration
    config.clear();
    
    // Read config
    unsigned int size;
    fread(&size, sizeof(unsigned int), 1, fp);
    
    for(unsigned int i = 0; i < size; i++)
    {
        unsigned int pqn;
        int kappa;
        unsigned int occupancy;

        // Read PQN, Kappa, occupancy
        fread(&pqn, sizeof(unsigned int), 1, fp);    
        fread(&kappa, sizeof(int), 1, fp);    
        fread(&occupancy, sizeof(unsigned int), 1, fp);    

        config[OrbitalInfo(pqn, kappa)] = occupancy;
    }
}

void ConfigList::Print()
{
    ConfigList::iterator it = begin();
    while(it != end())
    {   *outstream << it->Name() << ",";
        it++;
    }
    *outstream << std::endl;
}

ConfigurationPair::ConfigurationPair(const Configuration& aConfiguration, const double& aDouble)
{
    first = aConfiguration;
    second = aDouble;

}

void ConfigurationSet::Print()
{
    ConfigurationSet::iterator cs_it;
    for(cs_it = begin(); cs_it != end(); cs_it++)
    {
        *outstream << std::setw(20) << cs_it->first.Name() << "  " << std::setprecision(2) << cs_it->second << "%" << std::endl;
    }
}
