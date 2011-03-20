#include "Include.h"
#include "Configuration.h"
#include "HartreeFock/NonRelInfo.h"
#include "Universal/Constant.h"

Configuration::Configuration(const Configuration& other):
    Config(other.Config)
{
    if(other.AtEnd())
        it = Config.end();
    else
        it = Config.find(other.GetInfo());
}

Configuration::Configuration(const std::string& name)
{
    // Set num_states[L] to highest pqn
    unsigned int p = 0;
    const char* all_states = name.c_str();
    unsigned int pqn, L, occupancy;

    // Skip whitespace
    while(all_states[p] && isspace(all_states[p]))
        p++;

    while(all_states[p])
    {
        // Get pqn
        pqn = atoi(all_states + p);
        while(all_states[p] && isdigit(all_states[p]))
            p++;
        // Get L
        if(all_states[p])
        {
            // Get L
            for(L = 0; L < 10; L++)
            {   if(Constant::SpectroscopicNotation[L] == all_states[p])
                    break;
            }
            p++;
        }
        // Get occupancy
        occupancy = atoi(all_states + p);
        while(all_states[p] && isdigit(all_states[p]))
            p++;

        // Check everything is okay and skip any whitespace
        if((L >= 10) || (pqn < L+1) || (occupancy > 4*L + 2))
        {   *errstream << "Configuration() initialised with problematic string " << name << std::endl;
            exit(1);
        }
        SetOccupancy(NonRelInfo(pqn, L), occupancy);
        while(all_states[p] && isspace(all_states[p]))
            p++;
    }

    First();
}

void Configuration::SetIterator(const OrbitalInfo& info) const
{
    it = Config.find(info);
}

OrbitalInfo Configuration::GetInfo() const
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
    {   OrbitalInfo info = GetInfo();
        if(occupancy)
            Config[info] = occupancy;
        else
        {   Next();
            std::map<OrbitalInfo, unsigned int>::iterator kill_it = Config.find(info);
            Config.erase(kill_it);
        }
    }
}

unsigned int Configuration::GetOccupancy(const OrbitalInfo& info) const
{
    std::map<OrbitalInfo, unsigned int>::const_iterator r_it = Config.find(info);

     if(r_it != Config.end())
        return r_it->second;
    else
        return 0;
}

bool Configuration::RemoveSingleParticle(const OrbitalInfo& info)
{
    std::map<OrbitalInfo, unsigned int>::iterator r_it = Config.find(info);
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

bool Configuration::AddSingleParticle(const OrbitalInfo& info)
{
    std::map<OrbitalInfo, unsigned int>::iterator a_it = Config.find(info);
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

bool Configuration::SetOccupancy(const OrbitalInfo& info, unsigned int occupancy)
{
    if(occupancy)
    {   if(occupancy <= 4 * info.L() + 2)
            Config[info] = occupancy;
        else
            return false;
    }
    else
    {   std::map<OrbitalInfo, unsigned int>::iterator kill_it = Config.find(info);
        if(kill_it != Config.end())
            Config.erase(kill_it);
    }

    return true;
}

unsigned int Configuration::NumParticles() const
{
    unsigned int num = 0;
    std::map<OrbitalInfo, unsigned int>::const_iterator m_it = Config.begin();
    while(m_it != Config.end())
    {   num += m_it->second;
        m_it++;
    }
    return num;
}

Parity Configuration::GetParity() const
{
    std::map<OrbitalInfo, unsigned int>::const_iterator m_it = Config.begin();
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
    std::map<OrbitalInfo, unsigned int>::const_iterator m_it = Config.begin();
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
    std::map<OrbitalInfo, unsigned int>::const_iterator first = Config.begin();
    std::map<OrbitalInfo, unsigned int>::const_iterator second = other.Config.begin();

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
    std::map<OrbitalInfo, unsigned int>::const_iterator first = Config.begin();
    std::map<OrbitalInfo, unsigned int>::const_iterator second = other.Config.begin();

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

void Configuration::Write(FILE* fp) const
{
    // Write config
    unsigned int size = Config.size();
    fwrite(&size, sizeof(unsigned int), 1, fp);

    std::map<OrbitalInfo, unsigned int>::const_iterator cit;
    cit = Config.begin();

    while(cit != Config.end())
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
    Config.clear();
    
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

        SetOccupancy(OrbitalInfo(pqn, kappa), occupancy);
    }
    
    First();
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
