#include "Include.h"
#include "NonRelConfiguration.h"
#include "HartreeFock/NonRelInfo.h"
#include "Universal/MathConstant.h"
#include "HartreeFock/ConfigurationParser.h"

NonRelConfiguration::NonRelConfiguration(const std::string& name)
{
    *this = NonRelConfiguration(ConfigurationParser::ParseConfiguration<NonRelInfo, int>(name));
}

NonRelConfiguration::NonRelConfiguration(const RelativisticConfiguration& other)
{
    for(auto& element: other)
    {   int& occ = m_config[NonRelInfo(element.first)];
        occ += element.second;
        if(!occ)
            erase(NonRelInfo(element.first));
    }
}

bool NonRelConfiguration::RemoveSingleParticle(const NonRelInfo& info)
{
    auto r_it = m_config.find(info);
    if(r_it != m_config.end())
    {
        if(r_it->second <= -(4*info.L() + 2))
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

bool NonRelConfiguration::AddSingleParticle(const NonRelInfo& info)
{
    iterator a_it = m_config.find(info);
    if(a_it != m_config.end())
    {
        if(a_it->second >= 4*info.L() + 2) // maximum number of electrons
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

std::string NonRelConfiguration::Name(bool aSpaceFirst) const
{
    std::string name;
    if(aSpaceFirst)
        name = " ";

    name.append(BaseConfiguration::Name());
    return name;
}

std::string NonRelConfiguration::ShortName() const
{
    auto m_it = m_config.begin();
    std::string name;
    while(m_it != m_config.end())
    {
        name.append(NonRelInfo(m_it->first).Name());

        char buffer[20];
        sprintf(buffer, "%d", m_it->second);
        name.append(buffer);

        m_it++;
    }
    return name;
}

ConfigList::ConfigList(const RelativisticConfigList& rlist)
{
    // Generate non-relativistic configurations
    for(auto& rconfig: rlist)
    {
        m_list.push_back(NonRelConfiguration(rconfig));
    }

    m_list.sort(BaseComparator());
    unique();
}

std::ostream& operator<<(std::ostream& stream, const ConfigList& config_list)
{
    for(auto& config: config_list)
    {   stream << config.Name() << ",";
    }
    return stream;
}

ConfigurationPair::ConfigurationPair(const NonRelConfiguration& aConfiguration, const double& aDouble)
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
