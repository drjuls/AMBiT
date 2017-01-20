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

bool NonRelConfiguration::AddSingleParticle(const NonRelInfo& info)
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

std::string NonRelConfiguration::Name(bool aSpaceFirst) const
{
    std::string name;
    if(aSpaceFirst)
        name = " ";

    name.append(BaseConfiguration::Name());
    return name;
}

std::string NonRelConfiguration::NameNoSpaces() const
{
    std::string name = Name(false);
    std::replace(name.begin(), name.end(), ' ', '-');
    return name;
}

pRelativisticConfigList NonRelConfiguration::GenerateRelativisticConfigs()
{
    relconfiglist = std::make_shared<RelativisticConfigList>();
    RelativisticConfiguration rconfig;
    SplitNonRelInfo(begin(), rconfig);

    return relconfiglist;
}

void NonRelConfiguration::SplitNonRelInfo(NonRelConfiguration::const_iterator current_orbital, RelativisticConfiguration& relconfig)
{
    if(current_orbital == end())
    {   relconfiglist->push_back(relconfig);
        return;
    }

    const NonRelInfo& nrinfo(current_orbital->first);

    if(nrinfo.L() == 0)
    {   relconfig.insert(std::make_pair(current_orbital->first.GetFirstRelativisticInfo(), current_orbital->second));
        current_orbital++;
        SplitNonRelInfo(current_orbital, relconfig);
    }
    else
    {   // rinfo1 has kappa = -(L+1). rinfo2 has kappa = L.
        OrbitalInfo rinfo1 = nrinfo.GetFirstRelativisticInfo();
        OrbitalInfo rinfo2 = nrinfo.GetSecondRelativisticInfo();

        int num_electrons = current_orbital->second;
        int num_particles = abs(num_electrons);
        int start = mmax(0, num_particles - rinfo2.MaxNumElectrons());
        int end = mmin(num_particles, rinfo1.MaxNumElectrons());

        // Holes
        if(num_electrons < 0)
        {
            std::swap(start, end);
            start = -start;
            end = -end;
        }

        // Next orbital is the same for all loops
        current_orbital++;
        NonRelConfiguration::const_iterator next_orbital = current_orbital;

        for(int i=start; i<=end; i++)
        {
            RelativisticConfiguration new_rconfig(relconfig);

            if(i)
                new_rconfig.insert(std::make_pair(rinfo1, i));
            if(num_electrons - i)
                new_rconfig.insert(std::make_pair(rinfo2, num_electrons-i));

            SplitNonRelInfo(next_orbital, new_rconfig);
        }
    }
}

double NonRelConfiguration::CalculateConfigurationAverageEnergy(pOrbitalMapConst orbitals, pHFOperator one_body, pHartreeY two_body)
{
    if(std::isnan(config_average_energy))
    {
        double energy = 0.;
        num_levels = 0;

        if(relconfiglist == nullptr)
            GenerateRelativisticConfigs();

        for(const auto& rconfig: *relconfiglist)
        {
            int rconfig_num_levels = rconfig.GetNumberOfLevels();
            double rconfig_energy = rconfig.CalculateConfigurationAverageEnergy(orbitals, one_body, two_body);

            energy += rconfig_num_levels * rconfig_energy;
            num_levels += rconfig_num_levels;
        }

        config_average_energy = energy/num_levels;
    }

    return config_average_energy;
}

int NonRelConfiguration::GetTwiceMaxProjection() const
{
    int maximum_two_m = 0;

    // For each orbital, first fill largest values of M.
    // If odd number of electrons, last goes with spin up in last container.
    for(auto& it: m_config)
    {
        int N = abs(it.second);
        int L = it.first.L();
        if(N%2 == 0)
        {   // Two spins for each M cancel
            // max M = 2L + 2(L-1) + ...
            //       = 2[ nL - n(n-1)/2 ] where n = N/2
            maximum_two_m += N * (2 + 4*L - N)/2;
        }
        else
        {   // After (N-1) pairs (as above) add an odd electron
            // with M = L - (N-1)/2 and spin up.
            maximum_two_m += (N * (2 + 4*L - N) + 1)/2;
        }
    }

    return maximum_two_m;
}

int NonRelConfiguration::GetNumberOfLevels()
{
    if(!num_levels)
    {
        if(relconfiglist == nullptr)
            GenerateRelativisticConfigs();

        for(auto& rconfig: *relconfiglist)
            num_levels += rconfig.GetNumberOfLevels();
    }

    return num_levels;
}

std::ostream& operator<<(std::ostream& stream, const ConfigList& config_list)
{
    for(auto& config: config_list.first)
    {   stream << config.Name() << ",";
    }
    return stream;
}

void SortAndUnique(pConfigList& config_list)
{
    // Protect Nsmall
    config_list->second = mmin(config_list->second, config_list->first.size());

    // Sort 0 -> Nsmall
    std::list<NonRelConfiguration> smallside;
    auto itsmall = std::next(config_list->first.begin(), config_list->second);
    smallside.splice(smallside.begin(), config_list->first, config_list->first.begin(), itsmall);
    smallside.sort();
    smallside.unique();

    // Sort Nsmall -> config_list.end()
    config_list->first.sort();
    config_list->first.unique();
    config_list->second = smallside.size();
    config_list->first.splice(config_list->first.begin(), smallside);
}
