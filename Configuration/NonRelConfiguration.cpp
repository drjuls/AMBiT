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
    {
        int prev_occ = GetOccupancy(element.first);
        insert(std::make_pair(NonRelInfo(element.first), prev_occ+element.second));
    }
}

std::string NonRelConfiguration::Name(bool aSpaceFirst) const
{
    std::string name;
    if(aSpaceFirst)
        name = " ";

    name.append(BaseConfiguration::Name());
    return name;
}

void NonRelConfiguration::Write(FILE* fp) const
{
    // Write config
    unsigned int config_size = m_config.size();
    file_err_handler->fwrite(&config_size, sizeof(unsigned int), 1, fp);

    for(auto& pair: *this)
    {
        int pqn = pair.first.PQN();
        int kappa = pair.first.Kappa();
        int occupancy = pair.second;

        // Write PQN, Kappa, occupancy
        file_err_handler->fwrite(&pqn, sizeof(int), 1, fp);
        file_err_handler->fwrite(&kappa, sizeof(int), 1, fp);
        file_err_handler->fwrite(&occupancy, sizeof(int), 1, fp);
    }
}

void NonRelConfiguration::Read(FILE* fp)
{
    // Clear the current configuration
    clear();

    unsigned int config_size;
    fread(&config_size, sizeof(unsigned int), 1, fp);

    for(unsigned int i = 0; i < config_size; i++)
    {
        int pqn;
        int kappa;
        int occupancy;

        // Read PQN, Kappa, occupancy
        fread(&pqn, sizeof(int), 1, fp);
        fread(&kappa, sizeof(int), 1, fp);
        fread(&occupancy, sizeof(int), 1, fp);

        insert(std::make_pair(NonRelInfo(pqn, kappa), occupancy));
    }
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

double NonRelConfiguration::CalculateConfigurationAverageEnergy(pOrbitalMapConst orbitals, pHFIntegrals one_body, pSlaterIntegrals two_body)
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
    // Sort 0 -> Nsmall
    auto itsmall = std::next(config_list->first.begin(), config_list->second);
    std::sort(config_list->first.begin(), itsmall);
    auto last = std::unique(config_list->first.begin(), itsmall);
    itsmall = config_list->first.erase(last, itsmall);
    config_list->second = itsmall - config_list->first.begin();

    // Sort Nsmall -> config_list.end()
    std::sort(itsmall, config_list->first.end());
    last = std::unique(itsmall, config_list->first.end());
    config_list->first.resize(last-config_list->first.begin());
}
