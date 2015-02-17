#include "Include.h"
#include "RelativisticConfiguration.h"
#include "ConfigGenerator.h"
#include "Universal/Eigensolver.h"

RelativisticConfiguration::RelativisticConfiguration(const RelativisticConfiguration& other):
    BaseConfiguration(other), projections(other.projections), angular_data(other.angular_data)
{}

RelativisticConfiguration::RelativisticConfiguration(RelativisticConfiguration&& other):
    BaseConfiguration(other), projections(other.projections), angular_data(other.angular_data)
{}

const RelativisticConfiguration& RelativisticConfiguration::operator=(const RelativisticConfiguration& other)
{
    m_config = other.m_config;
    angular_data = other.angular_data;
    projections = other.projections;

    return *this;
}

RelativisticConfiguration& RelativisticConfiguration::operator=(RelativisticConfiguration&& other)
{
    m_config.swap(other.m_config);
    angular_data = other.angular_data;
    projections.swap(other.projections);

    return *this;
}

bool RelativisticConfiguration::GetProjections(pAngularDataLibrary data)
{
    angular_data = (*data)[*this];

    if(angular_data == nullptr)
        return false;

    auto it = angular_data->projection_begin();
    while(it != angular_data->projection_end())
    {   projections.emplace_back(*this, *it);
        it++;
    }

    return true;
}

unsigned int RelativisticConfiguration::NumCSFs() const
{
    if(!angular_data)
        return 0;
    else
        return angular_data->NumCSFs();
}

int RelativisticConfiguration::GetTwiceMaxProjection() const
{
    int maximum_two_m = 0;
    for(auto& it: m_config)
    {   // max M = J + (J-1) + (J-2) + ...
        //       = NJ - N(N-1)/2
        int N = abs(it.second);
        maximum_two_m += N * it.first.TwoJ() - N * (N-1);
    }

    return maximum_two_m;
}

void RelativisticConfiguration::Write(FILE* fp) const
{
    // Write config
    unsigned int config_size = m_config.size();
    fwrite(&config_size, sizeof(unsigned int), 1, fp);

    for(auto& pair: *this)
    {
        int pqn = pair.first.PQN();
        int kappa = pair.first.Kappa();
        int occupancy = pair.second;

        // Write PQN, Kappa, occupancy
        fwrite(&pqn, sizeof(int), 1, fp);
        fwrite(&kappa, sizeof(int), 1, fp);
        fwrite(&occupancy, sizeof(int), 1, fp);
    }
}

void RelativisticConfiguration::Read(FILE* fp)
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

        m_config[OrbitalInfo(pqn, kappa)] = occupancy;
    }
}

void RelativisticConfiguration::Print(bool include_CSFs) const
{
    bool print_CSFs = include_CSFs && angular_data;

    *outstream << Name();
    if(print_CSFs)
        *outstream << " (J = " << angular_data->GetTwoJ()/2. << ", M = " << angular_data->GetTwoM()/2. << ")";
    *outstream << ":\n";
    
    for(const auto& proj: projections)
    {
        *outstream << proj.Name() << " ";
    }
    *outstream << std::endl;

    if(print_CSFs)
    {
        const double* data = angular_data->GetCSFs();
        int num_CSFs = angular_data->NumCSFs();
        unsigned int num_projections = angular_data->projection_size();
    
        for(int i = 0; i < num_CSFs; i++)
        {   for(int j = 0; j < num_projections; j++)
                *outstream << data[j * num_CSFs + i] << " ";
            *outstream << std::endl;
        }
    }
}

unsigned int RelativisticConfigList::NumCSFs() const
{
    unsigned int total = 0;
    for(auto& rconfig: m_list)
    {
        total += rconfig.NumCSFs();
    }

    return total;
}

std::pair<RelativisticConfigList::iterator, int> RelativisticConfigList::operator[](unsigned int i)
{
    iterator it = begin();
    int csf = 0;
    for(unsigned int j = 0; j < i; j++)
    {   csf += it->NumCSFs();
        it++;
    }

    return std::pair<iterator, int>(it, csf);
}

std::pair<RelativisticConfigList::const_iterator, int> RelativisticConfigList::operator[](unsigned int i) const
{
    const_iterator it = begin();
    int csf = 0;
    for(unsigned int j = 0; j < i; j++)
    {   csf += it->NumCSFs();
        it++;
    }
    
    return std::pair<const_iterator, int>(it, csf);
}

RelativisticConfigList::const_projection_iterator RelativisticConfigList::projection_begin() const
{
    return const_projection_iterator(this);
}

RelativisticConfigList::const_projection_iterator RelativisticConfigList::projection_end() const
{
    const_iterator last = end();

    // Check size = 0 case
    if(begin() == last)
        return const_projection_iterator(this);

    last--;

    // NB: below we use const_projection_iterator(this, last->projection_end(), end(), 0) rather than
    //                  const_projection_iterator(this, last->projection_end(), end(), NumCSFs())
    //     even though the latter makes more sense. But the comparison is only made on
    //     const_projection_iterator.m_base, i.e. last->projection_end(), and so the CSF_index is meaningless
    return const_projection_iterator(this, last->projection_end(), end(), 0);
}

unsigned int RelativisticConfigList::projection_size() const
{
    unsigned int total = 0;
    for(auto& rconfig: m_list)
    {
        total += rconfig.projection_size();
    }

    return total;
}

void RelativisticConfigList::Read(FILE* fp)
{
    clear();
    unsigned int num_configs;
    fread(&num_configs, sizeof(unsigned int), 1, fp);

    for(unsigned int i = 0; i < num_configs; i++)
    {
        RelativisticConfiguration config;
        config.Read(fp);
        add(config);
    }
}

void RelativisticConfigList::Write(FILE* fp) const
{
    unsigned int num_configs = size();
    fwrite(&num_configs, sizeof(unsigned int), 1, fp);

    for(const auto& relconfig: *this)
    {
        relconfig.Write(fp);
    }
}
