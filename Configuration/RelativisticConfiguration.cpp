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
    {   projections.push_back(Projection(*this, *it));
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

unsigned int RelativisticConfigList::NumCSFs() const
{
    unsigned int total = 0;
    for(auto& rconfig: m_list)
    {
        total += rconfig.NumCSFs();
    }

    return total;
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

/*
void RelativisticConfiguration::Write(FILE* fp) const
{
    // Write config
    unsigned int size = config.size();
    fwrite(&size, sizeof(unsigned int), 1, fp);    

    std::map<OrbitalInfo, int>::const_iterator cit;
    cit = config.begin();

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

    // Write number of projections. Projections will be generated again upon read,
    // so this serves as a sanity check.
    size = projections.size();
    fwrite(&size, sizeof(unsigned int), 1, fp);

    if(size)
    {   // Write twoM (all projections should have the same value).
        int twoM = projections.front().GetTwoM();
        fwrite(&twoM, sizeof(int), 1, fp);

        // Write JCoefficients. (These can be costly to generate in some cases.)
        fwrite(&num_states, sizeof(unsigned int), 1, fp);
        if(num_states)
            fwrite(j_coefficients, sizeof(double), num_states * projections.size(), fp);
    }
}

void RelativisticConfiguration::Read(FILE* fp)
{
    // Clear the current configuration
    config.clear();
    if(num_states)
        delete[] j_coefficients;
    
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

    // Generate projections and check that the same number were stored.
    fread(&size, sizeof(unsigned int), 1, fp);
    
    if(size)
    {   // Read twoM.
        int twoM;
        fread(&twoM, sizeof(int), 1, fp);

        GenerateProjections(twoM);
        if(projections.size() != size)
        {   *errstream << "RelativisticConfiguration::Read(): " << Name()
                       << " has incorrect number of projections stored." << std::endl;
            exit(1);
        }

        // Read JCoefficients.
        fread(&num_states, sizeof(unsigned int), 1, fp);
        if(num_states)
        {   j_coefficients = new double[num_states * projections.size()];
            fread(j_coefficients, sizeof(double), num_states * projections.size(), fp);
        }
    }
}
 */
