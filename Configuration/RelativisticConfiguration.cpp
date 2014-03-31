#include "Include.h"
#include "RelativisticConfiguration.h"
#include "ConfigGenerator.h"
#include "Universal/Eigensolver.h"

RelativisticConfiguration::RelativisticConfiguration(const RelativisticConfiguration& other):
    config(other.config), projections(other.projections), angular_data(other.angular_data)
{}

int RelativisticConfiguration::ParticleNumber() const
{
    int num = 0;
    for(auto& value : config)
        num += value.second;

    return num;
}

/** Excitation number = number of electrons + number of holes. */
int RelativisticConfiguration::ExcitationNumber() const
{
    int num = 0;
    for(auto& value : config)
        num += value.second;

    return num;
}

/** Get occupancy of a particular single particle state (zero if absent). */
int RelativisticConfiguration::GetOccupancy(const OrbitalInfo& info) const
{
    auto it = find(info);
    if(it != end())
        return it->second;
    else
        return 0;
}

int& RelativisticConfiguration::operator[](const OrbitalInfo& info)
{
    return config[info];
}

RelativisticConfiguration::iterator RelativisticConfiguration::insert(const std::pair<OrbitalInfo, int>& val)
{
    auto it = find(val.first);
    if(it != end())
    {
        if(val.second)
        {   it->second = val.second;
            return it;
        }
        else
        {   config.erase(it);
            return end();
        }
    }

    // Not found already
    if(val.second)
    {   return config.insert(val).first;
    }
    else
        return end();
}

RelativisticConfiguration::iterator RelativisticConfiguration::find(const OrbitalInfo& info)
{
    for(auto it = config.begin(); it != config.end(); it++)
        if(it->first == info)
            return it;

    return config.end();
}

RelativisticConfiguration::const_iterator RelativisticConfiguration::find(const OrbitalInfo& info) const
{
    for(auto it = config.begin(); it != config.end(); it++)
        if(it->first == info)
            return it;

    return config.end();
}

void RelativisticConfiguration::clear()
{
    config.clear();
}

bool RelativisticConfiguration::empty() const
{
    return config.empty();
}

RelativisticConfiguration::iterator RelativisticConfiguration::erase(const_iterator position)
{
    return config.erase(position);
}

int RelativisticConfiguration::erase(const OrbitalInfo& info)
{
    auto it = find(info);
    if(it != end())
    {   erase(it);
        return 1;
    }
    else
        return 0;
}

RelativisticConfiguration::iterator RelativisticConfiguration::erase(const_iterator first, const_iterator last)
{
    return config.erase(first, last);
}

std::string RelativisticConfiguration::Name() const
{
    std::string ret;
    for(auto& element: config)
        ret += element.first.Name() + itoa(element.second);

    return ret;
}

bool RelativisticConfiguration::operator<(const RelativisticConfiguration& other) const
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
        if(abs(first->second) < abs(second->second))
            return true;
        else if(abs(first->second) > abs(second->second))
            return false;

        first++;
        second++;
    }

    if((first == config.end()) && (second != other.config.end()))
        return true;
    else return false;
}

bool RelativisticConfiguration::operator==(const RelativisticConfiguration& other) const
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
    for(auto& it: config)
    {   // max M = J + (J-1) + (J-2) + ...
        //       = NJ - N(N-1)/2
        int N = abs(it.second);
        maximum_two_m += N * it.first.TwoJ() - N * (N-1);
    }

    return maximum_two_m;
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
