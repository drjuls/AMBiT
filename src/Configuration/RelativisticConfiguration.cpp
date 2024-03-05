#include "Include.h"
#include "RelativisticConfiguration.h"
#include "ConfigGenerator.h"
#include "Universal/Eigensolver.h"
#include "HartreeFock/ConfigurationParser.h"

namespace Ambit
{
RelativisticConfiguration::RelativisticConfiguration(const std::string& name)
{
    *this = RelativisticConfiguration(ConfigurationParser::ParseConfiguration<OrbitalInfo, int>(name));
}

bool RelativisticConfiguration::GetProjections(pAngularDataLibrary data, const Symmetry& sym, int two_m)
{
    angular_data = data->GetData(*this, sym, two_m);

    if(angular_data == nullptr)
        return false;

    projections.clear();
    projections.reserve(angular_data->projection_size());
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

int RelativisticConfiguration::GetNumberOfLevels() const
{
    int result = 1;
    MathConstant* math = MathConstant::Instance();
    for(const auto& pair: *this)
    {
        result *= math->nChoosek(pair.first.MaxNumElectrons(), abs(pair.second));
    }
    return result;
}

double RelativisticConfiguration::CalculateConfigurationAverageEnergy(pOrbitalMapConst orbitals, pHFIntegrals one_body, pSlaterIntegrals two_body) const
{
    return Ambit::CalculateConfigurationAverageEnergy(*this, orbitals, one_body, two_body);
}

void RelativisticConfiguration::Write(FILE* fp) const
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

void RelativisticConfiguration::Read(FILE* fp)
{
    // Clear the current configuration
    clear();

    unsigned int config_size;
    file_err_handler->fread(&config_size, sizeof(unsigned int), 1, fp);

    for(unsigned int i = 0; i < config_size; i++)
    {
        int pqn;
        int kappa;
        int occupancy;

        // Read PQN, Kappa, occupancy
        file_err_handler->fread(&pqn, sizeof(int), 1, fp);
        file_err_handler->fread(&kappa, sizeof(int), 1, fp);
        file_err_handler->fread(&occupancy, sizeof(int), 1, fp);

        insert(std::make_pair(OrbitalInfo(pqn, kappa), occupancy));
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
}
