#include "LevelMap.h"
#include "Include.h"
#include <thread>

LevelMap::LevelMap(pAngularDataLibrary lib): angular_library(lib) {}

LevelMap::LevelMap(const std::string& file_id, pAngularDataLibrary lib):
    filename_prefix(file_id), angular_library(lib)
{}

LevelMap::LevelMap(pHamiltonianIDConst hamiltonian_example, const std::string& file_id, pAngularDataLibrary lib):
    filename_prefix(file_id), angular_library(lib)
{
    ReadLevelMap(hamiltonian_example);
}

LevelVector LevelMap::GetLevels(pHamiltonianID key)
{
    auto it = m_map.find(key);
    if(it == m_map.end())
        return LevelVector(key);
    else
        return it->second;
}

void LevelMap::Store(pHamiltonianID key, const LevelVector& level_vector)
{
    keys.insert(key);
    m_map.insert(std::make_pair(key, level_vector));

    if(ProcessorRank == 0 && level_vector.levels.size() && !filename_prefix.empty())
    {
        std::string filename = filename_prefix + ".levels";
        FILE* fp = fopen(filename.c_str(), "wb");
        
        // Write number of entries in LevelMap
        unsigned int num_hamiltonians = 0;
        for(auto& pair: m_map)
        {
            // Write all HamiltonianIDs that have been calculated
            if((pair.second.configs != nullptr) && pair.second.levels.size())
                num_hamiltonians++;
        }
        fwrite(&num_hamiltonians, sizeof(unsigned int), 1, fp);
        
        // Write all HamiltonianIDs and Level information
        for(auto& pair: m_map)
        {
            if((pair.second.configs != nullptr) && pair.second.levels.size())
            {
                // Write HamiltonianID
                pair.first->Write(fp);

                // Write configs
                pair.second.configs->Write(fp);

                // Write all level information
                unsigned int num_levels = pair.second.levels.size();
                fwrite(&num_levels, sizeof(unsigned int), 1, fp);
                
                for(auto& pl: pair.second.levels)
                {
                    const std::vector<double>& eigenvector = pl->GetEigenvector();
                    unsigned int N = eigenvector.size();
                    double eigenvalue = pl->GetEnergy();
                    double gfactor = pl->GetgFactor();
                    
                    fwrite(&eigenvalue, sizeof(double), 1, fp);
                    fwrite(&N, sizeof(unsigned int), 1, fp);
                    fwrite(eigenvector.data(), sizeof(double), N, fp);
                    fwrite(&gfactor, sizeof(double), 1, fp);
                }
            }
        }
        
        fclose(fp);
    }
}

void LevelMap::ReadLevelMap(pHamiltonianIDConst hamiltonian_example)
{
    std::string filename = filename_prefix + ".levels";
    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
        return;

    // Read number of entries in LevelMap
    unsigned int num_hamiltonians;
    fread(&num_hamiltonians, sizeof(unsigned int), 1, fp);

    // Read all HamiltonianIDs and Level information
    for(unsigned int i = 0; i < num_hamiltonians; i++)
    {
        // Read key and recover angular data
        pHamiltonianID key = hamiltonian_example->Clone();
        key->Read(fp);
        keys.insert(key);

        // Read configs and recover angular data
        LevelVector& levelvec = m_map[key];
        if(levelvec.hID == nullptr)
            levelvec.hID = key;
        levelvec.configs = std::make_shared<RelativisticConfigList>();
        levelvec.configs->Read(fp);

        for(auto& relconfig: *levelvec.configs)
            relconfig.GetProjections(angular_library, key->GetSymmetry(), key->GetTwoJ());

        // These should be generated already, but in case they weren't saved...
        angular_library->GenerateCSFs();

        // Read all level information
        unsigned int num_levels;
        fread(&num_levels, sizeof(unsigned int), 1, fp);

        levelvec.levels.reserve(num_levels);

        double eigenvalue;
        unsigned int N;
        std::vector<double> eigenvector;
        double gfactor;

        for(unsigned int index = 0; index < num_levels; index++)
        {
            fread(&eigenvalue, sizeof(double), 1, fp);
            fread(&N, sizeof(unsigned int), 1, fp);
            eigenvector.resize(N);
            fread(eigenvector.data(), sizeof(double), N, fp);
            fread(&gfactor, sizeof(double), 1, fp);

            levelvec.levels.emplace_back(std::make_shared<Level>(eigenvalue, eigenvector, key, gfactor));
        }
    }

    fclose(fp);
}

FileSystemLevelStore::FileSystemLevelStore(const std::string& file_prefix, pAngularDataLibrary lib):
    FileSystemLevelStore(".", file_prefix, lib)
{}

FileSystemLevelStore::FileSystemLevelStore(const std::string& dir_name, const std::string& file_prefix, pAngularDataLibrary lib):
    directory(dir_name), filename_prefix(file_prefix), angular_library(lib)
{
    while(!boost::filesystem::exists(directory) || !boost::filesystem::is_directory(directory))
    {
        if(ProcessorRank == 0)
        {   // Attempt to create, else fail
            if(!boost::filesystem::create_directory(directory))
            {   *errstream << "FileSystemLevelStore cannot create directory " << directory << std::endl;
                exit(1);
            }
        }
        else
        {   // Wait a short time while master rank makes directory
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }

    directory = boost::filesystem::canonical(directory);
}

LevelVector FileSystemLevelStore::GetLevels(pHamiltonianID key)
{
    LevelVector levelvec(key);
    std::string filename = filename_prefix + "." + key->Name() + ".levels";
    filename = (directory / filename).string();

    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
        return levelvec;

    // Read key and recover angular data
    pHamiltonianID read_key = key->Clone();
    read_key->Read(fp);

    if(*key == *read_key)
    {
        // Read configs and recover angular data
        levelvec.configs = std::make_shared<RelativisticConfigList>();
        levelvec.configs->Read(fp);

        for(auto& relconfig: *levelvec.configs)
            relconfig.GetProjections(angular_library, key->GetSymmetry(), key->GetTwoJ());

        // These should be generated already, but in case they weren't saved...
        angular_library->GenerateCSFs();

        // Read all level information
        unsigned int num_levels;
        fread(&num_levels, sizeof(unsigned int), 1, fp);

        levelvec.levels.reserve(num_levels);

        double eigenvalue;
        unsigned int N;
        std::vector<double> eigenvector;
        double gfactor;

        for(unsigned int index = 0; index < num_levels; index++)
        {
            fread(&eigenvalue, sizeof(double), 1, fp);
            fread(&N, sizeof(unsigned int), 1, fp);
            eigenvector.resize(N);
            fread(eigenvector.data(), sizeof(double), N, fp);
            fread(&gfactor, sizeof(double), 1, fp);

            levelvec.levels.emplace_back(std::make_shared<Level>(eigenvalue, eigenvector, read_key, gfactor));
        }
    }

    fclose(fp);
    return levelvec;
}

void FileSystemLevelStore::Store(pHamiltonianID key, const LevelVector& level_vector)
{
    keys.insert(key);
    if(ProcessorRank != 0 || level_vector.configs == nullptr || level_vector.levels.size() == 0)
        return;

    // Append level_vector to file
    std::string filename = filename_prefix + "." + key->Name() + ".levels";
    filename = (directory / filename).string();

    FILE* fp = fopen(filename.c_str(), "wb");
    if(!fp)
    {   *errstream << "FileSystemLevelStore::Store() cannot open " << filename << " for writing." << std::endl;
        exit(1);
    }

    // Write HamiltonianID and Level information
    key->Write(fp);
    level_vector.configs->Write(fp);

    unsigned int num_levels = level_vector.levels.size();
    fwrite(&num_levels, sizeof(unsigned int), 1, fp);

    for(auto& pl: level_vector.levels)
    {
        const std::vector<double>& eigenvector = pl->GetEigenvector();
        unsigned int N = eigenvector.size();
        double eigenvalue = pl->GetEnergy();
        double gfactor = pl->GetgFactor();

        fwrite(&eigenvalue, sizeof(double), 1, fp);
        fwrite(&N, sizeof(unsigned int), 1, fp);
        fwrite(eigenvector.data(), sizeof(double), N, fp);
        fwrite(&gfactor, sizeof(double), 1, fp);
    }

    fclose(fp);
}

void LevelVector::Print(double min_percentage) const
{   Print(min_percentage, false, 0.0);
}

void LevelVector::Print(double min_percentage, double max_energy) const
{   Print(min_percentage, true, max_energy);
}

void LevelVector::Print(double min_percentage, bool use_max_energy, double max_energy) const
{
    if(levels.size() == 0 || ProcessorRank != 0)
        return;

    if(!configs || !configs->size())
        return;

    bool use_min_percentage = (0. < min_percentage && min_percentage <= 100.);

    // Build map of non-rel configurations to percentage contributions
    std::map<NonRelConfiguration, double> percentages;
    if(use_min_percentage)
        for(auto& rconfig: *configs)
        {   percentages[rconfig] = 0.0; // Auto instantiate NonRelConfiguration from RelativisticConfiguration.
        }

    auto solution_it = levels.begin();
    *outstream << "Solutions for " << hID << " (N = " << (*solution_it)->GetEigenvectorLength() << "):\n";

    unsigned int index = 0;
    while(solution_it != levels.end() &&
          (!use_max_energy || (*solution_it)->GetEnergy() < max_energy))
    {
        double energy = (*solution_it)->GetEnergy();
        const std::vector<double>& eigenvector = (*solution_it)->GetEigenvector();

        *outstream << index << ": " << std::setprecision(8) << energy << "    "
                   << std::setprecision(12) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << " /cm\n";

        if(use_min_percentage)
        {
            // Clear percentages
            for(auto& pair: percentages)
                pair.second = 0.0;

            const double* eigenvector_csf = eigenvector.data();
            for(auto& rconfig: *configs)
            {
                double contribution = 0.0;
                for(unsigned int j = 0; j < rconfig.NumCSFs(); j++)
                {   contribution += (*eigenvector_csf) * (*eigenvector_csf) * 100.;
                    eigenvector_csf++;
                }

                percentages[rconfig] += contribution;
            }

            // Print important configurations.
            for(auto& pair: percentages)
            {
                if(pair.second > min_percentage)
                    *outstream << std::setw(20) << pair.first.Name() << "  "
                               << std::setprecision(2) << pair.second << "%\n";
            }
        }

        if(hID->GetTwoJ() != 0)
        {
            const double& gfactor = (*solution_it)->GetgFactor();
            if(!std::isnan(gfactor))
                *outstream << "    g-factor = " << std::setprecision(5) << gfactor << "\n";
        }

        *outstream << std::endl;
        solution_it++;
        index++;
    }
}

void LevelVector::PrintInline(bool print_leading_configuration, bool print_gfactors, std::string separator) const
{   PrintInline(false, 0.0, print_leading_configuration, print_gfactors, separator);
}

void LevelVector::PrintInline(double max_energy, bool print_leading_configuration, bool print_gfactors, std::string separator) const
{   PrintInline(true, max_energy, print_leading_configuration, print_gfactors, separator);
}

void LevelVector::PrintInline(bool use_max_energy, double max_energy, bool print_leading_configuration, bool print_gfactors, std::string separator) const
{
    if(levels.size() == 0 || ProcessorRank != 0)
        return;

    if(!configs || !configs->size())
        return;

    // This will exist if configs exists
    double J = hID->GetJ();
    Parity P = hID->GetParity();

    // Build map of non-rel configurations to percentage contributions
    std::map<NonRelConfiguration, double> percentages;
    if(print_leading_configuration)
        for(auto& rconfig: *configs)
        {   percentages[rconfig] = 0.0; // Auto instantiate NonRelConfiguration from RelativisticConfiguration.
        };

    auto solution_it = levels.begin();
    unsigned int index = 0;
    while(solution_it != levels.end() &&
          (!use_max_energy || (*solution_it)->GetEnergy() < max_energy))
    {
        double energy = (*solution_it)->GetEnergy();
        const std::vector<double>& eigenvector = (*solution_it)->GetEigenvector();

        *outstream << J << separator << Sign(P) << separator << index << separator
                   << std::setprecision(12) << energy;

        if(print_gfactors)
        {
            if(hID->GetTwoJ() != 0)
            {
                const double& gfactor = (*solution_it)->GetgFactor();
                *outstream << separator << std::setprecision(5) << gfactor;
            }
            else
                *outstream << separator << 0.0;
        }

        if(print_leading_configuration)
        {
            // Clear percentages
            for(auto& pair: percentages)
                pair.second = 0.0;

            const double* eigenvector_csf = eigenvector.data();
            for(auto& rconfig: *configs)
            {
                double contribution = 0.0;
                for(unsigned int j = 0; j < rconfig.NumCSFs(); j++)
                {   contribution += (*eigenvector_csf) * (*eigenvector_csf) * 100.;
                    eigenvector_csf++;
                }

                percentages[rconfig] += contribution;
            }

            // Get largest contribution
            auto largest_contribution
                = std::max_element(percentages.begin(), percentages.end(),
                    [](const std::pair<NonRelConfiguration, double>& p1, const std::pair<NonRelConfiguration, double>& p2)
                    { return p1.second < p2.second; });

            *outstream << separator << largest_contribution->first.NameNoSpaces()
                       << separator << std::setprecision(2) << largest_contribution->second;
        }

        *outstream << separator << hID->Name() << std::endl;

        solution_it++;
        index++;
    }
}
