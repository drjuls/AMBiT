#include "LevelMap.h"

pHamiltonianIDComparator key_comparator;

unsigned int LevelStore::count(const pHamiltonianIDConst& key) const
{
    if(std::binary_search(m_lib.begin(), m_lib.end(), key, key_comparator))
        return 1;
    else
        return 0;
}

LevelStore::const_iterator LevelStore::find(const pHamiltonianIDConst& key) const
{
    auto it = std::lower_bound(m_lib.begin(), m_lib.end(), key, key_comparator);
    if(it != m_lib.end() && **it == *key)
        return it;
    else
        return m_lib.end();
}

LevelStore::const_iterator LevelStore::insert(pHamiltonianID key)
{
    auto insert_point = std::lower_bound(m_lib.begin(), m_lib.end(), key, key_comparator);

    if(insert_point != m_lib.end() && **insert_point == *key)
    {   // Key already exists in library. Copy RelativisticConfigList.
        pRelativisticConfigList configs = key->GetRelativisticConfigList();
        if(configs)
            (*insert_point)->SetRelativisticConfigList(configs);
        return insert_point;
    }
    else
    {   // Insert into library
        return m_lib.insert(insert_point, key);
    }
}

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
        return LevelVector();
    else
        return it->second;
}

void LevelMap::Store(pHamiltonianID key, const LevelVector& level_vector)
{
    insert(key);
    m_map[key] = level_vector;

    if(ProcessorRank == 0 && level_vector.size() && !filename_prefix.empty())
    {
        std::string filename = filename_prefix + ".levels";
        FILE* fp = fopen(filename.c_str(), "wb");
        
        // Write number of entries in LevelMap
        unsigned int num_hamiltonians = 0;
        for(auto& pair: m_map)
        {
            // Write all HamiltonianIDs that have been calculated
            if((pair.first->GetRelativisticConfigList() != nullptr) && pair.second.size())
                num_hamiltonians++;
        }
        fwrite(&num_hamiltonians, sizeof(unsigned int), 1, fp);
        
        // Write all HamiltonianIDs and Level information
        for(auto& pair: m_map)
        {
            if((pair.first->GetRelativisticConfigList() != nullptr) && pair.second.size())
            {
                pair.first->Write(fp);
                
                // Write all level information
                unsigned int num_levels = pair.second.size();
                fwrite(&num_levels, sizeof(unsigned int), 1, fp);
                
                for(auto& pl: pair.second)
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
        insert(key);

        pRelativisticConfigList configs = key->GetRelativisticConfigList();
        for(auto& relconfig: *configs)
            relconfig.GetProjections(angular_library, key->GetSymmetry(), key->GetTwoJ());

        // These should be generated already, but in case they weren't saved...
        angular_library->GenerateCSFs();

        // Read all level information
        unsigned int num_levels;
        fread(&num_levels, sizeof(unsigned int), 1, fp);

        LevelVector& levels = m_map[key];
        levels.reserve(num_levels);

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

            levels.push_back(std::make_shared<Level>(eigenvalue, eigenvector, key, gfactor));
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
    if(!boost::filesystem::exists(directory) || !boost::filesystem::is_directory(directory))
    {
        // Attempt to create, else fail
        if(!boost::filesystem::create_directory(directory))
        {   *errstream << "FileSystemLevelStore cannot create directory " << directory << std::endl;
            exit(1);
        }
    }

    directory = boost::filesystem::canonical(directory);
}

LevelVector FileSystemLevelStore::GetLevels(pHamiltonianID key)
{
    LevelVector levels;
    std::string filename = filename_prefix + "." + key->Name() + ".levels";
    filename = (directory / filename).string();

    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
        return levels;

    // Read key and recover angular data
    pHamiltonianID read_key = key->Clone();
    read_key->Read(fp);

    if(*key == *read_key)
    {
        pRelativisticConfigList configs = read_key->GetRelativisticConfigList();
        for(auto& relconfig: *configs)
            relconfig.GetProjections(angular_library, read_key->GetSymmetry(), read_key->GetTwoJ());

        // These should be generated already, but in case they weren't saved...
        angular_library->GenerateCSFs();

        // Read all level information
        unsigned int num_levels;
        fread(&num_levels, sizeof(unsigned int), 1, fp);

        levels.reserve(num_levels);

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

            levels.push_back(std::make_shared<Level>(eigenvalue, eigenvector, read_key, gfactor));
        }
    }

    fclose(fp);
    return levels;
}

void FileSystemLevelStore::Store(pHamiltonianID key, const LevelVector& level_vector)
{
    insert(key);
    if(ProcessorRank != 0 || key->GetRelativisticConfigList() == nullptr || level_vector.size() == 0)
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

    unsigned int num_levels = level_vector.size();
    fwrite(&num_levels, sizeof(unsigned int), 1, fp);

    for(auto& pl: level_vector)
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

void Print(const LevelVector& levels, double min_percentage)
{   Print(levels, min_percentage, false, 0.0);
}

void Print(const LevelVector& levels, double min_percentage, double max_energy)
{   Print(levels, min_percentage, true, max_energy);
}

void Print(const LevelVector& levels, double min_percentage, bool use_max_energy, double max_energy)
{
    if(levels.size() == 0 || ProcessorRank != 0)
        return;

    pRelativisticConfigListConst configs = levels[0]->GetRelativisticConfigList();
    if(!configs || !configs->size())
        return;

    // This will exist if configs exists
    auto hID = levels[0]->GetHamiltonianID();

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

void PrintInline(const LevelVector& levels, bool print_leading_configuration, bool print_gfactors, std::string separator)
{   PrintInline(levels, false, 0.0, print_leading_configuration, print_gfactors, separator);
}

void PrintInline(const LevelVector& levels, double max_energy, bool print_leading_configuration, bool print_gfactors, std::string separator)
{   PrintInline(levels, true, max_energy, print_leading_configuration, print_gfactors, separator);
}

void PrintInline(const LevelVector& levels, bool use_max_energy, double max_energy, bool print_leading_configuration, bool print_gfactors, std::string separator)
{
    if(levels.size() == 0 || ProcessorRank != 0)
        return;

    pRelativisticConfigListConst configs = levels[0]->GetRelativisticConfigList();
    if(!configs || !configs->size())
        return;

    // This will exist if configs exists
    auto hID = levels[0]->GetHamiltonianID();
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
