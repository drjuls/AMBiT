#include "LevelMap.h"
#include "Include.h"
#include <thread>

namespace Ambit
{
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
        FILE* fp = file_err_handler->fopen(filename.c_str(), "wb");
        
        // Write number of entries in LevelMap
        unsigned int num_hamiltonians = 0;
        for(auto& pair: m_map)
        {
            // Write all HamiltonianIDs that have been calculated
            if((pair.second.configs != nullptr) && pair.second.levels.size())
                num_hamiltonians++;
        }
        file_err_handler->fwrite(&num_hamiltonians, sizeof(unsigned int), 1, fp);
        
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
                file_err_handler->fwrite(&num_levels, sizeof(unsigned int), 1, fp);
                
                for(auto& pl: pair.second.levels)
                {
                    const std::vector<double>& eigenvector = pl->GetEigenvector();
                    unsigned int N = eigenvector.size();
                    double eigenvalue = pl->GetEnergy();
                    double gfactor = pl->GetgFactor();
                    
                    file_err_handler->fwrite(&eigenvalue, sizeof(double), 1, fp);
                    file_err_handler->fwrite(&N, sizeof(unsigned int), 1, fp);
                    file_err_handler->fwrite(eigenvector.data(), sizeof(double), N, fp);
                    file_err_handler->fwrite(&gfactor, sizeof(double), 1, fp);
                }
            }
        }
        
        file_err_handler->fclose(fp);
    }
}

void LevelMap::ReadLevelMap(pHamiltonianIDConst hamiltonian_example)
{
    std::string filename = filename_prefix + ".levels";
    FILE* fp = file_err_handler->fopen(filename.c_str(), "rb");
    if(!fp)
        return;

    // Read number of entries in LevelMap
    unsigned int num_hamiltonians;
    file_err_handler->fread(&num_hamiltonians, sizeof(unsigned int), 1, fp);

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
        file_err_handler->fread(&num_levels, sizeof(unsigned int), 1, fp);

        levelvec.levels.reserve(num_levels);

        double eigenvalue;
        unsigned int N;
        std::vector<double> eigenvector;
        double gfactor;
        gfactors_needed = false; // Assume the file has g-factors unless one or more of them is NaN

        for(unsigned int index = 0; index < num_levels; index++)
        {
            file_err_handler->fread(&eigenvalue, sizeof(double), 1, fp);
            file_err_handler->fread(&N, sizeof(unsigned int), 1, fp);
            eigenvector.resize(N);
            file_err_handler->fread(eigenvector.data(), sizeof(double), N, fp);
            file_err_handler->fread(&gfactor, sizeof(double), 1, fp);
            
            // Check if we need to re-calculate this g-factor
            if(std::isnan(gfactor))
                gfactors_needed = true;

            levelvec.levels.emplace_back(std::make_shared<Level>(eigenvalue, eigenvector, key, gfactor));
        }
    }

    file_err_handler->fclose(fp);
}

FileSystemLevelStore::FileSystemLevelStore(const std::string& file_prefix, pAngularDataLibrary lib):
    FileSystemLevelStore(".", file_prefix, lib)
{}

FileSystemLevelStore::FileSystemLevelStore(const std::string& dir_name, const std::string& file_prefix, pAngularDataLibrary lib):
    directory(dir_name), filename_prefix(file_prefix), angular_library(lib)
{
    while(!std::filesystem::exists(directory) || !std::filesystem::is_directory(directory))
    {
        if(ProcessorRank == 0)
        {   // Attempt to create, else fail
            if(!std::filesystem::create_directory(directory))
            {   *errstream << "FileSystemLevelStore cannot create directory " << directory << std::endl;
                exit(1);
            }
        }
        else
        {   // Wait a short time while master rank makes directory
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }

    directory = std::filesystem::canonical(directory);
}

LevelVector FileSystemLevelStore::GetLevels(pHamiltonianID key)
{
    LevelVector levelvec(key);
    std::string filename = filename_prefix + "." + key->Name() + ".levels";
    filename = (directory / filename).string();

    FILE* fp = file_err_handler->fopen(filename.c_str(), "rb");
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
        file_err_handler->fread(&num_levels, sizeof(unsigned int), 1, fp);

        levelvec.levels.reserve(num_levels);

        double eigenvalue;
        unsigned int N;
        std::vector<double> eigenvector;
        double gfactor;
        gfactors_needed = false; // Assume the file has g-factors unless one or more of them is NaN

        for(unsigned int index = 0; index < num_levels; index++)
        {
            file_err_handler->fread(&eigenvalue, sizeof(double), 1, fp);
            file_err_handler->fread(&N, sizeof(unsigned int), 1, fp);
            eigenvector.resize(N);
            file_err_handler->fread(eigenvector.data(), sizeof(double), N, fp);
            file_err_handler->fread(&gfactor, sizeof(double), 1, fp);

            // Check if we need to re-calculate this g-factor
            if(std::isnan(gfactor))
                gfactors_needed = true;

            levelvec.levels.emplace_back(std::make_shared<Level>(eigenvalue, eigenvector, read_key, gfactor));
        }
    }

    file_err_handler->fclose(fp);
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

    FILE* fp = file_err_handler->fopen(filename.c_str(), "wb");
    if(!fp)
    {   *errstream << "FileSystemLevelStore::Store() cannot open " << filename << " for writing." << std::endl;
        exit(1);
    }

    // Write HamiltonianID and Level information
    key->Write(fp);
    level_vector.configs->Write(fp);

    unsigned int num_levels = level_vector.levels.size();
    file_err_handler->fwrite(&num_levels, sizeof(unsigned int), 1, fp);

    for(auto& pl: level_vector.levels)
    {
        const std::vector<double>& eigenvector = pl->GetEigenvector();
        unsigned int N = eigenvector.size();
        double eigenvalue = pl->GetEnergy();
        double gfactor = pl->GetgFactor();

        file_err_handler->fwrite(&eigenvalue, sizeof(double), 1, fp);
        file_err_handler->fwrite(&N, sizeof(unsigned int), 1, fp);
        file_err_handler->fwrite(eigenvector.data(), sizeof(double), N, fp);
        file_err_handler->fwrite(&gfactor, sizeof(double), 1, fp);
    }

    file_err_handler->fclose(fp);
}
}
