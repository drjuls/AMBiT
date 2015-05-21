#include "LevelMap.h"

void Print(const LevelVector& levels, double min_percentage)
{   Print(levels, min_percentage, false, 0.0);
}

void Print(const LevelVector& levels, double min_percentage, double max_energy)
{   Print(levels, min_percentage, true, max_energy);
}

void WriteLevelMap(const LevelMap& level_map, const std::string& filename)
{
    FILE* fp = fopen(filename.c_str(), "wb");

    // Write number of entries in LevelMap
    unsigned int num_hamiltonians = 0;
    for(auto& pair: level_map)
    {
        // Write all HamiltonianIDs that have been calculated
        if((pair.first->GetRelativisticConfigList() != nullptr) && pair.second.size())
            num_hamiltonians++;
    }
    fwrite(&num_hamiltonians, sizeof(unsigned int), 1, fp);

    // Write all HamiltonianIDs and Level information
    for(auto& pair: level_map)
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

pLevelMap ReadLevelMap(pHamiltonianIDConst hamiltonian_example, const std::string& filename, pAngularDataLibrary angular_library)
{
    pLevelMap level_map = std::make_shared<LevelMap>();

    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
        return level_map;

    // Read number of entries in LevelMap
    unsigned int num_hamiltonians;
    fread(&num_hamiltonians, sizeof(unsigned int), 1, fp);

    // Read all HamiltonianIDs and Level information
    for(unsigned int i = 0; i < num_hamiltonians; i++)
    {
        // Read key and recover angular data
        pHamiltonianID key = hamiltonian_example->Clone();
        key->Read(fp);

        pRelativisticConfigList configs = key->GetRelativisticConfigList();
        for(auto& relconfig: *configs)
            relconfig.GetProjections(angular_library, key->GetSymmetry(), key->GetTwoJ());

        // These should be generated already, but in case they weren't saved...
        angular_library->GenerateCSFs();

        // Read all level information
        unsigned int num_levels;
        fread(&num_levels, sizeof(unsigned int), 1, fp);

        LevelVector& levels = (*level_map)[key];
        levels.reserve(num_levels);

        for(unsigned int index = 0; index < num_levels; index++)
        {
            double eigenvalue;
            unsigned int N;
            std::vector<double> eigenvector;
            double gfactor;

            fread(&eigenvalue, sizeof(double), 1, fp);
            fread(&N, sizeof(unsigned int), 1, fp);

            eigenvector.resize(N);
            fread(eigenvector.data(), sizeof(double), N, fp);
            fread(&gfactor, sizeof(double), 1, fp);

            levels.push_back(std::make_shared<Level>(eigenvalue, eigenvector, key, gfactor));
        }
    }

    fclose(fp);
    return level_map;
}

void Print(const LevelVector& levels, double min_percentage, bool use_max_energy, double max_energy)
{
    if(levels.size() == 0)
        return;

    pRelativisticConfigListConst configs = levels[0]->GetRelativisticConfigList();
    if(!configs || !configs->size())
        return;

    // This will exist if configs exists
    auto hID = levels[0]->GetHamiltonianID();

    bool print_gFactors = false;
    if(hID->GetTwoJ() != 0)
    {
        auto it = levels.begin();
        while(it != levels.end())
        {
            if((*it)->GetgFactor())
            {   print_gFactors = true;
                break;
            }
            it++;
        }
    }

    // Build map of non-rel configurations to percentage contributions
    std::map<NonRelConfiguration, double> percentages;
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

        if(print_gFactors)
            *outstream << "    g-factor = " << std::setprecision(5) << (*solution_it)->GetgFactor() << "\n";

        *outstream << std::endl;
        solution_it++;
        index++;
    }
}
