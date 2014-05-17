#include "Level.h"
#include "Include.h"
#include "Universal/MathConstant.h"
#include "NonRelConfiguration.h"

bool LevelID::operator<(const LevelID& other) const
{
    if(m_twoJ != other.m_twoJ)
        return m_twoJ < other.m_twoJ;
    else if(m_parity != other.m_parity)
        return m_parity < other.m_parity;
    else
        return m_ID < other.m_ID;
}

Level::Level(const double& energy, const double* csf_eigenvector, pRelativisticConfigListConst configlist, unsigned int numCSFs):
    eigenvalue(energy), configs(configlist), N(numCSFs), gFactor(0.0)
{
    if(N == 0)
        N = configs->NumCSFs();
    eigenvector = new double[N];
    memcpy(eigenvector, csf_eigenvector, N * sizeof(double));
}

Level::Level(const Level& other):
    eigenvalue(other.eigenvalue), configs(other.configs), N(other.N), gFactor(other.gFactor)
{
    eigenvector = new double[N];
    memcpy(eigenvector, other.eigenvector, N * sizeof(double));
}

Level::Level(Level&& other):
    eigenvalue(other.eigenvalue), configs(other.configs), N(other.N), gFactor(other.gFactor)
{
    eigenvector = other.eigenvector;
    other.eigenvector = nullptr;
}

Level::~Level()
{
    if(eigenvector)
        delete[] eigenvector;
}

const Level& Level::operator=(const Level& other)
{
    if(eigenvector)
    {   delete[] eigenvector;
        eigenvector = nullptr;
    }

    configs = other.configs;
    eigenvalue = other.eigenvalue;
    N = other.N;
    gFactor = other.gFactor;

    if(other.eigenvector)
    {   eigenvector = new double[N];
        memcpy(eigenvector, other.eigenvector, N * sizeof(double));
    }

    return *this;
}

Level& Level::operator=(Level&& other)
{
    if(eigenvector)
    {   delete[] eigenvector;
        eigenvector = nullptr;
    }
    configs = other.configs;
    eigenvalue = other.eigenvalue;
    N = other.N;
    gFactor = other.gFactor;
    eigenvector = other.eigenvector;
    other.eigenvector = nullptr;

    return *this;
}

unsigned int LevelMap::size(const Symmetry& sym) const
{
    unsigned int total = 0;
    for_each(begin(), end(), [&](const Parent::value_type& pair){
        if(pair.first.GetSymmetry() == sym)
            total++;
    });

    return total;
}

std::set<Symmetry> LevelMap::GetSymmetries() const
{
    std::set<Symmetry> ret;

    for(auto pair: *this)
        ret.insert(pair.first.GetSymmetry());

    return ret;
}

void LevelMap::Print(const Symmetry& sym, double min_percentage, bool use_max_energy, double max_energy) const
{
    const_symmetry_iterator solution_it = begin(sym);
    if(solution_it == end(sym))
    {   *outstream << "No solutions found with J = " << sym.GetJ() << ", P = " << LowerName(sym.GetParity()) << std::endl;
        return;
    }

    *outstream << "Solutions for J = " << sym.GetJ() << ", P = " << LowerName(sym.GetParity()) << ":\n";

    bool print_gFactors = false;
    if(sym.GetJ() != 0)
    {
        while(solution_it != end(sym))
        {   if(solution_it->second->GetgFactor())
            {   print_gFactors = true;
                break;
            }
            solution_it++;
        }
    }

    // Build map of non-rel configurations to percentage contributions
    solution_it = begin(sym);
    std::map<NonRelConfiguration, double> percentages;
    for(auto& rconfig: *solution_it->second->GetRelativisticConfigList())
    {   percentages[rconfig] = 0.0; // Auto instantiate NonRelConfiguration from RelativisticConfiguration.
    }

    solution_it = begin(sym);
    while(solution_it != end(sym) &&
          (!use_max_energy || solution_it->second->GetEnergy() < max_energy))
    {
        double energy = solution_it->second->GetEnergy();
        const double* eigenvector = solution_it->second->GetEigenvector();

        *outstream << solution_it->first.GetID() << ": " << std::setprecision(8) << energy << "    "
                   << std::setprecision(12) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << " /cm\n";

        // Clear percentages
        for(auto& pair: percentages)
            pair.second = 0.0;

        const double* eigenvector_csf = eigenvector;
        for(auto& rconfig: *solution_it->second->GetRelativisticConfigList())
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
            *outstream << "    g-factor = " << std::setprecision(5) << solution_it->second->GetgFactor() << "\n";

        *outstream << std::endl;
        solution_it++;
    }
}

bool LevelMap::Read(const std::string& filename, const std::string& angular_directory)
{
    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
        return false;

    std::map<Symmetry, pRelativisticConfigList> sym_config_map;

    unsigned int num_symmetries;
    int two_j;
    Parity P;

    // Read symmetries and relativistic configurations
    fread(&num_symmetries, sizeof(unsigned int), 1, fp);

    for(unsigned int i = 0; i < num_symmetries; i++)
    {
        fread(&two_j, sizeof(int), 1, fp);
        fread(&P, sizeof(Parity), 1, fp);

        pRelativisticConfigList configs(new RelativisticConfigList());
        configs->Read(fp);
        sym_config_map[Symmetry(two_j, P)] = configs;

        // Recover angular data
        int particle_number = configs->front().ElectronNumber();
        pAngularDataLibrary angular_library(new AngularDataLibrary(particle_number, Symmetry(two_j, P), two_j, angular_directory));
        angular_library->Read();

        for(auto& relconfig: *configs)
            relconfig.GetProjections(angular_library);

        // These should be generated already, but in case they weren't saved...
        angular_library->GenerateCSFs();
    }

    // Read LevelIDs and Levels
    unsigned int num_levels;
    unsigned int ID;
    fread(&num_levels, sizeof(unsigned int), 1, fp);

    for(unsigned int i = 0; i < num_levels; i++)
    {
        // LevelID
        fread(&two_j, sizeof(int), 1, fp);
        fread(&P, sizeof(Parity), 1, fp);
        fread(&ID, sizeof(unsigned int), 1, fp);

        // Level
        pLevel level(new Level());
        fread(&level->eigenvalue, sizeof(double), 1, fp);
        fread(&level->N, sizeof(unsigned int), 1, fp);
        level->eigenvector = new double[level->N];
        fread(level->eigenvector, sizeof(double), level->N, fp);
        fread(&level->gFactor, sizeof(double), 1, fp);

        level->configs = sym_config_map[Symmetry(two_j, P)];

        // Add to this
        insert(std::make_pair(LevelID(two_j, P, ID), level));
    }

    fclose(fp);
    return true;
}

void LevelMap::Write(const std::string& filename) const
{
    FILE* fp = fopen(filename.c_str(), "wb");
    std::set<Symmetry> symmetries = GetSymmetries();

    // Write symmetries and relativistic configurations
    unsigned int num_symmetries = symmetries.size();
    fwrite(&num_symmetries, sizeof(unsigned int), 1, fp);

    for(const Symmetry& sym: symmetries)
    {
        int two_j = sym.GetTwoJ();
        Parity P = sym.GetParity();
        fwrite(&two_j, sizeof(int), 1, fp);
        fwrite(&P, sizeof(Parity), 1, fp);

        // Config list
        begin(sym)->second->GetRelativisticConfigList()->Write(fp);
    }

    // Write all LevelIDs and Level information
    unsigned int num_levels = size();
    fwrite(&num_levels, sizeof(unsigned int), 1, fp);

    for(auto& pair: *this)
    {
        // LevelID
        int two_j = pair.first.GetTwoJ();
        Parity P = pair.first.GetParity();
        unsigned int ID = pair.first.GetID();

        fwrite(&two_j, sizeof(int), 1, fp);
        fwrite(&P, sizeof(Parity), 1, fp);
        fwrite(&ID, sizeof(unsigned int), 1, fp);

        // Level
        pLevelConst level = pair.second;
        fwrite(&level->eigenvalue, sizeof(double), 1, fp);
        fwrite(&level->N, sizeof(unsigned int), 1, fp);
        fwrite(level->eigenvector, sizeof(double), level->N, fp);
        fwrite(&level->gFactor, sizeof(double), 1, fp);
    }

    fclose(fp);
}
