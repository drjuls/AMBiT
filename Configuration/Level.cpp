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
