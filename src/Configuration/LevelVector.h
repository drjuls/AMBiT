#ifndef LEVEL_VECTOR_PRINT_H
#define LEVEL_VECTOR_PRINT_H

#include "Level.h"
#include "NonRelConfiguration.h"
#include <iostream>

namespace Ambit
{
/** LevelVector consists of Hamiltonian that levels come from,
    common relativistic config list, and vector of levels.
 */
class LevelVector
{
public:
    LevelVector(pHamiltonianID hID = nullptr): hID(hID) {}
    LevelVector(pHamiltonianID hID, pRelativisticConfigList configs, pLevel level):
        hID(hID), configs(configs), levels(1, level) {}
    LevelVector(pRelativisticConfigList configs, pLevel level): configs(configs), levels(1, level)
    {
        if(level)
            hID = level->GetHamiltonianID();
    }

    pHamiltonianID hID {nullptr};
    pRelativisticConfigList configs {nullptr};
    std::vector<pLevel> levels;

    /** Print levels to outstream. min_percentages outside of the range (0, 100] will switch off the printing of configurations. */
    template<typename ConfigPrintType = NonRelConfiguration>
    inline void Print(double min_percentage = 1.0) const;

    template<typename ConfigPrintType = NonRelConfiguration>
    inline void Print(double min_percentage, double max_energy) const;

    /** Print one line per level: J, P, index, E, g, Largest configuration, HamiltonianID. */
    template<typename ConfigPrintType = NonRelConfiguration>
    inline void PrintInline(bool print_leading_configuration = true, bool print_gfactors = true, std::string separator = " ") const;

    template<typename ConfigPrintType = NonRelConfiguration>
    inline void PrintInline(double max_energy, bool print_leading_configuration = true, bool print_gfactors = true, std::string separator = " ") const;

protected:
    /** Print LevelVector to outstream, with all possible options for printing.
        All other print functions call this one.
     */
    template<typename ConfigPrintType = NonRelConfiguration>
    inline void Print(double min_percentage, bool use_max_energy, double max_energy) const;

    template<typename ConfigPrintType = NonRelConfiguration>
    inline void PrintInline(bool use_max_energy, double max_energy, bool print_leading_configuration, bool print_gfactors, std::string separator = " ") const;
};

template<typename ConfigPrintType>
void LevelVector::Print(double min_percentage) const
{   Print<ConfigPrintType>(min_percentage, false, 0.0);
}

template<typename ConfigPrintType>
void LevelVector::Print(double min_percentage, double max_energy) const
{   Print<ConfigPrintType>(min_percentage, true, max_energy);
}

template<typename ConfigPrintType>
void LevelVector::Print(double min_percentage, bool use_max_energy, double max_energy) const
{
    if(levels.size() == 0 || ProcessorRank != 0)
        return;

    if(!configs || !configs->size())
        return;

    bool use_min_percentage = (0. < min_percentage && min_percentage <= 100.);

    // Build map of non-rel configurations to percentage contributions
    std::map<ConfigPrintType, double> percentages;
    if(use_min_percentage)
        for(auto& rconfig: *configs)
        {   percentages[rconfig] = 0.0; // Auto instantiate NonRelConfiguration from RelativisticConfiguration.
        }

    auto solution_it = levels.begin();
    *outstream << "Solutions for " << hID << " (N = " << configs->NumCSFs();
    if(configs->NumCSFsSmall() != configs->NumCSFs())
        *outstream << " x " << configs->NumCSFsSmall();
    *outstream << "):\n";

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
                {   contribution += gsl_pow_2(*eigenvector_csf);
                    eigenvector_csf++;
                }

                percentages[rconfig] += contribution * 100.;
            }

            // Print important configurations.
            *outstream << std::fixed << std::setprecision(4);
            for(auto& pair: percentages)
            {
                if(pair.second > min_percentage)
                    *outstream << std::setw(20) << pair.first.Name() << "  "
                               << pair.second << "%\n";
            }
            outstream->unsetf(std::ios_base::floatfield);
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

template<typename ConfigPrintType>
void LevelVector::PrintInline(bool print_leading_configuration, bool print_gfactors, std::string separator) const
{   PrintInline<ConfigPrintType>(false, 0.0, print_leading_configuration, print_gfactors, separator);
}

template<typename ConfigPrintType>
void LevelVector::PrintInline(double max_energy, bool print_leading_configuration, bool print_gfactors, std::string separator) const
{   PrintInline<ConfigPrintType>(true, max_energy, print_leading_configuration, print_gfactors, separator);
}

template<typename ConfigPrintType>
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
    std::map<ConfigPrintType, double> percentages;
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

}
#endif
