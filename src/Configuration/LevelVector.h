#ifndef LEVEL_VECTOR_PRINT_H
#define LEVEL_VECTOR_PRINT_H

#include "HamiltonianID.h"
#include "NonRelConfiguration.h"
#include <Eigen/Eigen>
#include <iostream>

namespace Ambit
{
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMajorMatrix;

/** LevelVector consists of Hamiltonian that levels come from,
    common relativistic config list, and a set of levels including
    energies, coefficients, and g-factors.
    Internally the coefficients are stored as a matrix for efficient processing of operators.
 */
class LevelVector
{
public:
    LevelVector(pHamiltonianID hID = nullptr, pRelativisticConfigList rconfigs = nullptr):
        hID(hID), configs(rconfigs) {}
//    LevelVector(pHamiltonianID hID, pRelativisticConfigList configs, pLevel level):
//        hID(hID), configs(configs), levels(1, level) {}
//    LevelVector(pRelativisticConfigList configs, pLevel level): configs(configs), levels(1, level)
//    {
//        if(level)
//            hID = level->GetHamiltonianID();
//    }

    // All data is open, but don't break it. The number of rows in eigenvectors should match the number of elements
    // in eigenvectors and g_factors (unless g_factors is empty).
    // There are functions if you prefer to use them.
    pHamiltonianID hID {nullptr};
    pRelativisticConfigList configs {nullptr};
    RowMajorMatrix eigenvectors;
    std::vector<double> eigenvalues;
    std::vector<double> g_factors;

    inline void Resize(unsigned int num_levels, unsigned int num_csfs, bool allocate_gfactors = false)
    {
        eigenvectors.resize(num_levels, num_csfs);
        eigenvalues.resize(num_levels);
        if(allocate_gfactors)
            g_factors.resize(num_levels);
    }
    inline unsigned int NumLevels() const { return eigenvectors.rows(); }

    /** Return a new LevelVector including a subset of the current levels. */
    LevelVector getSubset(unsigned int start_index, unsigned int num_levels = 1) const
    {
        assert(start_index + num_levels <= eigenvalues.size());
        LevelVector subset(hID, configs);
        subset.Resize(num_levels, configs->NumCSFs());
        subset.eigenvectors = eigenvectors.middleRows(start_index, num_levels);
        std::copy(eigenvalues.begin()+start_index, eigenvalues.begin()+start_index+num_levels,subset.eigenvalues.begin());
        if(start_index + num_levels < g_factors.size())
        {
            subset.g_factors.resize(num_levels);
            std::copy(g_factors.begin()+start_index, g_factors.begin()+start_index+num_levels, subset.g_factors.begin());
        }
        return subset;
    }

    /** Print levels to outstream. min_percentages outside of the range (0, 100] will switch off the printing of configurations. */
    template<typename ConfigPrintType = NonRelConfiguration>
    inline void Print(double min_percentage = 1.0, std::optional<double> max_energy = std::nullopt) const;

    /** Print one line per level: J, P, index, E, g, Largest configuration, HamiltonianID. */
    template<typename ConfigPrintType = NonRelConfiguration>
    inline void PrintInline(bool print_leading_configuration = true, bool print_gfactors = true, std::string separator = " ", std::optional<double> max_energy = std::nullopt) const;
};

template<typename ConfigPrintType>
void LevelVector::Print(double min_percentage, std::optional<double> max_energy) const
{
    if(eigenvectors.size() == 0 || ProcessorRank != 0)
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

    *outstream << "Solutions for " << hID << " (N = " << configs->NumCSFs();
    if(configs->NumCSFsSmall() != configs->NumCSFs())
        *outstream << " x " << configs->NumCSFsSmall();
    *outstream << "):\n";

    for(unsigned int index = 0; index < NumLevels(); index++)
    {
        double energy = eigenvalues[index];
        if(max_energy && (energy > max_energy))
            break;

        *outstream << index << ": " << std::setprecision(8) << energy << "    "
                   << std::setprecision(12) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << " /cm\n";

        if(use_min_percentage)
        {
            // Clear percentages
            for(auto& pair: percentages)
                pair.second = 0.0;

            const double* eigenvector_csf = eigenvectors.row(index).data();
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

        if(hID->GetTwoJ() != 0 && index < g_factors.size())
        {
            if(!std::isnan(g_factors[index]))
                *outstream << "    g-factor = " << std::setprecision(5) << g_factors[index] << "\n";
        }

        *outstream << std::endl;
    }
}

template<typename ConfigPrintType>
void LevelVector::PrintInline(bool print_leading_configuration, bool print_gfactors, std::string separator, std::optional<double> max_energy) const
{
    if(eigenvectors.size() == 0 || ProcessorRank != 0)
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

    for(unsigned int index = 0; index < NumLevels(); index++)
    {
        double energy = eigenvalues[index];
        if(max_energy && (energy > max_energy))
            break;

        *outstream << J << separator << Sign(P) << separator << index << separator
                   << std::setprecision(12) << energy;

        if(print_gfactors)
        {
            if(hID->GetTwoJ() != 0 && g_factors.size())
            {
                *outstream << separator << std::setprecision(5) << g_factors[index];
            }
            else
                *outstream << separator << 0.0;
        }

        if(print_leading_configuration)
        {
            // Clear percentages
            for(auto& pair: percentages)
                pair.second = 0.0;

            const double* eigenvector_csf = eigenvectors.row(index).data();
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
    }
}

}
#endif
