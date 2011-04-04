#ifndef HAMILTONIAN_MATRIX_H
#define HAMILTONIAN_MATRIX_H

#include "RelativisticConfiguration.h"
#include "CIIntegrals.h"
#include "Basis/ExcitedStates.h"
#include "Configuration/Solution.h"
#include "MBPT/Sigma3Calculator.h"
#include "Universal/Enums.h"
#include "Universal/Matrix.h"

#include "ConfigFileGenerator.h"
#include "Eigenstates.h"

class SolutionMap;

class HamiltonianMatrix
{
public:
    HamiltonianMatrix(const CIIntegrals& coulomb_integrals, ConfigGenerator* config_generator);
    virtual ~HamiltonianMatrix(void)
    {   if(M)
            delete M;
    }

    void UpdateIntegrals();
    virtual void GenerateMatrix();
    virtual void WriteToFile(const std::string& filename);
    virtual void PollMatrix();

    /** Solve the matrix that has been generated.
        If gFactors are required, set boolean to true.
        min_percentage is the threshold whereby leading configurations are printed
     */
    virtual void SolveMatrix(unsigned int num_solutions, Eigenstates& eigenstates, SolutionMap* aSolutionMapPointer, bool gFactors = false, double min_percentage = 1.);

    virtual void GetEigenvalues(const Eigenstates& eigenstates) const;

public:
    /** Include Sigma3 in the leading configurations. */
    void IncludeSigma3(Sigma3Calculator* sigma3)
    {   if(sigma3 && (configs->front().NumParticles() >= 3))
            include_sigma3 = true;
        sigma3calc = sigma3;
    }

protected:
    /** Get the Hamiltonian matrix element between two projections. */
    double GetProjectionH(const Projection& first, const Projection& second) const;

    /** Get the Coulomb matrix element < e1, e2 | 1/r | e3, e4 >. */
    double CoulombMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4) const;

    /** Get the SMS matrix element between two projections. */
    double GetProjectionSMS(const Projection& first, const Projection& second) const;

    /** Get the SMS matrix element < e1, e2 | p.p | e3, e4 >. */
    double SMSMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4) const;

    /** Calculate the Lande g-factors for eigenfunctions of the Hamiltonian.
        PRE: gFactors.size >= NumSolutions.
     */
    virtual void GetgFactors(const Eigenstates& eigenstates, double* g_factors) const;

    /** Get average z-component of spin for single electron. */
    double GetSz(const ElectronInfo& e) const;
    double GetSz(const ElectronInfo& e1, const ElectronInfo& e2) const;

protected:
    bool include_sigma3;
    Sigma3Calculator* sigma3calc;
    std::set<Configuration> leading_configs;

    /** Return value of Sigma3 for matrix element (added in GetProjectionH).
        Checks to see that first and second are leading configurations.
      */
    double GetSigma3(const Projection& first, const Projection& second) const;

    /** Get Sigma3(e1 e2 e3 -> e4 e5 e6) including all permutations.
        Sigma3 does the pair matching (matching e1 to e4, e5, and then e6).
     */
    double Sigma3(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
                  const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const;

    /** This function does the line permutations, putting the pairs on different levels
        of the three-body interaction.
     */
    inline double Sigma3LinePermutations(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
                  const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const;

protected:
    /** ConfigGenerator* confgen stores leading configs and the configuration list. */
    ConfigGenerator* confgen;
    /** RelativisticConfigList* configs is just a pointer to the list stored in confgen. */
    const RelativisticConfigList* configs;
    const CIIntegrals& integrals;

    Matrix* M;      // Hamiltonian Matrix
    unsigned int N; // Matrix M dimension = N x N
};

#endif
