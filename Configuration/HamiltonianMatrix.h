#ifndef HAMILTONIAN_MATRIX_H
#define HAMILTONIAN_MATRIX_H

#include "RelativisticConfiguration.h"
#include "Basis/ExcitedStates.h"
#include "Matrix.h"

class HamiltonianMatrix
{
public:
    HamiltonianMatrix(const ExcitedStates& excited_states, const RelativisticConfigList& rconfigs);
    virtual ~HamiltonianMatrix(void) 
    {   delete M;
        if(NumSolutions)
        {   delete[] E;
            delete[] V;
        }
    }

    void UpdateIntegrals();
    virtual void GenerateMatrix();
    virtual void PollMatrix();

    /** Solve the matrix that has been generated.
        If gFactors are required, set boolean to true.
     */
    virtual void SolveMatrix(unsigned int num_solutions, unsigned int two_j, bool gFactor = false);

    virtual void GetEigenvalues() const;

    void IncludeSMS_V2(bool include);

protected:
    /** Get the Hamiltonian matrix element between two projections. */
    double GetProjectionH(const Projection& first, const Projection& second) const;

    /** Get the Coulomb matrix element < e1, e2 | 1/r | e3, e4 >. */
    double CoulombMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4) const;

    /** Get the SMS matrix element between two projections. */
    double GetProjectionSMS(const Projection& first, const Projection& second) const;

    /** Get single particle SMS interaction with the core. */
    double CoreSMS(const ElectronInfo& e1, const ElectronInfo& e2) const;

    /** Get the SMS matrix element < e1, e2 | p.p | e3, e4 >. */
    double SMSMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4) const;

    double GetOneElectronIntegral(const StateInfo& s1, const StateInfo& s2) const;
    double GetTwoElectronIntegral(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4) const;

    /** Calculate the Lande g-factors for eigenfunctions of the Hamiltonian.
        PRE: gFactors.size >= NumSolutions.
     */
    virtual void GetgFactors(unsigned int two_j, double* g_factors) const;

    /** Get average z-component of spin for single electron. */
    double GetSz(const ElectronInfo& e) const;
    double GetSz(const ElectronInfo& e1, const ElectronInfo& e2) const;

protected:
    bool include_sms_v2;

    const RelativisticConfigList& configs;
    const ExcitedStates& states;

    unsigned int NumStates;
    std::map<StateInfo, unsigned int> state_index;

    // Storage for one and two electron integrals.
    // If these are null, it means that there is not enough space in memory to store them,
    // they must be generated as needed.

    // OneElectronIntegrals(i, j) = OneElectronIntegrals(i * NumStates + j) = <i|H|j>
    std::map<unsigned int, double> OneElectronIntegrals;

    // TwoElectronIntegrals(k, i, j, l, m) = R_k(ij, lm): i->l, j->m
    std::map<unsigned int, double> TwoElectronIntegrals;

    // SMSIntegrals(i, j) = <i|p|j>
    std::map<unsigned int, double> SMSIntegrals;

    Matrix* M;      // Hamiltonian Matrix
    unsigned int N; // Matrix M dimension = N x N

    double* E;      // Eigenvalues
    double* V;      // Eigenvectors
    unsigned int NumSolutions;
};

#endif
