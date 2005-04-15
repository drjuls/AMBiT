#ifndef HAMILTONIAN_MATRIX_H
#define HAMILTONIAN_MATRIX_H

#include "RelativisticConfiguration.h"
#include "CIIntegrals.h"
#include "Basis/ExcitedStates.h"
#include "Universal/Matrix.h"

class HamiltonianMatrix
{
public:
    HamiltonianMatrix(const CIIntegrals& coulomb_integrals, const RelativisticConfigList& rconfigs);
    virtual ~HamiltonianMatrix(void) 
    {   if(M)
            delete M;
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
    virtual void GetgFactors(unsigned int two_j, double* g_factors) const;

    /** Get average z-component of spin for single electron. */
    double GetSz(const ElectronInfo& e) const;
    double GetSz(const ElectronInfo& e1, const ElectronInfo& e2) const;

protected:
    const RelativisticConfigList& configs;
    const CIIntegrals& integrals;

    Matrix* M;      // Hamiltonian Matrix
    unsigned int N; // Matrix M dimension = N x N

    double* E;      // Eigenvalues
    double* V;      // Eigenvectors
    unsigned int NumSolutions;
};

#endif
