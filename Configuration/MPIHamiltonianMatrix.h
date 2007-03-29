#ifndef MPI_HAMILTONIAN_MATRIX_H
#define MPI_HAMILTONIAN_MATRIX_H

#include "HamiltonianMatrix.h"

class MPIHamiltonianMatrix: public HamiltonianMatrix
{
public:
    MPIHamiltonianMatrix(const CIIntegrals& coulomb_integrals, ConfigGenerator* config_generator):
        HamiltonianMatrix(coulomb_integrals, config_generator)
    {}
    virtual ~MPIHamiltonianMatrix(void) {}

public:
    virtual void GenerateMatrix();
    virtual void WriteToFile(const std::string& filename);
    virtual void PollMatrix();
    virtual void SolveMatrix(unsigned int num_solutions, Eigenstates& eigenstates, bool gFactors = false);

    virtual void GetEigenvalues(const Eigenstates& eigenstates) const;

public:
    /** SolveScalapack works differently to the SolveMatrix:
        - All eigenvalues and eigenvectors are calculated
        - All eigenvectors with eigenvalues less than eigenvalue_limit (atomic units) from
          the lowest eigenvalue are processed
        - Eigenvectors are processed in batches of size NumSolutions.
     */
    virtual void SolveScalapack(const std::string& filename, double eigenvalue_limit, Eigenstates& eigenstates, bool gFactors = false);

protected:
    virtual void GetgFactors(const Eigenstates& eigenstates, double* g_factors) const;
    virtual void GetgFactors(unsigned int NumSolutions, const double* V, unsigned int two_j, double* g_factors) const;
};

#endif
