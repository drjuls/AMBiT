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
    virtual void SolveMatrix(unsigned int num_solutions, unsigned int two_j, bool gFactors = false);

    virtual void GetEigenvalues() const;

public:
    /** SolveScalapack works differently to the SolveMatrix:
        - All eigenvalues and eigenvectors are calculated
        - All eigenvectors with eigenvalues less than eigenvalue_limit from
          the lowest eigenvalue are processed
        - Eigenvectors are processed in batches of size NumSolutions.
     */
    virtual void SolveScalapack(const std::string& filename, double eigenvalue_limit, unsigned int two_j, bool gFactors = false);

protected:
    virtual void GetgFactors(unsigned int two_j, double* g_factors) const;
};

#endif
