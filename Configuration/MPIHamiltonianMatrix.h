#ifndef MPI_HAMILTONIAN_MATRIX_H
#define MPI_HAMILTONIAN_MATRIX_H

#include "HamiltonianMatrix.h"
#include "Configuration/Solution.h"

class MPIHamiltonianMatrix: public HamiltonianMatrix
{
public:
    MPIHamiltonianMatrix(const CIIntegrals& coulomb_integrals, ConfigGenerator* config_generator):
        HamiltonianMatrix(coulomb_integrals, config_generator)
    {}
    virtual ~MPIHamiltonianMatrix(void) {}

public:
    virtual void GenerateMatrix();
    virtual void WriteToFile(const std::string& filename, bool binary = true);
    virtual void PollMatrix();
    virtual void SolveMatrix(unsigned int num_solutions, Eigenstates& eigenstates, SolutionMap* aSolutionMapPointer, bool gFactors = false, bool TruncateDisplayAtMaxEnergy = false, double min_percentage = 1., double DavidsonMaxEnergy = 0.0);

    virtual void GetEigenvalues(const Eigenstates& eigenstates) const;

public:
#ifdef _SCALAPACK
    /** SolveScalapack works differently to the SolveMatrix:
        - All eigenvalues and eigenvectors are calculated
        - All eigenvectors with eigenvalues less than eigenvalue_limit (atomic units) are processed
          until max_num_solutions is reached. Use max_num_solutions = 0 to remove the limit.
        - Eigenvectors are processed in batches of size NumSolutions.
     */
    virtual void SolveScalapack(const std::string& filename, double eigenvalue_limit, Eigenstates& eigenstates, bool gFactors = false, unsigned int max_num_solutions = 50);
#endif

protected:
    virtual void GetgFactors(const Eigenstates& eigenstates, double* g_factors) const;
    virtual void GetgFactors(unsigned int NumSolutions, const double* V, unsigned int two_j, double* g_factors) const;
};

#endif
