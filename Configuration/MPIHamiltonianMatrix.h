#ifndef MPI_HAMILTONIAN_MATRIX_H
#define MPI_HAMILTONIAN_MATRIX_H

#include "HamiltonianMatrix.h"

class MPIHamiltonianMatrix: public HamiltonianMatrix
{
public:
    MPIHamiltonianMatrix(const ExcitedStates& excited_states, const RelativisticConfigList& rconfigs);
    virtual ~MPIHamiltonianMatrix(void) {}

public:
    virtual void GenerateMatrix();
    virtual void PollMatrix();
    virtual void SolveMatrix(unsigned int num_solutions, unsigned int two_j, bool gFactor = false);

    virtual void GetEigenvalues() const;

protected:
    virtual void GetgFactors(unsigned int two_j, double* g_factors) const;

protected:
    RelativisticConfigList::const_iterator config_start, config_end;
    unsigned int M_start, M_end;
};

#endif
