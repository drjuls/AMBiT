#ifndef MPI_MATRIX_H
#define MPI_MATRIX_H

#include "Universal/Matrix.h"
#include "RelativisticConfiguration.h"

class MPIMatrix: public Matrix
{
    /** Storage container for a distributed symmetrical matrix.
        Only the upper triangle is stored in memory.
        Each matrix stores a number of rows dependent on the RelativisticConfigList.
        Only rows that correspond to configs[i] are stored, where
            proc = i%(2*NumProcessors) - NumProcessors
        and 
           (proc == ProcessorRank) || (proc == -i - ProcessorRank)
     */
public:
    MPIMatrix(unsigned int size, const RelativisticConfigList& rlist);
    virtual ~MPIMatrix(void);

    virtual void WriteMode(bool write);
    virtual void Clear();

    // PRE: (i, j) is within the range of what this matrix stores.
    //      See class description and constructor for how to make sure of this.
    virtual double& At(unsigned int i, unsigned int j);

    virtual void MatrixMultiply(int m, double* b, double* c) const;
    virtual void GetDiagonal(double* diag) const;

protected:
    const RelativisticConfigList& configs;
    double** M;  // Pointers to matrix rows for this processor
};

#endif
