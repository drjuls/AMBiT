#ifndef MPI_MATRIX_H
#define MPI_MATRIX_H

#include "Matrix.h"

class MPIMatrix: public Matrix
{
    /** Storage container for a distributed symmetrical matrix.
        Only the upper triangle is stored in memory.
     */
public:
    /** Stores rows from start to (end - 1). */
    MPIMatrix(unsigned int size, unsigned int row_start, unsigned int row_end);
    virtual ~MPIMatrix(void);

    virtual unsigned int StartRow() { return start; }
    virtual unsigned int EndRow() { return end; }

    virtual void WriteMode(bool write);
    virtual void Clear();
    virtual double& At(unsigned int i, unsigned int j);

    virtual void MatrixMultiply(int m, double* b, double* c);

protected:
    unsigned int start, end;
    double** M;  // Pointers to matrix rows (end - start rows)
};

#endif
