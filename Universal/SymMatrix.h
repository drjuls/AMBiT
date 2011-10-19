#ifndef SYM_MATRIX_H
#define SYM_MATRIX_H

#include "Matrix.h"

class SymMatrix: public Matrix
{
    /** Storage container for a symmetrical matrix.
        Only the upper triangle is stored in memory.
     */
public:
    SymMatrix(unsigned int size);
    virtual ~SymMatrix(void);

    virtual void WriteMode(bool write);
    virtual void Clear();
    virtual double& At(unsigned int i, unsigned int j);

    virtual void MatrixMultiply(int m, double* b, double* c) const;
    virtual void GetDiagonal(double* diag) const;

public:
    double** M;  // Pointers to matrix rows (N rows)
};

#endif
