#ifndef SMALL_MATRIX_H
#define SMALL_MATRIX_H

#include "Matrix.h"

class SmallMatrix: public Matrix
{
    /** Storage container for a small square matrix, stored in memory. */
public:
    SmallMatrix(unsigned int size);
    virtual ~SmallMatrix(void);

    virtual void WriteMode(bool write);
    virtual void Clear();
    virtual double& At(unsigned int i, unsigned int j);

    virtual void MatrixMultiply(int m, double* b, double* c);

    double* GetMatrix() { return M; }

protected:
    /** Make the matrix symmetrical. Copy the upper triangle to the lower one. */
    void Symmetrise();

protected:
    double* M;  // Entire matrix
};

#endif
