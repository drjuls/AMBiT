#ifndef MATRIX_H
#define MATRIX_H

class Matrix
{
    /** Storage container for a square, symmetrical matrix of size N*N. */
public:
    Matrix(unsigned int size): N(size) {}
    virtual ~Matrix(void) {}

    unsigned int GetSize() { return N; }

    /** Set/Clear write mode.
        WriteMode means that the matrix may still be undergoing changes.
        When WriteMode is turned off, Matrix should clean itself up and prepare for
        MatrixMultiply.
     */
    virtual void WriteMode(bool write) = 0;

    /** Reset the matrix to zero */
    virtual void Clear() = 0;

    /** Get reference to element. */
    virtual double& At(unsigned int i, unsigned int j) = 0;

    /** Multiply matrix by another matrix, size N*M.
        PRE: b and c are N * M matrices; c is initialised to zero.
        POST: c = A * b, where A is *this.
        Since this routine is intended for use with fortran code,
        the indices of b and c are in (column, row) order.
     */
    virtual void MatrixMultiply(int m, double* b, double* c) = 0;

protected:
    unsigned int N;
    bool write_mode;
};

#endif