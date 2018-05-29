#ifndef MATRIX_H
#define MATRIX_H

namespace Ambit
{
/** Storage container for a square, symmetrical matrix of size N to be used with Davidson method.
    We need to use polymorphism rather than templating because Eigensolver must pass the method
    MatrixMultiply to the Davidson method wrapped in an "extern C" function.
 */
class Matrix
{
public:
    Matrix(): N(0) {}
    virtual ~Matrix() {}

    virtual unsigned int size() const { return N; }

    /** Multiply matrix by another matrix, size N*M.
        PRE: b and c are N * M matrices; WriteMode() returns false.
        POST: c = A * b, where A is *this.
        Since this routine is intended for use with fortran code,
        b and c are column-major.
     */
    virtual void MatrixMultiply(int m, double* b, double* c) const = 0;

    /** Get the diagonal of the matrix.
        PRE: diag = double[N], where N = this->GetSize().
     */
    virtual void GetDiagonal(double* diag) const = 0;

protected:
    unsigned int N;
};

}
#endif
