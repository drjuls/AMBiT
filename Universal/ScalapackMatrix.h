#ifndef SCALAPACK_MATRIX_H
#define SCALAPACK_MATRIX_H

#include "Matrix.h"

class ScalapackMatrix: public Matrix
{
    /** Storage container for a matrix decomposed as per ScaLAPACK requirements.
        The entire matrix is stored in memory.
        The local array is stored in column-first (fortran style) format.
     */
public:
    ScalapackMatrix(unsigned int size);
    virtual ~ScalapackMatrix(void);

    virtual void WriteMode(bool write);
    virtual void Clear();
    virtual double& At(unsigned int i, unsigned int j);

    virtual void MatrixMultiply(int m, double* b, double* c) const;
    virtual void GetDiagonal(double* diag) const;

public:
    /** Read symmetric matrix from stored triangle.
        Upper triangles should be stored by row, which is identical to
        lower triangles stored by column.
     */
    virtual void ReadTriangle(const std::string& filename);

    /** Read lower part of symmetric matrix from stored triangle.
        This is the same as ReadTriangle(), but only stores the lower half,
        the upper half is set to zero.
     */
    virtual void ReadLowerTriangle(const std::string& filename);

    /** Write entire matrix to file, column by column.*/
    virtual void WriteToFile(const std::string& filename) const;

    /** Diagonalise matrix using scalapack.
        PRE: eigenvalues is of size N (this.GetSize()).
        POST: The matrix now stores the eigenvector instead of the original matrix.
     */
    virtual void Diagonalise(double* eigenvalues);

    /** Copy row from matrix.
        PRE: row is of size N (this.GetSize()).
     */
    virtual void GetRow(unsigned int row_number, double* row) const;

    /** Copy column from matrix.
        PRE: col is of size N (this.GetSize()).
     */
    virtual void GetColumn(unsigned int col_number, double* col) const;

protected:
    /** Checks whether the columns of the matrix V are eigenvalues by
        multiplying by M and calculating the norm and ratio to original.
        V is distributed in the same fashion as M.
        Writes results to logstream.
     */
    virtual void TestEigenvalues(const double* eigenvalues, const double* V) const;

protected:
    double* M;  // Pointer to local array
    
    // Descriptors for local array within global array
    unsigned int M_rows, M_cols;    // Size parameters
    unsigned int* M_row_numbers;    // Index of row and column numbers in global array
    unsigned int* M_col_numbers;

    // ScaLAPACK variables
    int ICTXT;      // BLACS context identifying the created process grid
    int DESC[9];    // Array descriptor
    int num_proc_rows, num_proc_cols;   // Number of processor rows and columns
    int proc_row, proc_col;             // Processor row and column of this node
};

#endif