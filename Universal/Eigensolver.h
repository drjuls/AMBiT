#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include "Matrix.h"

class Eigensolver
{
public:
    Eigensolver() {}
    ~Eigensolver() {}

    /** Solve a double symmetric matrix (using lapack routine "dsyev").
        PRE: matrix[N*N]
             eigenvalues[N]
        POST: matrix[i][j], j=(0, N-1) is the eigenvector of the original matrix
              with eigenvalue "eigenvalues[i]".
              Eigenvalues are sorted in ascending order.
     */
    void SolveSmallSymmetric(double* matrix, double* eigenvalues, unsigned int N);

    /** Solve a double symmetric matrix using Davidson algorithm.
        PRE: Matrix.GetSize() == N
             eigenvalues[num_solutions], eigenvectors[num_solutions * N]
             Only calculates lowest num_solutions
        POST: eigenvectors[i*N + j], i=(0, num_solutions-1), j=(0, N-1) is the eigenvector 
              of the original matrix with eigenvalue "eigenvalues[i]".
              Eigenvalues are sorted in ascending order.
     */
    void SolveLargeSymmetric(Matrix* matrix, double* eigenvalues, double* eigenvectors, unsigned int N, unsigned int num_solutions);

    /** Solve a double symmetric matrix using Davidson algorithm on a distributed architecture.
        PRE: Matrix.GetSize() == N
             eigenvalues[num_solutions], eigenvectors[num_solutions * N]
             Only calculates lowest num_solutions
        POST: eigenvectors[i*N + j], i=(0, num_solutions-1), j=(0, N-1) is the eigenvector 
              of the original matrix with eigenvalue "eigenvalues[i]".
              Eigenvalues are sorted in ascending order.
     */
    void MPISolveLargeSymmetric(Matrix* matrix, double* eigenvalues, double* eigenvectors, unsigned int N, unsigned int num_solutions);

    /** Solve a matrix equation in the form A*x = B, using lapack routine "dgesv".
        PRE: A = matrix[N][N]
             B = vector[N]
        POST: if(return == true)
                  vector[N] = x is the solution of the equations.
              otherwise routine failed (possibly no solution exists).
     */
    bool SolveSimultaneousEquations(double* matrix, double* vector, unsigned int N);

    /** Solve a matrix equation in the form A*x = e*B*x, using lapack routine "dgysv".
        PRE: A = matrix[N][N]
             B = matrix[N][N]
        POST: if(return == true)
                  e = vector[N] is a set of eigenvalues of the equation.
              matrix[i][j], j=(0, N-1) is the eigenvector (x) of the equation
                  with eigenvalue e[i].
              Eigenvalues are sorted in ascending order.                  
     */
    bool SolveMatrixEquation(double* A_matrix, double* B_matrix, double* eigenvalues, unsigned int N);
};

#endif
