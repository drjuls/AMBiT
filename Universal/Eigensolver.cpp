#include "Eigensolver.h"
#include "Configuration/SmallMatrix.h"
#include <iostream>

#include "f2c.h"
#include "clapack.h"

#define SMALL_LIM 400

static Matrix* aa;

extern "C"{
/** Method in Davidson.lib for calculation of some eigenvalues of a large (sparse) matrix.
    See file Davidson/Acpz.c for details on usage.
 */
extern int dvdson_(int (*op)(int*, int*, doublereal*, doublereal*),
    integer *n, integer *lim, doublereal *diag, integer *ilow, integer *ihigh, integer *iselec,
    integer *niv, integer *mblock,
	doublereal *crite, doublereal *critc, doublereal *critr, doublereal *ortho, integer *maxiter,
    doublereal *work, integer *iworksz, integer* iwork, integer* iiwsz,
    logical* hiend, integer *nloops, integer *nmv, integer *ierr);

/** User supplied routine used by Davidson method.
    PRE: b and c are n x m matrices
    POST: c = A x b, where A is the matrix that is being diagonalised.
 */
int op(int *n, int *m, doublereal* b, doublereal* c)
{
#ifdef WIN32
    std::cerr << ".";
#endif

    doublereal* cp = &c[0];
    const doublereal* cend = c + (*n)*(*m);
    while(cp < cend)
        *cp++ = 0.;

    aa->MatrixMultiply(*m, reinterpret_cast<double*>(b), reinterpret_cast<double*>(c));

    return 0;
}
}

void Eigensolver::SolveSmallSymmetric(double* matrix, double* eigenvalues, unsigned int N)
{
    if(N)
    {   double* work = new double[3*N];
        char jobs='V', uplo='U';
        integer order = N,
                leading_dimension = N,
                work_size = 3*N,
                info;

        dsyev_(&jobs, &uplo, &order, matrix, &leading_dimension, eigenvalues,
            &work[0], &work_size, &info);

        delete[] work;

        if(info != 0)
        {   std::cerr << "dsyev failed, info = " << info << std::endl;
            return;
        }
    }
}

void Eigensolver::SolveLargeSymmetric(Matrix* matrix, double* eigenvalues, double* eigenvectors, unsigned int N, unsigned int num_solutions)
{
    static integer n, lim;
    static doublereal *diag;
    static integer ilow, ihigh, *iselec;
    static integer niv, mblock, maxiter;
    static integer nloops, nmv, ierr;
    static logical hiend;
    static integer worksize, intworksize;
    static doublereal *work;
    static integer *intwork;
    static doublereal crite, critc, critr, ortho;

    SmallMatrix* sm = dynamic_cast<SmallMatrix*>(matrix);
    
    if(sm && (N <= SMALL_LIM))
    {
        double* M = sm->GetMatrix();
        double* E = new double[N];

        SolveSmallSymmetric(M, E, N);

        unsigned int i;
        for(i = 0; i < num_solutions; i++)
            eigenvalues[i] = E[i];
        for(i = 0; i < num_solutions*N; i++)
            eigenvectors[i] = M[i];

        delete[] E;
        return;
    }

    int i, j;

    aa = matrix;
    
    n = N;
    lim = min(N, num_solutions+20);
    diag = new doublereal[n];
    for(i=0; i<n; i++)
        diag[i] = aa->At(i, i);

    ilow = 1;
    ihigh = num_solutions;  // num eigenvalues wanted
    iselec = new integer[lim];    // unused if lowest eigenvalues are wanted
    for(i=0; i<lim; i++)
        iselec[i] = 0;
    niv = 0;
    mblock = num_solutions;
    maxiter = 2000;
    hiend = false;
    
    worksize = 2*n*lim + lim*lim + (ihigh+10)*lim + ihigh;
    work = new doublereal[worksize];

    intworksize = 7*lim;
    intwork = new integer[intworksize];
    
    crite = 1.e-14;  // Algorithm stops if ANY of these are satisfied
    critc = 1.e-12;  //1.0e-8;
    critr = 1.e-10;  //1.0e-10;
    ortho = 1.e-9;   //1.0e-6;

    dvdson_(op, &n, &lim, diag, &ilow, &ihigh, iselec, &niv, &mblock,
            &crite, &critc, &critr, &ortho, &maxiter,
            work, &worksize, intwork, &intworksize, &hiend, &nloops, &nmv, &ierr);

    if(ierr != 0)
    {   std::cerr << "dvdson failed, ierr = " << ierr << std::endl;
    }
    else
    {   for(i=0; i<num_solutions; i++)
        {   eigenvalues[i] = work[ihigh * n + i];
            for(j=0; j<n; j++)
                eigenvectors[i*n + j] = work[i*n + j];
        }
    }
    std::cout << "    nloops=" << nloops << std::endl;;

    delete[] diag;
    delete[] iselec;
    delete[] work;
    delete[] intwork;
}

bool Eigensolver::SolveSimultaneousEquations(double* matrix, double* vector, unsigned int N)
{
    if(N)
    {
        integer order = N,
                num_equations = 1,
                info;
        integer* work = new integer[N];

        // Transpose the matrix
        double* tmatrix = new double[N*N];
        for(unsigned int i=0; i<N; i++)
            for(unsigned int j=0; j<N; j++)
                tmatrix[i*N + j] = matrix[j*N + i];

        dgesv_(&order, &num_equations, tmatrix, &order, &work[0], vector, &order, &info);

        delete[] work;
        delete[] tmatrix;

        if(info == 0)
            return true;
        else
        {   std::cerr << "dgesv failed, info = " << info << std::endl;
        }
    }
    return false;
}

bool Eigensolver::SolveMatrixEquation(double* A_matrix, double* B_matrix, double* eigenvalues, unsigned int N)
{
    if(N)
    {
        integer itype = 1;
        char jobs='V', uplo='U';
        integer order = N,
                leading_dim_A = N,
                leading_dim_B = N,
                info;
        integer worksize = 3*N;
        double* work = new double[worksize];

        dsygv_(&itype, &jobs, &uplo, &order, &A_matrix[0], &leading_dim_A,
               &B_matrix[0], &leading_dim_B, &eigenvalues[0], &work[0], &worksize, &info);

        delete[] work;

        if(info == 0)
            return true;
        else
        {   std::cerr << "dsygv failed, info = " << info << std::endl;
        }
    }
    return false;
}
