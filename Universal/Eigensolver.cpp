#include "Include.h"
#include "Eigensolver.h"
#include "Configuration/SmallMatrix.h"

#define SMALL_LIM 400

static Matrix* aa;

extern "C"{
/** Method in Davidson.lib for calculation of some eigenvalues of a large (sparse) matrix.
    See file Davidson/Acpz.c for details on usage.
 */
extern int dvdson_(int (*op)(int*, int*, double*, double*),
    int *n, int *lim, double *diag, int *ilow, int *ihigh, int *iselec,
    int *niv, int *mblock,
	double *crite, double *critc, double *critr, double *ortho, int *maxiter,
    double *work, int *iworksz, int* iwork, int* iiwsz,
    bool* hiend, int *nloops, int *nmv, int *ierr);

/** User supplied routine used by Davidson method.
    PRE: b and c are n x m matrices
    POST: c = A x b, where A is the matrix that is being diagonalised.
 */
int op(int *n, int *m, double* b, double* c)
{
#ifdef WIN32
    std::cerr << ".";
#endif

    double* cp = &c[0];
    const double* cend = c + (*n)*(*m);
    while(cp < cend)
        *cp++ = 0.;

    aa->MatrixMultiply(*m, reinterpret_cast<double*>(b), reinterpret_cast<double*>(c));

    return 0;
}

/** Lapack routines */
void dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);

void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);

void dsygv_(int*, char*, char*, int*, double*, int*, double*, int*, double*, double*, int*, int*);
}

void Eigensolver::SolveSmallSymmetric(double* matrix, double* eigenvalues, unsigned int N)
{
    if(N)
    {   double* work = new double[3*N];
        char jobs='V', uplo='U';
        int order = N,
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
    static int n, lim;
    static double *diag;
    static int ilow, ihigh, *iselec;
    static int niv, mblock, maxiter;
    static int nloops, nmv, ierr;
    static bool hiend;
    static int worksize, intworksize;
    static double *work;
    static int *intwork;
    static double crite, critc, critr, ortho;

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
    lim = mmin(N, num_solutions+20);
    diag = new double[n];
    for(i=0; i<n; i++)
        diag[i] = aa->At(i, i);

    ilow = 1;
    ihigh = num_solutions;  // num eigenvalues wanted
    iselec = new int[lim];    // unused if lowest eigenvalues are wanted
    for(i=0; i<lim; i++)
        iselec[i] = 0;
    niv = 0;
    mblock = num_solutions;
    maxiter = 2000;
    hiend = false;
    
    worksize = 2*n*lim + lim*lim + (ihigh+10)*lim + ihigh;
    work = new double[worksize];

    intworksize = 7*lim;
    intwork = new int[intworksize];
    
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
        int order = N,
                num_equations = 1,
                info;
        int* work = new int[N];

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
        int itype = 1;
        char jobs='V', uplo='U';
        int order = N,
                leading_dim_A = N,
                leading_dim_B = N,
                info;
        int worksize = 3*N;
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
