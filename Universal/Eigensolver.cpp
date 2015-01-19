#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "Eigensolver.h"
#include "SmallMatrix.h"
#include "Configuration/MPIMatrix.h"

#define SMALL_LIM 1000

#if !(_FUS)
    #define dvdson_ dvdson
    #define dsyev_  dsyev
    #define dgesv_  dgesv
    #define dsygv_  dsygv
#endif

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
    aa->MatrixMultiply(*m, b, c);
    return 0;
}

/** Lapack routines */
void dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);
void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
void dsygv_(int*, char*, char*, int*, double*, int*, double*, int*, double*, double*, int*, int*);
}

#ifdef _MPI
MPI::Intracomm comm_world;
double* c_copy;

extern "C" {
int MPI_op(int *n, int *m, double* b, double* c)
{
    comm_world.Bcast(m, 1, MPI::INT, 0);
    comm_world.Bcast(b, (*n)*(*m), MPI::DOUBLE, 0);

    aa->MatrixMultiply(*m, b, c_copy);

    comm_world.Reduce(c_copy, c, (*n)*(*m), MPI::DOUBLE, MPI::SUM, 0);

    return 0;
}
}
#endif

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
        {   *errstream << "dsyev failed, info = " << info << std::endl;
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

    int i, j;

    aa = matrix;
    
    n = N;
    lim = mmin(N, num_solutions+20);
    diag = new double[n];
    aa->GetDiagonal(diag);

    ilow = 1;
    ihigh = num_solutions;  // num eigenvalues wanted
    iselec = new int[lim];    // unused if lowest eigenvalues are wanted
    for(i=0; i<lim; i++)
        iselec[i] = 0;
    niv = 0;
    mblock = num_solutions;
    maxiter = 20000;
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
    {   *errstream << "dvdson failed, ierr = " << ierr << std::endl;
    }
    else
    {   for(i=0; i<num_solutions; i++)
        {   eigenvalues[i] = work[ihigh * n + i];
            for(j=0; j<n; j++)
                eigenvectors[i*n + j] = work[i*n + j];
        }
    }
    *outstream << "    nloops=" << nloops << std::endl;;

    delete[] diag;
    delete[] iselec;
    delete[] work;
    delete[] intwork;
}

#ifdef _MPI
void Eigensolver::MPISolveLargeSymmetric(Matrix* matrix, double* eigenvalues, double* eigenvectors, unsigned int N, unsigned int num_solutions)
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

    int i, j;

    n = N;
    comm_world = MPI::COMM_WORLD;
    c_copy = new double[num_solutions * N]; 

    // Get diagonal
    diag = new double[n];
    double* my_diag = new double[n];
    matrix->GetDiagonal(my_diag);
    comm_world.Reduce(my_diag, diag, n, MPI::DOUBLE, MPI::SUM, 0);

    if(ProcessorRank == 0)
    {
        aa = matrix;
        lim = mmin(N, num_solutions+20);

        // Start davidson
        ilow = 1;
        ihigh = num_solutions;  // num eigenvalues wanted
        iselec = new int[lim];    // unused if lowest eigenvalues are wanted
        for(i=0; i<lim; i++)
            iselec[i] = 0;
        niv = 0;
        mblock = num_solutions;
        maxiter = 20000;
        hiend = false;
        
        worksize = 2*n*lim + lim*lim + (ihigh+10)*lim + ihigh;
        work = new double[worksize];
        for(i=0; i<worksize; i++)
            work[i] = 0.;

        intworksize = 7*lim;
        intwork = new int[intworksize];
        for(i=0; i<intworksize; i++)
            intwork[i] = 0;
        
        crite = 1.e-14;  // Algorithm stops if ANY of these are satisfied
        critc = 1.e-12;  //1.0e-8;
        critr = 1.e-10;  //1.0e-10;
        ortho = 1.e-9;   //1.0e-6;

        dvdson_(MPI_op, &n, &lim, diag, &ilow, &ihigh, iselec, &niv, &mblock,
                &crite, &critc, &critr, &ortho, &maxiter,
                work, &worksize, intwork, &intworksize, &hiend, &nloops, &nmv, &ierr);

        // send finish (m = 0)
        int finish_m = 0;
        comm_world.Bcast(&finish_m, 1, MPI::INT, 0);

        // send success
        comm_world.Bcast(&ierr, 1, MPI::INT, 0);

        if(ierr != 0)
        {   *errstream << "dvdson failed, ierr = " << ierr << std::endl;
            exit(1);
        }
        else
        {   for(i=0; i<num_solutions; i++)
            {   eigenvalues[i] = work[ihigh * n + i];
                for(j=0; j<n; j++)
                    eigenvectors[i*n + j] = work[i*n + j];
            }

            // broadcast results
            comm_world.Bcast(eigenvalues, num_solutions, MPI::DOUBLE, 0);
            comm_world.Bcast(eigenvectors, num_solutions * n, MPI::DOUBLE, 0);
        }

        delete[] iselec;
        delete[] work;
        delete[] intwork;
    }
    else    // worker nodes
    {   
        // Multiplications in davidson method
        double* b = new double[num_solutions * N];
        double* c = new double[num_solutions * N];

        nloops = 0;
        int m = 1;
        while(m != 0)
        {
            comm_world.Bcast(&m, 1, MPI::INT, 0);
            if(m != 0)
            {   nloops++;
                comm_world.Bcast(b, m * N, MPI::DOUBLE, 0);

                matrix->MatrixMultiply(m, b, c_copy);

                comm_world.Reduce(c_copy, c, m * N, MPI::DOUBLE, MPI::SUM, 0);                
            }
        }

        delete[] b;
        delete[] c;

        comm_world.Bcast(&ierr, 1, MPI::INT, 0);

        if(ierr != 0)
        {   *errstream << "dvdson failed, ierr = " << ierr << std::endl;
            exit(1);
        }
        else
        {   // get broadcast results
            comm_world.Bcast(eigenvalues, num_solutions, MPI::DOUBLE, 0);
            comm_world.Bcast(eigenvectors, num_solutions * N, MPI::DOUBLE, 0);
        }
    }

    delete[] c_copy;
    delete[] my_diag;
    delete[] diag;
    *outstream << "    nloops=" << nloops << std::endl;;
}
#endif

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
        {   *errstream << "dgesv failed, info = " << info << std::endl;
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
        {   *errstream << "dsygv failed, info = " << info << std::endl;
        }
    }
    return false;
}
