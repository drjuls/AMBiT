#ifdef AMBIT_USE_MPI
#ifdef AMBIT_USE_SCALAPACK
#include <mpi.h>
#include "Include.h"
#include "ScalapackMatrix.h"
#include "MathConstant.h"
#include "Atom/Debug.h"

// Diagonaliser can use pdsyev or pdsyevd. pdsyevd is faster and
// guarantees orthogonality between eigenvectors, but requires more memory.
#define USE_PDSYEVD false

#if !(_FUS)
    #define sl_init_        sl_init
    #define blacs_gridinfo_ blacs_gridinfo
    #define numroc_         numroc
    #define descinit_       descinit
    #define blacs_gridexit_ blacs_gridexit
    #define pdsyevd_        pdsyevd
    #define pdsyev_         pdsyev
#endif

// Blocking parameters for division of global matrix.
// MB(NB) is the number of rows(columns) in the blocks used to distribute the global matrix.
const static int MB = 32;
const static int NB = 32;
// RSRC(CSRC) is the process row(col) over which the first row(col) of the global matrix is distributed.
// Set these to one so that root processor has the smallest local array possible.
// Unfortunately pdsyevd only works when they are set to zero.
const static int RSRC = 0;
const static int CSRC = 0;

extern "C"{
/** ScaLAPACK and BLACS routines */
void sl_init_(int*, const int*, const int*);
void blacs_gridinfo_(const int*, int*, int*, int*, int*);
int  numroc_(const int*, const int*, const int*, const int*, const int*);
void descinit_(int*, const int*, const int*, const int*, const int*, const int*, const int*,
               const int*, const int*, int*);
void blacs_gridexit_(const int*);
void pdsyevd_(const char* JOBZ, const char* UPLO, const int* N, double* A, const int* IA, const int* JA,
              const int* DESCA, double* W, double* Z, const int* IZ, const int* JZ, const int* DESCZ,
              double* WORK, const int* LWORK, int* IWORK, const int* LIWORK, int* INFO);
void pdsyev_(const char* JOBZ, const char* UPLO, const int* N, double* A, const int* IA, const int* JA,
             const int* DESCA, double* W, double* Z, const int* IZ, const int* JZ, const int* DESCZ,
             double* WORK, const int* LWORK, int* INFO);
}

ScalapackMatrix::ScalapackMatrix(unsigned int size):
    Matrix()
{
    N = size;

    // Get number of processor rows and columns
    num_proc_rows = (int)sqrt(double(NumProcessors));
    num_proc_cols = NumProcessors/num_proc_rows;

    while(num_proc_cols * num_proc_rows != NumProcessors)
    {   num_proc_rows--;
        num_proc_cols = NumProcessors/num_proc_rows;
    }
    *logstream << "ScalapackMatrix decomposition: " << num_proc_rows << " * " << num_proc_cols << std::endl;

    // Initialise process grid
    sl_init_(&ICTXT, &num_proc_rows, &num_proc_cols);
    blacs_gridinfo_(&ICTXT, &num_proc_rows, &num_proc_cols, &proc_row, &proc_col);

    // Get dimensions of local array
    M_rows = (unsigned int)numroc_((const int*)&N, &MB, &proc_row, &RSRC, &num_proc_rows);
    M_cols = (unsigned int)numroc_((const int*)&N, &NB, &proc_col, &CSRC, &num_proc_cols);

    int info;

    // Initialise array descriptor (DESC) for dense matrix
    descinit_(DESC, (const int*)&N, (const int*)&N, &MB, &NB, &RSRC, &CSRC, &ICTXT, (const int*)&M_rows, &info);

    M = new double[M_rows * M_cols];
    M_row_numbers = new unsigned int[M_rows];
    M_col_numbers = new unsigned int[M_cols];

    // Initialise row numbers
    unsigned int i, k;
    unsigned int M_i;
    i = (int(proc_row + num_proc_rows) - int(RSRC))%num_proc_rows * MB;

    if(DebugOptions.LogScalapack())
        *logstream << "ScalapackMatrix:\n Num rows/cols: " << M_rows << " " << M_cols
                   << "\n Row numbers: " << std::endl;

    M_i = 0; k = 0;
    while(M_i < M_rows)
    {
        M_row_numbers[M_i] = i + k;

        if(DebugOptions.LogScalapack())
            *logstream << " " << i + k;

        M_i++; k++;
        if(k >= MB)
        {   i = i + num_proc_rows*MB;
            k = 0;
        }
    }

    // Initialise column numbers
    i = (int(proc_col + num_proc_cols) - int(CSRC))%num_proc_cols * NB;

    if(DebugOptions.LogScalapack())
        *logstream << "\n Column numbers: " << std::endl;

    M_i = 0; k = 0;
    while(M_i < M_cols)
    {
        M_col_numbers[M_i] = i + k;

        if(DebugOptions.LogScalapack())
            *logstream << " " << i + k;

        M_i++; k++;
        if(k >= NB)
        {   i = i + num_proc_cols*NB;
            k = 0;
        }
    }

    if(DebugOptions.LogScalapack())
        *logstream << std::endl;
}

ScalapackMatrix::~ScalapackMatrix(void)
{
    // Release process grid
    blacs_gridexit_(&ICTXT);

    delete[] M;
    delete[] M_row_numbers;
    delete[] M_col_numbers;
}

void ScalapackMatrix::Clear()
{
    memset(M, 0, M_rows*M_cols*sizeof(double));
}

/** Multiply matrix by another matrix, size N*M.
    PRE: b and c are N * M matrices; c is initialised to zero.
    POST: c = A * b, where A is *this.
    Since this routine is intended for use with fortran code,
    the indices of b and c are in (column, row) order.
*/
void ScalapackMatrix::MatrixMultiply(int m, double* b, double* c) const
{
    MPI::Intracomm& comm_world = MPI::COMM_WORLD;

    double* buffer = new double[N];

    // Loop over columns of b
    for(unsigned int k = 0; k < m; k++)
    {
        // Fill buffer
        memset(buffer, 0, sizeof(double) * N);

        unsigned int M_i, M_j;
        M_i = 0;

        for(unsigned int i = 0; i < N; i++)
        {
            if(M_row_numbers[M_i] == i)
            {
                M_j = 0;
                unsigned int j_count = 0;
                unsigned int j = M_col_numbers[M_j];

                while(M_j < M_cols)
                {
                    buffer[i] += M[M_i + M_j*M_rows] * b[k*N + j + j_count];
                    M_j++; j_count++;
                    if(j_count >= NB)
                    {   j = j + num_proc_cols*NB;
                        j_count = 0;
                    }
                }

                M_i++;
            }

        }

        comm_world.Allreduce(buffer, &c[k * N], N, MPI::DOUBLE, MPI::SUM);
    }
}

void ScalapackMatrix::GetDiagonal(double* diag) const
{
}

void ScalapackMatrix::ReadTriangle(const std::string& filename)
{
    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
    {   *errstream << "ScalapackMatrix::ReadTriangle: Couldn't open " << filename << std::endl;
        exit(1);
    }

    // Check size of stored matrix matches
    unsigned int size;
    fread(&size, sizeof(unsigned int), 1, fp);

    if(size != N)
    {   *errstream << "ScalapackMatrix::ReadTriangle: Matrix size mismatch: " << filename
                   << " N = " << size << std::endl;
        exit(1);
    }

    // Read line by line and store our parts.
    unsigned int i, j;      // Position in global array
    unsigned int j_end;     // End of current block in global array
    unsigned int M_i, M_j;  // Position of start of line local array
    double* M_pos;          // Running position in local array

    double* buffer = new double[N];

    i = 0; M_i = 0; M_j = 0;
    while((i < N) && (M_i < M_rows) && (M_j < M_cols))
    {
        fread(buffer+i, sizeof(double), N-i, fp);

        // Copy row
        if(M_row_numbers[M_i] == i)
        {
            // Get first start/end point
            j = M_col_numbers[M_j];
            j_end = mmin((j/NB + 1)*NB, N);

            // First position in local array
            M_pos = &M[M_j * M_rows + M_i];

            while(j < j_end)
            {
                // need to use blas here
                for(unsigned int count =0; count < j_end-j; count++)
                    *(M_pos+count*M_rows) = buffer[j+count];
                M_pos += (j_end - j) * M_rows;
                j = j_end + (num_proc_cols - 1) * NB;
                j_end = mmin((j + NB), N);
            }
            // Increase M_i after checking columns
        }

        // Copy column
        if(M_col_numbers[M_j] == i)
        {
            // Get first start/end point
            j = M_row_numbers[M_i];
            j_end = mmin((j/MB + 1)*MB, N);

            M_pos = &M[M_j*M_rows + M_i];

            while(j < j_end)
            {
                memcpy(M_pos, &buffer[j], sizeof(double) * (j_end - j));
                M_pos += (j_end - j);
                j = j_end + (num_proc_rows - 1) * MB;
                j_end = mmin((j + MB), N);
            }

            M_j++;
        }

        if(M_row_numbers[M_i] == i)
            M_i++;

        i++;
    }

    delete[] buffer;
}

void ScalapackMatrix::ReadLowerTriangle(const std::string& filename)
{
    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
    {   *errstream << "ScalapackMatrix::ReadLowerTriangle: Couldn't open " << filename << std::endl;
        exit(1);
    }

    // Check size of stored matrix matches
    unsigned int size;
    fread(&size, sizeof(unsigned int), 1, fp);

    if(size != N)
    {   *errstream << "ScalapackMatrix::ReadLowerTriangle: Matrix size mismatch: " << filename
                   << " N = " << size << std::endl;
        exit(1);
    }

    // Read line by line and store our parts.
    unsigned int i, j;      // Position in global array
    unsigned int j_end;     // End of current block in global array
    unsigned int M_i, M_j;  // Position of start of line local array
    double* M_pos;          // Running position in local array

    double* buffer = new double[N];

    i = 0; M_i = 0; M_j = 0;
    while((i < N) && (M_i < M_rows) && (M_j < M_cols))
    {
        fread(buffer+i, sizeof(double), N-i, fp);

        // Copy column
        if(M_col_numbers[M_j] == i)
        {
            // Get first start/end point
            j = M_row_numbers[M_i];
            j_end = mmin((j/MB + 1)*MB, N);

            M_pos = &M[M_j*M_rows + M_i];

            while(j < j_end)
            {
                memcpy(M_pos, &buffer[j], sizeof(double) * (j_end - j));
                M_pos += (j_end - j);
                j = j_end + (num_proc_rows - 1) * MB;
                j_end = mmin((j + MB), N);
            }

            M_j++;
        }

        if(M_row_numbers[M_i] == i)
            M_i++;

        i++;
    }

    delete[] buffer;
}

void ScalapackMatrix::WriteToFile(const std::string& filename) const
{
    // Send columns to root node, which writes them sequentially.
    double* buffer = new double[N];
    double* writebuf;

    FILE* fp;
    if(ProcessorRank == 0)
    {
        fp = fopen(filename.c_str(), "wt");

        // Write size of matrix.
        fwrite(&N, sizeof(unsigned int), 1, fp);
        //fprintf(fp, "%i\n", N);
        writebuf = new double[N];
    }

    MPI::Intracomm& comm_world = MPI::COMM_WORLD;

    unsigned int i, j;
    unsigned int i_end;
    unsigned int M_i, M_j;
    double* M_pos;
    M_j = 0;

    for(j = 0; j < N; j++)
    {
        // Fill buffer
        memset(buffer, 0, sizeof(double) * N);

        if((M_j < M_cols) && (M_col_numbers[M_j] == j))
        {
            // Get first start/end point
            M_i = 0;
            i = M_row_numbers[M_i];
            i_end = mmin((i/MB + 1)*MB, N);
            M_pos = &M[M_j*M_rows + M_i];

            while(i < i_end)
            {
                memcpy(&buffer[i], M_pos, sizeof(double) * (i_end - i));
                M_pos += (i_end - i);
                i = i_end + (num_proc_rows - 1) * MB;
                i_end = mmin((i + MB), N);
            }

            M_i++;
        }

        comm_world.Reduce(buffer, writebuf, N, MPI::DOUBLE, MPI::SUM, 0);

        if(ProcessorRank == 0)
        {   fwrite(writebuf, sizeof(double), N, fp);
            //for(int count=0; count<N; count++)
            //        fprintf(fp, "%11.3e", writebuf[count]);
            //fprintf(fp, "\n");
        }

        M_j++;
    }

    delete[] buffer;
    if(ProcessorRank == 0)
    {   delete[] writebuf;
        fclose(fp);
    }
}

void ScalapackMatrix::Diagonalise(double* eigenvalues)
{
    char jobz = 'V';
    char uplo = 'L';    // Use upper or lower part of matrix: 'U' or 'L'
    int size = N;
    int IA = 1;         // Start position of submatrix
    int ierr;

    // Allocate the output (eigenvector) matrix
    int DESC_V[9];
    descinit_(DESC_V, (const int*)&N, (const int*)&N, &MB, &NB, &RSRC, &CSRC, &ICTXT, (const int*)&M_rows, &ierr);
    double* V = new double[M_rows * M_cols];

    double worksize;
    int lwork = -1;
    int liwork = 1;

#if USE_PDSYEVD
    // Calculate workspace size
    pdsyevd_(&jobz, &uplo, &size, M, &IA, &IA, DESC, eigenvalues, V, &IA, &IA,
             DESC_V, &worksize, &lwork, &liwork, &liwork, &ierr);
    lwork = (int)worksize;

    double* work = new double[lwork];
    int* iwork = new int[liwork];

    pdsyevd_(&jobz, &uplo, &size, M, &IA, &IA, DESC, eigenvalues, V, &IA, &IA,
            DESC_V, work, &lwork, iwork, &liwork, &ierr);

    delete[] work;
    delete[] iwork;
#else
    // Calculate workspace size
    pdsyev_(&jobz, &uplo, &size, M, &IA, &IA, DESC, eigenvalues, V, &IA, &IA,
            DESC_V, &worksize, &lwork, &ierr);
    lwork = (int)worksize;

    double* work = new double[lwork];

    pdsyev_(&jobz, &uplo, &size, M, &IA, &IA, DESC, eigenvalues, V, &IA, &IA,
            DESC_V, work, &lwork, &ierr);

    delete[] work;
#endif

    if(ierr != 0)
    {   *errstream << "scalapack::pdsyev failed, ierr = " << ierr << std::endl;
        exit(1);
    }

    if(DebugOptions.LogScalapack())
    {   // Restore matrix and test eigenvalues.
        ReadTriangle("MgI001.matrix");
        TestEigenvalues(eigenvalues, V);
    }

    // Make M the eigenvector matrix
    delete[] M;
    memcpy(DESC, DESC_V, sizeof(int) * 9);
    M = V;
}

void ScalapackMatrix::GetRow(unsigned int row_number, double* row) const
{
    if(row_number >= N)
        return;

    double* buffer = new double[N];

    MPI::Intracomm& comm_world = MPI::COMM_WORLD;

    unsigned int i, j;
    unsigned int M_i, M_j;

    // Reset buffer
    memset(buffer, 0, sizeof(double) * N);

    for(M_i = 0; M_i < M_rows; M_i++)
    {
        if(M_row_numbers[M_i] == row_number)
            break;
    }

    if(M_i < M_rows)
    {   // Fill our part of buffer
        M_j = 0;
        unsigned int k = 0;
        j = M_col_numbers[M_j];

        while(M_j < M_cols)
        {   buffer[j+k] = M[M_i + M_j*M_rows];
            M_j++; k++;
            if(k >= NB)
            {   j = j + num_proc_cols*NB;
                k = 0;
            }
        }
    }

    comm_world.Allreduce(buffer, row, N, MPI::DOUBLE, MPI::SUM);

    delete[] buffer;
}

void ScalapackMatrix::GetColumn(unsigned int col_number, double* col) const
{
    if(col_number >= N)
        return;

    double* buffer = new double[N];

    MPI::Intracomm& comm_world = MPI::COMM_WORLD;

    unsigned int i, j;
    unsigned int i_end;
    unsigned int M_i, M_j;
    double* M_pos;

    // Reset buffer
    memset(buffer, 0, sizeof(double) * N);

    for(M_j = 0; M_j < M_cols; M_j++)
    {
        if(M_col_numbers[M_j] == col_number)
            break;
    }

    // Fill our part of buffer
    if(M_j < M_cols)
    {
        // Get first start/end point
        M_i = 0;
        i = M_row_numbers[M_i];
        i_end = mmin((i/MB + 1)*MB, N);
        M_pos = &M[M_j*M_rows + M_i];

        while(i < i_end)
        {
            memcpy(&buffer[i], M_pos, sizeof(double) * (i_end - i));
            M_pos += (i_end - i);
            i = i_end + (num_proc_rows - 1) * MB;
            i_end = mmin((i + MB), N);
        }
    }

    comm_world.Allreduce(buffer, col, N, MPI::DOUBLE, MPI::SUM);

    delete[] buffer;
}

void ScalapackMatrix::TestEigenvalues(const double* eigenvalues, const double* V) const
{
    *logstream << "\nTesting Eigenvalues:" << std::endl;
    double* buffer = new double[N];
    double* column = new double[N];
    double* MtimesV = new double[N];

    MPI::Intracomm& comm_world = MPI::COMM_WORLD;

    for(unsigned int j = 0; j < N; j++)
    {
        // Get column j from V
        unsigned int i;
        unsigned int i_end;
        unsigned int V_i, V_j;
        const double* V_pos;

        memset(buffer, 0, sizeof(double) * N);

        for(V_j = 0; V_j < M_cols; V_j++)
        {
            if(M_col_numbers[V_j] == j)
                break;
        }

        if(V_j < M_cols)
        {
            // Get first start/end point
            V_i = 0;
            i = M_row_numbers[V_i];
            i_end = mmin((i/MB + 1)*MB, N);
            V_pos = &V[V_j*M_rows + V_i];

            while(i < i_end)
            {
                memcpy(&buffer[i], V_pos, sizeof(double) * (i_end - i));
                V_pos += (i_end - i);
                i = i_end + (num_proc_rows - 1) * MB;
                i_end = mmin((i + MB), N);
            }
        }

        comm_world.Allreduce(buffer, column, N, MPI::DOUBLE, MPI::SUM);

        // Mulitiply by M
        memset(MtimesV, 0, sizeof(double) * N);
        MatrixMultiply(1, column, MtimesV);

        // Test norm of M.v
        unsigned int k;

        double norm = 0.0;
        for(k=0; k < N; k++)
            norm += MtimesV[k]*MtimesV[k];
        norm = sqrt(norm);

        double average = 0.0;
        int av_count = 0;
        for(k=0; k < N; k++)
        {   if(fabs(column[k]) > 1.e-15)
            {   av_count++;
                average += MtimesV[k]/column[k];
            }
        }
        if(av_count > 0)
            average = average/av_count;

        double deviation = 0.0;
        for(k=0; k < N; k++)
        {   if(fabs(column[k]) > 1.e-15)
                deviation += pow(MtimesV[k]/column[k] - eigenvalues[j], 2.);
        }

        *logstream << j << ":\t" << eigenvalues[j] << " \tnorm = " << norm
                        << " \tav = " << average << " \tdev = " << deviation << std::endl;
    }

    delete[] buffer;
    delete[] column;
    delete[] MtimesV;
}

#endif
#endif
