#include "Include.h"
#include "MPIMatrix.h"

MPIMatrix::MPIMatrix(unsigned int size, unsigned int row_start, unsigned int row_end): Matrix(size), start(row_start), end(row_end)
{
    unsigned int i;

    try
    {   M = new double*[end - start];

        for(i = start; i < end; i++)
            M[i-start] = new double[N-i];
    }
    catch(std::bad_alloc& ba)
    {   *errstream << "MPIMatrix (N = " << N << "): " << ba.what() << std::endl;
        i -= start;
        while(i > 0)
        {   i--;
            delete[] M[i];
        }
        delete[] M;
        throw;
    }

    Clear();
}

MPIMatrix::~MPIMatrix(void)
{
    for(unsigned int i=start; i<end; i++)
        delete[] M[i-start];
    delete[] M;
}

void MPIMatrix::WriteMode(bool write)
{   write_mode = write;
}

void MPIMatrix::Clear()
{
    for(unsigned int i = 0; i < (end-start); i++)
        for(unsigned int j = 0; j < (N-i-start); j++)
            M[i][j] = 0.;
}

double& MPIMatrix::At(unsigned int i, unsigned int j)
{
#ifdef _DEBUG
    if((i < start) || (i >= end) || (j < i) || (j >= N))
        *errstream << "M(" << i << ", " << j << ")" << std::endl;
#endif
    if(j >= i)
        return M[i-start][j-i];
    else
        return M[j-start][i-j];
}

void MPIMatrix::MatrixMultiply(int m, double* b, double* c)
{
    unsigned int i, j;
    double temp;
    double B;
    double *rowp, *bp, *cp;   // Row iterators

    cp = &c[0];
    const double* cend = c + N * m;
    while(cp < cend)
        *cp++ = 0.;

    /* The basic plan for this follows Stathopoulos and Fischer's advice in their
       Davidson algorithm paper. The matrix is only accessed once.
       Upper half is dot product:
         c(i, j) += Sum_(k >= i) M(i, k) * b(k, j)
       Lower half is vector times scalar:
         c(k, j) += b(i, j) * M(i, k), for k > i

       PROOF:
         We transform the second one: i->k and k->i
             c(i, j) += b(k, j) * M(k, i), for i > k
         Then adding it to the upper half:
             c(i, j) = Sum(k >= i) M(i, k) * b(k, j) + Sum (k < i) b(k, j) * M(k, i)
                     = Sum(k) M(i, k) * b(k, j)
         as required.
     */
    for(i = start; i < end; i++)
    {
        double* row = M[i-start];     // current matrix row
        const double* rowend = row + (N - i);

        for(j = 0; j < (unsigned int)m; j++)
        {
            // Upper half is dot product:
            //   c(i, j) += Sum_(k >= i) M(i, k) * b(k, j)
            temp = c[i + j*N];
            rowp = &row[0];
            bp = &b[i+j*N];
            while(rowp < rowend)    // for(k = 0; k < N-i; k++)
                temp += *rowp++ * *bp++;
            c[i + j*N] = temp;

            // Lower half is vector times scalar:
            //   c(k, j) += b(i, j) * M(i, k), for k > i
            B = b[i + j*N];
            cp = &c[i+1+j*N];
            rowp = &row[1];
            while(rowp < rowend)    // for(k = 1; k < N-i; k++)
                *cp++ += B * *rowp++;
        }
    }
}
