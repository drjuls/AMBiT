#include "Include.h"
#include "SmallMatrix.h"

SmallMatrix::SmallMatrix(unsigned int size): Matrix(size)
{
    try
    {   M = new double[N*N];
    }
    catch(std::bad_alloc& ba)
    {   *errstream << "Small Matrix (N = " << N << "): " << ba.what() << std::endl;
        throw;
    }

    Clear();
}

SmallMatrix::~SmallMatrix(void)
{   delete[] M;
}

void SmallMatrix::WriteMode(bool write)
{
    // If the write mode is being switched off, symmetrise the matrix.
    if(write == false && write_mode == true)
        Symmetrise();

    write_mode = write;
}

void SmallMatrix::Clear()
{
    for(unsigned int i = 0; i < N*N; i++)
        M[i] = 0.;
}

double& SmallMatrix::At(unsigned int i, unsigned int j)
{   return M[i*N + j];
}

void SmallMatrix::MatrixMultiply(int m, double* b, double* c) const
{
    unsigned int i, j, k;
    double temp;

    for(i = 0; i < N; i++)
    {   for(j = 0; j < m; j++)
        {
            temp = 0.;
            for(k = 0; k < N; k++)
            {
                temp += M[i * N + k] * b[j * N + k];
            }
            c[j * N + i] = temp;
        }
    }
}

void SmallMatrix::Symmetrise()
{
    unsigned int i, j;

    for(i=0; i<N; i++)
        for(j=i+1; j<N; j++)
            M[j*N + i] = M[i*N + j];
}

void SmallMatrix::GetDiagonal(double* diag) const
{
    for(unsigned int i=0; i<N; i++)
        diag[i] = M[i*N + i];
}
