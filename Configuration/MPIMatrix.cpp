#include "Include.h"
#include "MPIMatrix.h"

MPIMatrix::MPIMatrix(unsigned int size, const RelativisticConfigList& rlist): Matrix(size), configs(rlist)
{
    // proc runs from (-NumProcessors) to (NumProcessors-1) cyclically.
    // Our part of the matrix corresponds to proc == rank1 or rank2
    int proc, i;
    int rank1 = int(ProcessorRank);
    int rank2 = -1 - int(ProcessorRank);

    try
    {   M = new double*[N];
        for(i = 0; i < N; i++)
            M[i] = NULL;

        RelativisticConfigList::const_iterator it = configs.begin();
        proc = - int(NumProcessors);
        i = 0;
        while(it != configs.end())
        {
            if((proc == rank1) || (proc == rank2))
            {   for(unsigned int j = 0; j < it->NumJStates(); j++)
                {   M[i] = new double[N-i];
                    i++;
                }
            }
            else
                i += it->NumJStates();

            proc++;
            if(proc == int(NumProcessors))
                proc = -proc;
            
            it++;
        }
    }
    catch(std::bad_alloc& ba)
    {   *errstream << "MPIMatrix (N = " << N << "): " << ba.what() << std::endl;
        while(i > 0)
        {   i--;
            if(M[i])
                delete[] M[i];
        }
        delete[] M;
        throw;
    }

    Clear();
}

MPIMatrix::~MPIMatrix(void)
{
    for(unsigned int i=0; i<N; i++)
        if(M[i])
            delete[] M[i];
    delete[] M;
}

void MPIMatrix::WriteMode(bool write)
{   write_mode = write;
}

void MPIMatrix::Clear()
{
    for(unsigned int i = 0; i < N; i++)
        if(M[i])
            for(unsigned int j = 0; j < (N-i); j++)
                M[i][j] = 0.;
}

double& MPIMatrix::At(unsigned int i, unsigned int j)
{
#ifdef _DEBUG
    if((i < 0) || (i >= N) || (j < i) || (j >= N))
        *errstream << "M(" << i << ", " << j << ")" << std::endl;
#endif
    if(j >= i)
        return M[i][j-i];
    else
        return M[j][i-j];
}

void MPIMatrix::MatrixMultiply(int m, double* b, double* c) const
{
    unsigned int i, j;
    double temp;
    double B;
    double *row, *rowend, *rowp, *bp, *cp;   // Row iterators

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
    for(i = 0; i < N; i++)
    {
        if(row = M[i])  // current matrix row is part of our process
        {
            rowend = row + (N - i);

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
}

void MPIMatrix::GetDiagonal(double* diag) const
{
    for(unsigned int i=0; i<N; i++)
    {   if(M[i])
            diag[i] = *M[i];
        else
            diag[i] = 0.;
    }
}
