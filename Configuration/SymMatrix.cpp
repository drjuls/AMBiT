#include "Include.h"
#include "SymMatrix.h"

SymMatrix::SymMatrix(unsigned int size): Matrix(size)
{
    unsigned int i;

    try
    {   M = new double*[N];

        for(i=0; i<N; i++)
            M[i] = new double[N-i];
    }
    catch(std::bad_alloc& ba)
    {   std::cout << "SymMatrix (N = " << N << "): " << ba.what() << std::endl;
        while(i > 0)
        {   i--;
            delete[] M[i];
        }
        delete[] M;
        throw;
    }

    Clear();
}

SymMatrix::~SymMatrix(void)
{
    for(unsigned int i=0; i<N; i++)
        delete[] M[i];
    delete[] M;
}

void SymMatrix::WriteMode(bool write)
{   write_mode = write;
}

void SymMatrix::Clear()
{
    for(unsigned int i = 0; i < N; i++)
        for(unsigned int j = 0; j<N-i; j++)
            M[i][j] = 0.;
}

double& SymMatrix::At(unsigned int i, unsigned int j)
{
    if(j >= i)
        return M[i][j-i];
    else
        return M[j][i-j];
}

void SymMatrix::MatrixMultiply(int m, double* b, double* c)
{
    unsigned int i, j;
    double temp;
    double B;
    double *rowp, *bp, *cp;   // Row iterators

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
    for(i=0; i < N; i++)
    {
        double* row = M[i];     // current matrix row
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
