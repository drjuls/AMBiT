#include "Include.h"
#include "Matrix.h"

#define MATRIX_EXISTS false

Matrix::Matrix(unsigned int size):
    N(size), write_mode(true)
{
    if(N < BIG_LIMIT)
    {
        try
        {   M = new double[N * N];
        }
        catch(std::bad_alloc& ba)
        {   std::cout << ba.what() << " N = " << N << std::endl;
            exit(1);
        }

        for(unsigned int i=0; i < N * N; i++)
            M[i] = 0.;
        M_start = 0;
        M_end = N * N;
        fp = NULL;
    }
    else
    {
        try
        {   M = new double[CHUNK_SIZE];
        }
        catch(std::bad_alloc& ba)
        {   std::cout << ba.what() << " Chunk = " << CHUNK_SIZE << std::endl;
            exit(1);
        }

        if(MATRIX_EXISTS)
        {   fp = fopen("temp.matrix", "r+b");
            fread(M, sizeof(double), CHUNK_SIZE, fp);
        }
        else
        {   for(unsigned int i=0; i < CHUNK_SIZE; i++)
                M[i] = 0.;
            fp = fopen("temp.matrix", "w+b");

            unsigned int amount_written = 0;
            while(amount_written < N*N)
            {
                if(amount_written < N*N - CHUNK_SIZE)
                {   fwrite(M, sizeof(double), CHUNK_SIZE, fp);
                    amount_written += CHUNK_SIZE;
                }
                else
                {   fwrite(M, sizeof(double), N*N - amount_written, fp);
                    amount_written = N*N;
                }
            }
        }

        M_start = 0;
        M_end = CHUNK_SIZE;
    }
}

Matrix::~Matrix()
{
    if(fp)
        fclose(fp);
    delete[] M;
}

void Matrix::GetChunk(unsigned int i)
{
    long start_point = ftell(fp);
    if(start_point == -1L)
        start_point = N*N;
    else
        start_point = start_point/sizeof(double);

    if((i > start_point) && (i < N*N - CHUNK_SIZE))
    {   start_point = (i - start_point) * sizeof(double);
        fseek(fp, start_point, SEEK_CUR);
    }
    else
    {   start_point = mmin(i, N*N - CHUNK_SIZE);
        start_point = start_point * sizeof(double);
        fseek(fp, start_point, SEEK_SET);
    }

    M_start = mmin(i, N*N - CHUNK_SIZE);

    fread(M, sizeof(double), CHUNK_SIZE, fp);

    M_end = M_start + CHUNK_SIZE;
}

void Matrix::StoreChunk()
{
    long start_point = long(M_start) * sizeof(double);

    fseek(fp, start_point, SEEK_SET);
    fwrite(M, sizeof(double), M_end - M_start, fp);
}

void Matrix::Symmetrise()
{
    unsigned int i, j;

    if(N < BIG_LIMIT)
    {
        for(i=0; i<N; i++)
            for(j=i+1; j<N; j++)
                M[j*N + i] = M[i*N + j];
        return;
    }

    double* upper_square = new double[1000 * 1000];
    double* lower_square = new double[1000 * 1000];

    unsigned int upper_start_row, upper_start_column;
    unsigned int& lower_start_row = upper_start_column;
    unsigned int& lower_start_column = upper_start_row;

    unsigned int upper_width, upper_height;
    unsigned int& lower_width = upper_height;
    unsigned int& lower_height = upper_width;

    long start_pos;

    // Loop through upper squares
    for(upper_start_row = 0; upper_start_row < N; upper_start_row+=1000)
    {
        upper_height = mmin(N - upper_start_row, (unsigned int)(1000));

        for(upper_start_column = upper_start_row; upper_start_column < N; upper_start_column += 1000)
        {
            upper_width = mmin(N - upper_start_column, (unsigned int)(1000));

            // Read upper_square
            for(i=0; i<upper_height; i++)
            {   start_pos = (upper_start_row + i) * N + upper_start_column;
                start_pos = start_pos * sizeof(double);
                
                if(fseek(fp, start_pos, SEEK_SET))
                {   printf("fseek failed\n");
                    getchar();
                    exit(1);
                }
                
                if(fread(&upper_square[i*upper_width], sizeof(double), upper_width, fp) != upper_width)
                {   printf("fread failed\n");
                    getchar();
                    exit(1);
                }
            }

            // Upper square == lower square (diagonal)
            if(upper_start_row == upper_start_column)
            {   
                // NOTE: upper_height == upper_width
                for(i=0; i<upper_height; i++)
                    for(j=i; j<upper_width; j++)
                    {   lower_square[i*lower_width + j] = upper_square[i*upper_width + j];
                        lower_square[j*lower_width + i] = upper_square[i*upper_width + j];
                    }
            }
            else
            {   for(i=0; i<upper_height; i++)
                    for(j=0; j<upper_width; j++)
                        lower_square[j*lower_width + i] = upper_square[i*upper_width + j];
            }

            // Write lower_square
            for(i=0; i<lower_height; i++)
            {   start_pos = (lower_start_row + i)*N + lower_start_column;
                start_pos = start_pos * sizeof(double);
                fseek(fp, start_pos, SEEK_SET);

                fwrite(&lower_square[i*lower_width], lower_width, sizeof(double), fp);
            }
        }
    }

    delete[] upper_square;
    delete[] lower_square;
}
