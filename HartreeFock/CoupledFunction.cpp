#include "Include.h"
#include "CoupledFunction.h"

void CoupledFunction::Write(FILE* fp) const
{
    unsigned int size = Size();
    fwrite(&size, sizeof(unsigned int), 1, fp);

    // Copy each vector to a buffer and then write in one hit.
    // This is important when the code is run on the supercomputer.
    double* buffer = new double[size];
    for(unsigned int i=0; i<4; i++)
    {
        const std::vector<double>* v;
        switch(i)
        {
        case 0:
            v = &f;
            break;
        case 1:
            v = &g;
            break;
        case 2:
            v = &df;
            break;
        case 3:
            v = &dg;
            break;
        }

        double* pbuffer = buffer;
        std::vector<double>::const_iterator it = v->begin();
        while(it != v->end())
        {   *pbuffer = *it;
            pbuffer++;
            it++;
        }
        fwrite(buffer, sizeof(double), size, fp);
    }

    delete[] buffer;
}

void CoupledFunction::Read(FILE* fp)
{
    unsigned int size;
    fread(&size, sizeof(unsigned int), 1, fp);
    ReSize(size);

    // Copy each vector to a buffer and then write in one hit.
    // This is important when the code is run on the supercomputer.
    double* buffer = new double[size];
    for(unsigned int i=0; i<4; i++)
    {
        fread(buffer, sizeof(double), size, fp);
        std::vector<double>* v;
        switch(i)
        {
        case 0:
            v = &f;
            break;
        case 1:
            v = &g;
            break;
        case 2:
            v = &df;
            break;
        case 3:
            v = &dg;
            break;
        }

        double* pbuffer = buffer;
        std::vector<double>::iterator it = v->begin();
        while(it != v->end())
        {   *it = *pbuffer;
            pbuffer++;
            it++;
        }
    }

    delete[] buffer;
}
