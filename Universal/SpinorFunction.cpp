#include "Include.h"
#include "SpinorFunction.h"

SpinorFunction::SpinorFunction(int kappa, unsigned int size):
    kappa(kappa)
{
    if(size)
        ReSize(size);
}

SpinorFunction::SpinorFunction(const SpinorFunction& other):
    f(other.f), g(other.g), dfdr(other.dfdr), dgdr(other.dgdr), kappa(other.kappa)
{}

void SpinorFunction::ReSize(unsigned int size)
{   f.resize(size);
    g.resize(size);
    dfdr.resize(size);
    dgdr.resize(size);
}

void SpinorFunction::Clear()
{   f.clear();
    g.clear();
    dfdr.clear();
    dgdr.clear();
}

const SpinorFunction& SpinorFunction::operator=(const SpinorFunction& other)
{
    f = other.f;
    g = other.g;
    dfdr = other.dfdr;
    dgdr = other.dgdr;

    kappa = other.kappa;

    return *this;
}

const SpinorFunction& SpinorFunction::operator*=(double scale_factor)
{
    if(scale_factor != 1.)
        for(unsigned int i=0; i<Size(); i++)
        {   f[i] = f[i] * scale_factor;
            g[i] = g[i] * scale_factor;
            dfdr[i] = dfdr[i] * scale_factor;
            dgdr[i] = dgdr[i] * scale_factor;
        }

    return *this;
}

SpinorFunction SpinorFunction::operator*(double scale_factor) const
{
    SpinorFunction ret(*this);
    return (ret *= scale_factor);
}

const SpinorFunction& SpinorFunction::operator+=(const SpinorFunction& other)
{
    if(Size() < other.Size())
        ReSize(other.Size());

    for(unsigned int i = 0; i < other.Size(); i++)
    {   f[i] += other.f[i];
        g[i] += other.g[i];
        dfdr[i] += other.dfdr[i];
        dgdr[i] += other.dgdr[i];
    }

    return *this;
}

const SpinorFunction& SpinorFunction::operator-=(const SpinorFunction& other)
{
    (*this) += other * (-1.0);
    return *this;
}

SpinorFunction SpinorFunction::operator+(const SpinorFunction& other) const
{
    SpinorFunction ret(*this);
    return ret += other;
}

SpinorFunction SpinorFunction::operator-(const SpinorFunction& other) const
{
    SpinorFunction ret(*this);
    return ret -= other;
}

const SpinorFunction& SpinorFunction::TimesEqualsVector(const std::vector<double>& chi, const std::vector<double>& dchidr)
{
    // Outside range of chi, chi is assumed to be zero.
    if(chi.size() < Size())
        ReSize(chi.size());

    for(unsigned int i = 0; i < Size(); i++)
    {   f[i] *= chi[i];
        g[i] *= chi[i];
        dfdr[i] = f[i] * dchidr[i] + dfdr[i] * chi[i];
        dgdr[i] = g[i] * dchidr[i] + dgdr[i] * chi[i];
    }

    return *this;
}

SpinorFunction SpinorFunction::TimesVector(const std::vector<double>& chi, const std::vector<double>& dchidr) const
{
    SpinorFunction ret(*this);
    return ret.TimesEqualsVector(chi, dchidr);
}

void SpinorFunction::Write(FILE* fp) const
{
    fwrite(&kappa, sizeof(int), 1, fp);

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
            v = &dfdr;
            break;
        case 3:
            v = &dgdr;
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

void SpinorFunction::Read(FILE* fp)
{
    fread(&kappa, sizeof(int), 1, fp);

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
            v = &dfdr;
            break;
        case 3:
            v = &dgdr;
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
