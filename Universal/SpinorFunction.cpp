#include "Include.h"
#include "SpinorFunction.h"
#include "Interpolator.h"

SpinorFunction::SpinorFunction(int kappa, unsigned int size):
    kappa(kappa)
{
    if(size)
        resize(size);
}

SpinorFunction::SpinorFunction(const SpinorFunction& other):
    f(other.f), g(other.g), dfdr(other.dfdr), dgdr(other.dgdr), kappa(other.kappa)
{}

SpinorFunction::SpinorFunction(SpinorFunction&& other):
    f(other.f), g(other.g), dfdr(other.dfdr), dgdr(other.dgdr), kappa(other.kappa)
{}

void SpinorFunction::resize(unsigned int size)
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

SpinorFunction& SpinorFunction::operator=(SpinorFunction&& other)
{
    swap(other);
    kappa = other.kappa;

    return *this;
}

void SpinorFunction::swap(SpinorFunction& other)
{
    f.swap(other.f);
    g.swap(other.g);
    dfdr.swap(other.dfdr);
    dgdr.swap(other.dgdr);
}

const SpinorFunction& SpinorFunction::operator*=(double scale_factor)
{
    if(scale_factor != 1.)
        for(unsigned int i=0; i<size(); i++)
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
    if(size() < other.size())
        resize(other.size());

    for(unsigned int i = 0; i < other.size(); i++)
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

const SpinorFunction& SpinorFunction::operator*=(const RadialFunction& chi)
{
    // Outside range of chi, chi is assumed to be zero.
    if(chi.size() < size())
        resize(chi.size());

    for(unsigned int i = 0; i < size(); i++)
    {   f[i] *= chi.f[i];
        g[i] *= chi.f[i];
        dfdr[i] = f[i] * chi.dfdr[i] + dfdr[i] * chi.f[i];
        dgdr[i] = g[i] * chi.dfdr[i] + dgdr[i] * chi.f[i];
    }

    return *this;
}

SpinorFunction SpinorFunction::operator*(const RadialFunction& chi) const
{
    SpinorFunction ret(*this);
    return (ret *= chi);
}

RadialFunction SpinorFunction::GetDensity() const
{
    RadialFunction ret(size());

    for(unsigned int i = 0; i < ret.size(); i++)
    {   ret.f[i] = f[i] * f[i] + g[i] * g[i];
        ret.dfdr[i] = 2. * (f[i] * dfdr[i] + g[i] * dgdr[i]);
    }

    return ret;
}

RadialFunction SpinorFunction::GetDensity(const SpinorFunction& other) const
{
    RadialFunction ret(mmin(size(), other.size()));

    for(unsigned int i = 0; i < ret.size(); i++)
    {   ret.f[i] = f[i] * other.f[i] + g[i] * other.g[i];
        ret.dfdr[i] = f[i] * other.dfdr[i] + dfdr[i] * other.f[i]
                     +g[i] * other.dgdr[i] + dgdr[i] * other.g[i];
    }

    return ret;
}

void SpinorFunction::Write(FILE* fp) const
{
    fwrite(&kappa, sizeof(int), 1, fp);

    unsigned int my_size = size();
    fwrite(&my_size, sizeof(unsigned int), 1, fp);

    // Write each vector in one hit.
    // This is important when the code is run on the supercomputer.
    fwrite(f.data(), sizeof(double), my_size, fp);
    fwrite(g.data(), sizeof(double), my_size, fp);
    fwrite(dfdr.data(), sizeof(double), my_size, fp);
    fwrite(dgdr.data(), sizeof(double), my_size, fp);
}

void SpinorFunction::Read(FILE* fp)
{
    fread(&kappa, sizeof(int), 1, fp);

    unsigned int my_size;
    fread(&my_size, sizeof(unsigned int), 1, fp);
    resize(my_size);

    fread(f.data(), sizeof(double), my_size, fp);
    fread(g.data(), sizeof(double), my_size, fp);
    fread(dfdr.data(), sizeof(double), my_size, fp);
    fread(dgdr.data(), sizeof(double), my_size, fp);
}

RadialFunction::RadialFunction(const std::vector<double>& pf, const std::vector<double>& pdfdr)
{
    f = pf;
    dfdr = pdfdr;

    dfdr.resize(f.size());
}

RadialFunction::RadialFunction(unsigned int size)
{   Clear();
    resize(size);
}

const RadialFunction& RadialFunction::operator=(const RadialFunction& other)
{   f = other.f;
    dfdr = other.dfdr;
    return *this;
}

RadialFunction& RadialFunction::operator=(RadialFunction&& other)
{   f.swap(other.f);
    dfdr.swap(other.dfdr);
    return *this;
}

void RadialFunction::resize(unsigned int size)
{   f.resize(size);
    dfdr.resize(size);
}

void RadialFunction::Clear()
{   f.clear();
    dfdr.clear();
}

const RadialFunction& RadialFunction::operator*=(double scale_factor)
{
    for(unsigned int i = 0; i < f.size(); i++)
    {
        f[i] *= scale_factor;
        dfdr[i] *= scale_factor;
    }

    return *this;
}

RadialFunction RadialFunction::operator*(double scale_factor) const
{
    RadialFunction ret(*this);
    return (ret *= scale_factor);
}

const RadialFunction& RadialFunction::operator+=(const RadialFunction& other)
{
    if(size() < other.size())
        resize(other.size());
    
    for(unsigned int i = 0; i < other.size(); i++)
    {   f[i] += other.f[i];
        dfdr[i] += other.dfdr[i];
    }
    
    return *this;
}

const RadialFunction& RadialFunction::operator-=(const RadialFunction& other)
{
    if(size() < other.size())
        resize(other.size());
    
    for(unsigned int i = 0; i < other.size(); i++)
    {   f[i] -= other.f[i];
        dfdr[i] -= other.dfdr[i];
    }
    
    return *this;
}

RadialFunction RadialFunction::operator+(const RadialFunction& other) const
{
    RadialFunction ret(*this);
    return (ret += other);
}

RadialFunction RadialFunction::operator-(const RadialFunction& other) const
{
    RadialFunction ret(*this);
    return (ret -= other);
}

const RadialFunction& RadialFunction::operator*=(const RadialFunction& other)
{
    if(size() > other.size())
        resize(other.size());

    for(unsigned int i = 0; i < size(); i++)
    {
        f[i] *= other.f[i];
        dfdr[i] = f[i] * other.dfdr[i] + dfdr[i] * other.f[i];
    }

    return *this;
}

RadialFunction RadialFunction::operator*(const RadialFunction& other) const
{
    RadialFunction ret(*this);
    return (ret *= other);
}

bool RadialFunction::Print(const std::string& filename, pLattice lattice) const
{
    FILE* fp = fopen(filename.c_str(), "wt");

    if(fp)
    {   bool success = Print(fp, lattice);
        fclose(fp);
        return success;
    }
    else
        return false;
}

bool RadialFunction::Print(FILE* fp, pLattice lattice) const
{
    unsigned int i;
    if(lattice)
        for(i = 0; i < size(); i++)
        {
            fprintf(fp, "%12.6e %12.6e %12.6e\n", lattice->R(i), f[i], dfdr[i]);
        }
    else
        for(i = 0; i < size(); i++)
        {
            fprintf(fp, "%d %12.6e %12.6e\n", i, f[i], dfdr[i]);
        }

    return true;
}
