#ifndef COUPLED_FUNCTION_H
#define COUPLED_FUNCTION_H

#include <stdio.h>
#include <vector>

class CoupledFunction
{
public:
    CoupledFunction(unsigned int size = 0)
    {   Clear();
        if(size)
            ReSize(size);
    }
    CoupledFunction(const CoupledFunction& other): f(other.f), g(other.g), df(other.df), dg(other.dg) {}
    virtual ~CoupledFunction() {}

    std::vector<double> f, g,   // coupled functions (eg: upper and lower components)
                        df, dg; // function derivatives, ie. (df/dR).dR

    unsigned int Size() const { return static_cast<unsigned int>(f.size()); }

    void ReSize(unsigned int size)
    {   f.resize(size);
        g.resize(size);
        df.resize(size);
        dg.resize(size);
    }

    void Clear()
    {   f.clear();
        g.clear();
        df.clear();
        dg.clear();
    }

    /** Multiply all points of all vectors (f, g, df, dg) by the scale factor.
     */
    void Scale(double scale_factor)
    {   if(scale_factor != 1.)
            for(unsigned int i=0; i<Size(); i++)
            {   f[i] = f[i] * scale_factor;
                g[i] = g[i] * scale_factor;
                df[i] = df[i] * scale_factor;
                dg[i] = dg[i] * scale_factor;
            }
    }

    const CoupledFunction& operator=(const CoupledFunction& other)
    {   f = other.f;
        g = other.g;
        df = other.df;
        dg = other.dg;
        return *this;
    }

    /** Store the coupled function vectors.
        PRE: File pointer fp must be open and binary writable.
     */
    virtual void Write(FILE* fp) const;

    /** Read in a previously stored coupled function.
        PRE: fp must point to the start of a coupled function record.
     */
    virtual void Read(FILE* fp);
};

#endif
