#ifndef SPINOR_FUNCTION_H
#define SPINOR_FUNCTION_H

#include <stdio.h>
#include <vector>

/** Spinor function is the radial part of a Dirac function of the form
      psi = ( f(r) Omega_{kappa m} )
            (ig(r) Omega_{-kappa m})
    where Omega_{kappa m} contains the angular and spin variables.
 */
class SpinorFunction
{
public:
    SpinorFunction(int kappa, unsigned int size = 0);
    SpinorFunction(const SpinorFunction& other);
    virtual ~SpinorFunction() {}

    std::vector<double> f, g,       // Upper and lower components of the radial functions
                        dfdr, dgdr; // Function derivatives, ie. df/dr, dg/dr

    virtual unsigned int Size() const { return static_cast<unsigned int>(f.size()); }
    int Kappa() const { return kappa; }
    unsigned int L() const;
    double J() const;
    unsigned int TwoJ() const;

    /** Resize the functions. Pad with zeros if necessary. */
    virtual void ReSize(unsigned int size);
    virtual void Clear();

    const SpinorFunction& operator=(const SpinorFunction& other);

    /** Multiply all points of all vectors (f, g, df, dg) by the scale factor. */
    const SpinorFunction& operator*=(double scale_factor);
    SpinorFunction operator*(double scale_factor) const;

    /** Adding or subtracting two spinor functions can only occur if both have same angular part. */
    const SpinorFunction& operator+=(const SpinorFunction& other);
    const SpinorFunction& operator-=(const SpinorFunction& other);
    SpinorFunction operator+(const SpinorFunction& other) const;
    SpinorFunction operator-(const SpinorFunction& other) const;

    /** Multiply spinor function by a vector (assumed zero outside range).
        PRE: chi.size() == dchidr.size()
     */
    const SpinorFunction& TimesEqualsVector(const std::vector<double>& chi, const std::vector<double>& dchidr);
    SpinorFunction TimesVector(const std::vector<double>& chi, const std::vector<double>& dchidr) const;

    /** Store the coupled function vectors.
        PRE: File pointer fp must be open and binary writable.
     */
    virtual void Write(FILE* fp) const;

    /** Read in a previously stored coupled function.
        PRE: fp must point to the start of a coupled function record.
     */
    virtual void Read(FILE* fp);

protected:
    int kappa;
};

inline unsigned int SpinorFunction::L() const
{   if (kappa > 0)
        return (unsigned int)(kappa);
    else
        return (unsigned int)(-kappa-1);
}

inline double SpinorFunction::J() const
{   return (double)abs(kappa) - 0.5;
}

inline unsigned int SpinorFunction::TwoJ() const
{   return (2*abs(kappa) - 1);
}

#endif
