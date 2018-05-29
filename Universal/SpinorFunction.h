#ifndef SPINOR_FUNCTION_H
#define SPINOR_FUNCTION_H

#include "Lattice.h"
#include "Enums.h"
#include <memory>
#include <stdio.h>
#include <vector>

class RadialFunction;

/** SpinorFunction is the radial part of a Dirac function of the form
    \f[ \psi = {1 \over r} \left({ f(r) \Omega_{\kappa m} \atop
                                  ig(r) \Omega_{-\kappa m} }\right) \f]
    where \f$ \Omega_{\kappa m} \f$ contains the angular and spin variables.
 */
class SpinorFunction : public std::enable_shared_from_this<SpinorFunction>
{
public:
    SpinorFunction(int kappa, unsigned int size = 0);
    virtual ~SpinorFunction() = default;

    std::vector<double> f, g,       // Upper and lower components of the radial functions
                        dfdr, dgdr; // Function derivatives, ie. df/dr, dg/dr

    virtual unsigned int size() const { return static_cast<unsigned int>(f.size()); }
    int Kappa() const { return kappa; }
    int L() const;
    double J() const;
    int TwoJ() const;
    inline Parity GetParity() const;

    /** Resize the functions. Pad with zeros if necessary. */
    virtual void resize(unsigned int size);
    virtual void Clear();

    virtual void SetKappa(int new_kappa) { kappa = new_kappa; }
    virtual void swap(SpinorFunction& other);

    /** Multiply all points of all vectors (f, g, df, dg) by the scale factor. */
    const SpinorFunction& operator*=(double scale_factor);
    SpinorFunction operator*(double scale_factor) const;

    /** Adding or subtracting two spinor functions can only occur if both have same angular part. */
    const SpinorFunction& operator+=(const SpinorFunction& other);
    const SpinorFunction& operator-=(const SpinorFunction& other);
    SpinorFunction operator+(const SpinorFunction& other) const;
    SpinorFunction operator-(const SpinorFunction& other) const;

    /** Multiply spinor function by another radial function (assumed zero outside range). */
    const SpinorFunction& operator*=(const RadialFunction& chi);
    SpinorFunction operator*(const RadialFunction& chi) const;

    /** Get density
            rho(r) = (f^2 + g^2)
     */
    RadialFunction GetDensity() const;

    /** Get density
            rho(r) = (f_i f_j + g_i g_j)
     */
    RadialFunction GetDensity(const SpinorFunction& other) const;

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

typedef std::shared_ptr<SpinorFunction> pSpinorFunction;
typedef std::shared_ptr<const SpinorFunction> pSpinorFunctionConst;

/** RadialFunction is a vector function and its derivative,
    similar to one component of a SpinorFunction.
 */
class RadialFunction
{
public:
    RadialFunction(const std::vector<double>& pf, const std::vector<double>& pdfdr);
    RadialFunction(unsigned int size = 0);
    virtual ~RadialFunction() = default;

    std::vector<double> f, dfdr;
    
    virtual unsigned int size() const { return static_cast<unsigned int>(f.size()); }
    /** Resize the functions. Pad with zeros if necessary. */
    virtual void resize(unsigned int size);
    virtual void Clear();

    /** Multiply all points by the scale factor. */
    const RadialFunction& operator*=(double scale_factor);
    RadialFunction operator*(double scale_factor) const;

    const RadialFunction& operator+=(const RadialFunction& other);
    const RadialFunction& operator-=(const RadialFunction& other);
    RadialFunction operator+(const RadialFunction& other) const;
    RadialFunction operator-(const RadialFunction& other) const;
    
    /** Multiply radial function by another (assumed zero outside range). */
    const RadialFunction& operator*=(const RadialFunction& other);
    RadialFunction operator*(const RadialFunction& other) const;

    /** Print state to file, optionally printing lattice. Return success. */
    bool Print(const std::string& filename, pLattice lattice = pLattice()) const;
    bool Print(FILE* fp, pLattice lattice = pLattice()) const;
};

inline int SpinorFunction::L() const
{   if (kappa > 0)
    return (kappa);
else
    return (-kappa-1);
}

inline double SpinorFunction::J() const
{   return (double)abs(kappa) - 0.5;
}

inline int SpinorFunction::TwoJ() const
{   return (2*abs(kappa) - 1);
}

inline Parity SpinorFunction::GetParity() const
{   if(L()%2 == 0)
    return Parity::even;
else
    return Parity::odd;
}

#endif
