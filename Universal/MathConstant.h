#ifndef MATH_CONSTANT_H
#define MATH_CONSTANT_H

#include <map>
#include <stdlib.h>
#include <boost/math/special_functions.hpp>

/** Class of mathematical constants, following the Singleton pattern. 
    Note that AMBiT uses atomic units, hbar = m_e = e = 1, so
    Bohr radius = 1, energy of ground state in hydrogen is 1/2.
 */
class MathConstant
{
public:
    static MathConstant* Instance();

    /** Reset and free used memory. Useful for testing. */
    void Reset();

    // Mathematical constants and conversion factors
    double Pi() const;

    // Conversion factors
    double BohrRadiusInFermi() const;     //!< Fermi is 10^(-15) m.
    double HartreeEnergyIneV() const;     //!< Hartree energy natural unit in code, convert to eV.
    double HartreeEnergyInInvCm() const;  //!< Hartree energy natural unit in code, convert to inverse cm (wavenumber).
    double SpeedOfLightSI() const;        //!< Speed of light in SI units (m/s).
    double InvCmInMHz() const;            //!< Multiply by speed of light (cm/10^6 s).
    double BohrRadiusSI() const;          //!< Atomic unit in meters.
    double RydbergConstantSI() const;     //!< In SI units (inverse m).
    double AtomicFrequencySI() const;     //!< Multiply a frequency in atomic units by this factor to get frequency in SI units (Hz).

    char GetSpectroscopicNotation(unsigned int l) const;
    unsigned int GetL(char spectroscopic_notation) const;

    // Functions
    double Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3) const;
    double Wigner6j(double j1, double j2, double j3, double j4, double j5, double j6) const;

    /** Calculate/store/retrieve 3j symbol.
        Assumes m3 = - m1 - m2.
     */
    double Wigner3j(int twoj1, int twoj2, int twoj3, int twom1, int twom2);
    
    /** Calculate 3j symbol where j1, j2 are half integer and k is integer.
        Assumes q (projection of k) = - m1 - m2.
        ( j1  j2  k )
        ( m1  m2  q )
     */
    double Electron3j(int twoj1, int twoj2, int k, int twom1, int twom2);

    /** Calculate 3j symbol
        ( j1   j2   k )
        ( 1/2 -1/2  0 )
        where j1 and j2 are half integer and k is integer.
     */
    double Electron3j(int twoj1, int twoj2, int k);

    /** Get number of stored Electron3j symbols. */
    unsigned int GetStorageSize() const;

    /** Defined by Johnson as
        <kappa_1 || C^k || kappa_2> = (-1)^{j_1 + 1/2} [j_1, j_2]^{1/2} \xi(l_1 + l_2 + k) ( j_1  j_2 k )
                                                                                           ( -1/2 1/2 0 )
        Has symmetry property
        <kappa_1 || C^k || kappa_2> = (-1)^{j_1 + j2 + 1} <kappa_2 || C^k || kappa_1>
     */
    double SphericalTensorReducedMatrixElement(int kappa1, int kappa2, int rank);

    // Helper functions
    inline int is_half_integer(int twoj) const { return twoj%2 == 1; }

    /** Return +1 if power is even, -1 if odd. */
    inline int minus_one_to_the_power(int power) const { return (abs(power)%2)? -1 : 1; }

    /** Return true if triangular condition is satisfied, false otherwise. */
    inline bool triangular_condition(int a, int b, int c) const { return (abs(a-b) <= c) && ((a+b) >= c); }

    /** Return true if (a + b + c) is even. */
    inline bool sum_is_even(int a, int b, int c) const { return (a + b + c)%2 == 0; }

    /** Spherical bessel function j_v(x) */
    inline double sph_bessel(unsigned int v, double x) const { return boost::math::sph_bessel(v, x); }

    /** Derivative of spherical bessel function j_v'(x)
            j_n'(z) = d/dz j_n(z) = (n/z) j_n(z) - j_{n+1}(z)
     */
    inline double sph_bessel_prime(unsigned int v, double x) const
    {   return (double)v/x * boost::math::sph_bessel(v, x) - boost::math::sph_bessel(v + 1, x);
    }

    /** Spherical bessel function j_v(x) in the limit of small x:
            j_v(x) = x^v / (2v + 1)!!
     */
    inline double sph_bessel_small_limit(unsigned int v, double x) const
    {   return pow(x, v)/boost::math::double_factorial<double>(2 * v + 1);
    }

    /** Return number of combinations of k elements in n slots. */
    unsigned int nChoosek(unsigned int n, unsigned int k) const;

protected:
    MathConstant();
    ~MathConstant();

protected:
    char* SpectroscopicNotation;
    std::map<int, double> Symbols3j;

    unsigned int MaxStoredTwoJ;
    unsigned int MSize;
    int HashWigner3j(int twoj1, int twoj2, int twoj3, int twom1, int twom2) const;
    
    /** Calculate the logarithm of a fraction where the numerator and denominator are factorials
     log( n!/d! )
     */
    double LogFactorialFraction(unsigned int num, unsigned int denom) const;
};

#endif
