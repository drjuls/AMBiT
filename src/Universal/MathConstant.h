#ifndef MATH_CONSTANT_H
#define MATH_CONSTANT_H

#include <map>
#include <stdlib.h>
#include <memory>
#include <boost/math/special_functions.hpp>
#include <absl/container/flat_hash_map.h>
#include <gsl/gsl_math.h>
#include "Enums.h"
#include "Include.h"

#define BIGGEST_K_EVER 100

namespace Ambit
{
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
    double SpeedOfLightAU() const;
    double InvCmInMHz() const;            //!< Multiply by speed of light (cm/10^6 s).
    double BohrRadiusSI() const;          //!< Atomic unit in meters.
    double RydbergConstantSI() const;     //!< In SI units (inverse m).
    double AtomicFrequencySI() const;     //!< Multiply a frequency in atomic units by this factor to get frequency in SI units (Hz).
    double AtomicFrequencyMHz() const;    //!< Atomic frequency in MHz.
    double NuclearMagneton() const;       //!< Nuclear magneton 1/(2 m_p c) in atomic units.
    double Barn() const;                  //!< 1 barn = 10^(-24) cm^2 in atomic units (square Bohr radius)
    constexpr static double ElectronMassInEV {510998.9461};

    char GetSpectroscopicNotation(unsigned int l) const;    //!< Returns 0 if not found
    int GetL(char spectroscopic_notation) const;            //!< Returns -1 if not found

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

    /** Convert half-integer J and parity to kappa = (-1)^(l+j+1/2) (j+1/2) */
    static inline int convert_to_kappa(int twoj, Parity P);

    /** Find minimum k that satisfies triangle(a.J(), b.J(), kmin) and
     * (a.L() + b.L() + kmin) is even unless include_off_parity is set.
     * Types can be anything with TwoJ(), L(), GetParity() defined. */
    template <typename A, typename B>
    static int kmin(const A& a, const B& b, bool include_off_parity = false);

    /** Find minimum k connecting a -> c and b -> d. */
    template <typename A, typename B, typename C, typename D>
    static int kmin(const A& a, const B& b, const C& c, const D& d, bool include_off_parity = false);

    /** Get maximum k. */
    template <typename A, typename B>
    static int kmax(const A& a, const B& b, bool include_off_parity = false);

    template <typename A, typename B, typename C, typename D>
    static int kmax(const A& a, const B& b, const C& c, const D& d, bool include_off_parity = false);

    template <typename A, typename B>
    static std::pair<int, int> klimits(const A& a, const B& b, bool include_off_parity = false);

    template <typename A, typename B, typename C, typename D>
    static std::pair<int, int> klimits(const A& a, const B& b, const C& c, const D& d, bool include_off_parity = false);

    /** Return +1 if power is even, -1 if odd. */
    inline int minus_one_to_the_power(int power) const { return (abs(power)%2)? -1 : 1; }

    /** Return true if triangular condition is satisfied, false otherwise. */
    inline bool triangular_condition(int a, int b, int c) const { return (abs(a-b) <= c) && ((a+b) >= c); }

    /** Return true if (a + b + c) is even. */
    inline bool sum_is_even(int a, int b, int c) const { return (a + b + c)%2 == 0; }
    inline bool sum_is_even(int a, int b, int c, int d) const { return (a + b + c + d)%2 == 0; }

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
    {   return gsl_pow_int(x, v)/boost::math::double_factorial<double>(2 * v + 1);
    }

    /** Derivative of spherical bessel function j_v(x) in the limit of small x:
            j_v'(x) = v x^(v-1) / (2v + 1)!!
     */
    inline double sph_bessel_small_limit_prime(unsigned int v, double x) const
    {   return v * gsl_pow_int(x, int(v)-1)/boost::math::double_factorial<double>(2 * v + 1);
    }

    /** Return number of combinations of k elements in n slots. */
    unsigned int nChoosek(unsigned int n, unsigned int k) const;

protected:
    MathConstant();

public:
    ~MathConstant();

protected:
    const std::string SpectroscopicNotation;
    absl::flat_hash_map<int, double> Symbols3j;

    unsigned int MaxStoredTwoJ;
    unsigned int MSize;
    int HashWigner3j(int twoj1, int twoj2, int twoj3, int twom1, int twom2) const;
};

int MathConstant::convert_to_kappa(int twoj, Parity P)
{
    int kappa = (twoj + 1)/2;
    int L = kappa;

    // If inconsistent, swap kappa
    if((L%2 == 0 && P == Parity::odd) || (L%2 == 1 && P == Parity::even))
        kappa = -kappa;

    return kappa;
}

template <typename A, typename B>
int MathConstant::kmin(const A& a, const B& b, bool include_off_parity)
{
    // ensure triangle(a.J(), b.J(), kmin)
    int k = abs(a.TwoJ() - b.TwoJ())/2;

    // ensure eps(a.L() + b.L() + kmin)
    if(!include_off_parity && (a.L() + b.L() + k)%2 != 0)
        ++k;

    return k;
}

template <typename A, typename B, typename C, typename D>
int MathConstant::kmin(const A& a, const B& b, const C& c, const D& d, bool include_off_parity)
{
    if((a.L() + b.L() + c.L() + d.L())%2 != 0)
        return BIGGEST_K_EVER;

    int opt1 = kmin(a, c, include_off_parity);
    int opt2 = kmin(b, d, include_off_parity);

    return mmax(opt1, opt2);
}

template <typename A, typename B>
int MathConstant::kmax(const A& a, const B& b, bool include_off_parity)
{
    // ensure triangle(a.J(), b.J(), kmin)
    int k = (a.TwoJ() + b.TwoJ())/2;

    // ensure eps(a.L() + b.L() + kmax)
    if(!include_off_parity && (a.L() + b.L() + k)%2 != 0)
        --k;

    return k;
}

template <typename A, typename B, typename C, typename D>
int MathConstant::kmax(const A& a, const B& b, const C& c, const D& d, bool include_off_parity)
{
    if((a.L() + b.L() + c.L() + d.L())%2 != 0)
        return -1;

    int opt1 = kmax(a, c, include_off_parity);
    int opt2 = kmax(b, d, include_off_parity);

    return mmin(opt1, opt2);
}

template <typename A, typename B>
std::pair<int, int> MathConstant::klimits(const A& a, const B& b, bool include_off_parity)
{
    return std::make_pair(kmin(a, b, include_off_parity), kmax(a, b, include_off_parity));
}

template <typename A, typename B, typename C, typename D>
std::pair<int, int> MathConstant::klimits(const A& a, const B& b, const C& c, const D& d, bool include_off_parity)
{
    return std::make_pair(kmin(a, b, c, d, include_off_parity), kmax(a, b, c, d, include_off_parity));
}

}
#endif
