#ifndef MATH_CONSTANT_H
#define MATH_CONSTANT_H

#include <map>

/** Class of mathematical constants, following the Singleton pattern. 
    Note that AMBiT uses atomic units, hbar = m_e = e = 1, so
    Bohr radius = 1, energy of ground state in hydrogen is 1/2.
 */
class MathConstant
{
public:
    static MathConstant* Instance();

    // Mathematical constants
    double Pi() const;
    double Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3);
    double Wigner6j(double j1, double j2, double j3, double j4, double j5, double j6);
    
    /** Calculate 3j symbol where j1, j2 are half integer and k is integer.
        Assumes q (projection of k) = - m1 - m2.
        ( j1  j2  k )
        ( m1  m2  q )
     */
    double Electron3j(unsigned int twoj1, unsigned int twoj2, unsigned int k, int twom1, int twom2);

    /** Calculate 3j symbol
        ( j1   j2   k )
        ( 1/2 -1/2  0 )
        where j1 and j2 are half integer and k is integer.
     */
    double Electron3j(unsigned int twoj1, unsigned int twoj2, unsigned int k);

    /** Get number of stored Electron3j symbols. */
    unsigned int GetStorageSize() const;

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

protected:
    MathConstant();
    ~MathConstant();

protected:
    char* SpectroscopicNotation;
    std::map<int, double> Symbols3j;

    unsigned int MaxStoredTwoJ;
    unsigned int MSize;
    std::size_t HashElectron3j(unsigned int twoj1, unsigned int twoj2, unsigned int k, int twom1, int twom2);
    
    /** Calculate the logarithm of a fraction where the numerator and denominator are factorials
     log( n!/d! )
     */
    double LogFactorialFraction(unsigned int num, unsigned int denom);
};

#endif
