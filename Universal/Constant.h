#ifndef CONSTANT_H
#define CONSTANT_H

#include <map>

class Constant
{
public:
    Constant(void) {}
    ~Constant(void) {}

    // Natural variables
    const static double Pi;
    static double Alpha;
    static double AlphaSquared;
    static double Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3);
    static double Wigner6j(double j1, double j2, double j3, double j4, double j5, double j6);
    
    /** Calculate 3j symbol where j1, j2 are half integer and k is integer.
        Assumes q (projection of k) = - m1 - m2.
     */
    static double Electron3j(unsigned int twoj1, unsigned int twoj2, unsigned int k, int twom1, int twom2);

    const static double NucleonElectronMassRatio;

    // Conversion factors
    const static double AtomicToFermi;
    const static double HartreeEnergy_eV;
    const static double HartreeEnergy_cm;
    const static double SpeedOfLight;       // In SI units (m/s)
    const static double InvCmToMHz;         // Multiply by speed of light

    static const char SpectroscopicNotation[10];

private:
    const static unsigned int MaxStoredTwoJ;
    static std::map<int, double> Symbols3j;
};

#endif