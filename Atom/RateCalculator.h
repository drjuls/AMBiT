#ifndef RATE_CALCULATOR_H
#define RATE_CALCULATOR_H

#include <complex>
#include "Atom.h"
#include "Universal/Enums.h"

class Atom;
class Symmetry;

class RateCalculator
{
public:
    RateCalculator(pExcitedStates basis);
    virtual ~RateCalculator() {}

    /** Calculate E1 matrix elements to all states that lie below solution1 of sym1.
        Optionally print radiative rates (inverse seconds) and/or generalised oscillator strengths.
        Return total radiative rate (sum of all rates, in inverse seconds).
     */
    double CalculateAllDipoleStrengths(Atom* A, Symmetry sym1, unsigned int solution1, bool print_rates = true, bool print_oscillator_strengths = false);
    double CalculateDipoleStrength(Atom* A, Symmetry sym1, unsigned int solution1, Symmetry sym2, unsigned int solution2);
    double CalculateMultipoleStrength(MultipolarityType::Enum MType, unsigned int J, Atom* A, Symmetry sym1, unsigned int solution1, Symmetry sym2, unsigned int solution2, TransitionGaugeType::Enum aGaugeType = TransitionGaugeType::Length, TransitionCalculationMethod::Enum aMethod = TransitionCalculationMethod::WalterJohnson);

    /** Calculate Auger rates of num_solutions levels above the continuum energy.
        The continuum energy will be the energy of the atom ground state plus the ionization potential.
     */
    double CalculateAugerRate(Atom* atom, Symmetry sym1, unsigned int solution1, double continuum_energy);

    /** Calculate Auger rates and radiative rates for states above the continuum;
        from them calculate the DR cross-section and print the results to a file.
        This function is made to be modified.
     */
    void DielectronicRecombination(Atom* A);
    void AutoionisationRates(Atom* A);

protected:
    pExcitedStates excited;

    unsigned int NumStates;
    std::map<OrbitalInfo, unsigned int> state_index;

    // E1Integrals(i, j) = <i|r|j>
    std::map<unsigned int, double> E1Integrals;

    double GetE1MatrixElement(const ElectronInfo& e1, const ElectronInfo& e2);
    double GetM1MatrixElement(const ElectronInfo& e1, const ElectronInfo& e2);

    /** Get the Hamiltonian matrix element <first + continuum | H | second> */
    double GetProjectionH(const Projection& first, const Projection& second, pContinuumWaveConst cs, const ElectronInfo& cs_electron) const;

    /** Get the Coulomb matrix element < e1, e2 | 1/r | e3, e4 >.
        "sign" provides an additional phase which doesn't affect the results but is multiplied
        by the angular factor when printing debug info.
     */
    double CoulombMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4, pOrbitalConst cs, int sign = 1) const;
    double SubtractionDiagram(pContinuumWaveConst sa, pOrbitalConst sb, pOrbitalConst sc, pOrbitalConst sd, unsigned int k) const;
    
    double GetEJMatrixElement(unsigned int J, const ElectronInfo& e1, const ElectronInfo& e2, TransitionGaugeType::Enum aGaugeType = TransitionGaugeType::Length);
    std::complex<double> GetEJMatrixElementZenonas(unsigned int J, const ElectronInfo& e1, const ElectronInfo& e2, TransitionGaugeType::Enum aGaugeType = TransitionGaugeType::Length);
    double GetMJMatrixElement(Atom* A, unsigned int J, const ElectronInfo& e1, const ElectronInfo& e2);
};

#endif
