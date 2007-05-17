#ifndef RATE_CALCULATOR_H
#define RATE_CALCULATOR_H

#include "Atom.h"

class RateCalculator
{
public:
    RateCalculator(ExcitedStates* basis): excited(basis) {}
    virtual ~RateCalculator() {}

    void CalculateAllDipoleStrengths(Atom* A, Symmetry sym1, unsigned int solution1);
    double CalculateDipoleStrength(Atom* A, Symmetry sym1, unsigned int solution1, Symmetry sym2, unsigned int solution2);

    double CalculateAugerRate(Atom* ion, Symmetry gs, Atom* atom, Symmetry sym1, unsigned int solution1);
    //void UpdateOneElectronIntegrals();

protected:
    ExcitedStates* excited;

    // E1Integrals(i, j) = <i|r|j>
    //std::map<unsigned int, double> E1Integrals;

    double GetE1MatrixElement(const ElectronInfo& e1, const ElectronInfo& e2) const;

    /** Get the Hamiltonian matrix element <first + continuum | H | second> */
    double GetProjectionH(const Projection& first, const Projection& second, const ContinuumState* cs, const ElectronInfo& cs_electron) const;

    /** Get the Coulomb matrix element < e1, e2 | 1/r | e3, e4 >. */
    double CoulombMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4, const ContinuumState* cs) const;
};

#endif
