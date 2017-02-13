#ifndef TWO_ELECTRON_COULOMB_OPERATOR_H
#define TWO_ELECTRON_COULOMB_OPERATOR_H

#include "Configuration/ElectronInfo.h"
#include "MBPT/SlaterIntegrals.h"

/** Holds two-electron radial integrals (which may have MBPT) and adds angular part to give two-body matrix elements.
    Option include_off_parity will include "off-parity" matrix elements if they are found in the radial integrals
    (these diagams are found in MBPT and Breit interactions).
 */
class TwoElectronCoulombOperator
{
public:
    TwoElectronCoulombOperator(pSlaterIntegrals ci_integrals, bool include_off_parity):
        integrals(ci_integrals), include_off_parity(include_off_parity)
    {}

    TwoElectronCoulombOperator(pSlaterIntegrals ci_integrals):
        integrals(ci_integrals)
    {   include_off_parity = ci_integrals->OffParityExists();
    }

    double GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4) const;

    /** < e1, e2 || V^k || e3, e4 > = [j1, j2, j3, j4]^(1/2) (k  j1  j3 ) (k  j2  j4 ) R^k(12,34)
                                                             (0 -1/2 1/2) (0 -1/2 1/2)
            * \xi(l1 + l3 + k) \xi(l2 + l4 + k)
     */
    double GetReducedMatrixElement(int k, const OrbitalInfo& e1, const OrbitalInfo& e2, const OrbitalInfo& e3, const OrbitalInfo& e4) const;

protected:
    bool include_off_parity;
    pSlaterIntegrals integrals;
};

typedef std::shared_ptr<TwoElectronCoulombOperator> pTwoElectronCoulombOperator;

#endif
