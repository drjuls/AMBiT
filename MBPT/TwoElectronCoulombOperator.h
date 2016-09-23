#ifndef TWO_ELECTRON_COULOMB_OPERATOR_H
#define TWO_ELECTRON_COULOMB_OPERATOR_H

#include "Include.h"
#include "HartreeFock/SpinorOperator.h"
#include "HartreeFock/HartreeY.h"
#include "Configuration/ElectronInfo.h"
#include "Universal/MathConstant.h"
#include "MBPT/SlaterIntegrals.h"

/** Holds two-electron radial integrals (which may have MBPT) and adds angular part to give two-body matrix elements.
    Option include_off_parity will include "off-parity" matrix elements if they are found in the radial integrals
    (these diagams are found in MBPT and Breit interactions).
 */
template <class pTwoElectronIntegralType>
class TwoElectronCoulombOperator
{
public:
    TwoElectronCoulombOperator(pTwoElectronIntegralType ci_integrals, bool include_off_parity):
        integrals(ci_integrals), include_off_parity(include_off_parity)
    {}

    TwoElectronCoulombOperator(pTwoElectronIntegralType ci_integrals):
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
    pTwoElectronIntegralType integrals;
};

typedef std::shared_ptr<TwoElectronCoulombOperator<pSlaterIntegrals>> pTwoElectronCoulombOperator;

template <class TwoElectronIntegralType>
double TwoElectronCoulombOperator<TwoElectronIntegralType>::GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4) const
{
    if((e1.L() + e2.L() + e3.L() + e4.L())%2)
        return 0.;

    int two_q = e1.TwoM() - e3.TwoM();
    if(two_q != - e2.TwoM() + e4.TwoM())
        return 0.;

    int k = mmax(abs(e1.L() - e3.L()), abs(e2.L() - e4.L()));
    if((fabs(e1.J() - e3.J()) > double(k)) || (fabs(e2.J() - e4.J()) > double(k)))
        k += 2;

    int kmax = mmin(e1.L() + e3.L(), e2.L() + e4.L());
    if((e1.J() + e3.J() < double(kmax)) || (e2.J() + e4.J() < double(kmax)))
        kmax -= 2;

    double q = double(two_q)/2.;

    double total = 0.;
    double sqrt_multiplicity = sqrt(double(e1.MaxNumElectrons() * e2.MaxNumElectrons() *
                                           e3.MaxNumElectrons() * e4.MaxNumElectrons()));

    MathConstant* constants = MathConstant::Instance();

    while(k <= kmax)
    {
        double coeff = 0.;
        if(fabs(q) <= k)
            coeff = constants->Electron3j(e1.TwoJ(), e3.TwoJ(), k, -e1.TwoM(), e3.TwoM()) *
                    constants->Electron3j(e2.TwoJ(), e4.TwoJ(), k, -e2.TwoM(), e4.TwoM());

        if(coeff)
            coeff = coeff * constants->Electron3j(e1.TwoJ(), e3.TwoJ(), k, 1, -1) *
                    constants->Electron3j(e2.TwoJ(), e4.TwoJ(), k, 1, -1);

        if(coeff)
        {
            if(((two_q - e1.TwoM() - e2.TwoM())/2 + 1)%2)
                coeff = - coeff;

            double radial = integrals->GetTwoElectronIntegral(k, e1, e2, e3, e4);

            total += coeff * radial * sqrt_multiplicity;
        }

        k = k+2;
    }

    // Include the box diagrams with "wrong" parity.
    if(include_off_parity)
    {
        k = mmax(abs(e1.L() - e3.L()), abs(e2.L() - e4.L())) + 1;
        kmax = mmin(e1.L() + e3.L(), e2.L() + e4.L()) - 1;

        while(k <= kmax)
        {
            double radial = integrals->GetTwoElectronIntegral(k, e1, e2, e3, e4);

            if(radial)
            {
                double coeff = 0.;
                if(fabs(q) <= k)
                    coeff = constants->Electron3j(e1.TwoJ(), e3.TwoJ(), k, -e1.TwoM(), e3.TwoM()) *
                            constants->Electron3j(e2.TwoJ(), e4.TwoJ(), k, -e2.TwoM(), e4.TwoM());

                if(coeff)
                    coeff = coeff * constants->Electron3j(e1.TwoJ(), e3.TwoJ(), k, 1, -1) *
                            constants->Electron3j(e2.TwoJ(), e4.TwoJ(), k, 1, -1);

                if(coeff)
                {
                    if(((two_q - e1.TwoM() - e2.TwoM())/2 + 1)%2)
                        coeff = - coeff;

                    total += coeff * radial * sqrt_multiplicity;
                }
            }

            k = k+2;
        }
    }

    return total;
}

template <class TwoElectronIntegralType>
double TwoElectronCoulombOperator<TwoElectronIntegralType>::GetReducedMatrixElement(int k, const OrbitalInfo& e1, const OrbitalInfo& e2, const OrbitalInfo& e3, const OrbitalInfo& e4) const
{
    MathConstant* math = MathConstant::Instance();

    if(!include_off_parity)
    {
        if(!math->sum_is_even(e1.L(), e3.L(), k) || !math->sum_is_even(e2.L(), e4.L(), k))
            return 0.0;
    }

    double total = integrals->GetTwoElectronIntegral(k, e1, e2, e3, e4);

    if(total)
    {
        total *= math->Electron3j(e1.TwoJ(), e3.TwoJ(), k) *
                 math->Electron3j(e2.TwoJ(), e4.TwoJ(), k);

        if(total)
            total *= sqrt(double(e1.MaxNumElectrons() * e2.MaxNumElectrons() *
                                 e3.MaxNumElectrons() * e4.MaxNumElectrons()));
    }

    return total;
}

#endif
