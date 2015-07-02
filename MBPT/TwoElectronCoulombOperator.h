#ifndef TWO_ELECTRON_COULOMB_OPERATOR_H
#define TWO_ELECTRON_COULOMB_OPERATOR_H

#include "Include.h"
#include "HartreeFock/SpinorOperator.h"
#include "HartreeFock/HartreeY.h"
#include "Configuration/ElectronInfo.h"
#include "Universal/MathConstant.h"
#include "MBPT/SlaterIntegrals.h"

//#define INCLUDE_EXTRA_BOX_DIAGRAMS 1

template <class pTwoElectronIntegralType>
class TwoElectronCoulombOperator
{
public:
    TwoElectronCoulombOperator(pTwoElectronIntegralType ci_integrals): integrals(ci_integrals) {}

    double GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4) const;

    /** < e1, e2 || V^k || e3, e4 > = [j1, j2, j3, j4]^(1/2) (k  j1  j3 ) (k  j2  j4 ) R^k(12,34)
                                                             (0 -1/2 1/2) (0 -1/2 1/2)
            * \xi(l1 + l3 + k) \xi(l2 + l4 + k)
     */
    double GetReducedMatrixElement(int k, const OrbitalInfo& e1, const OrbitalInfo& e2, const OrbitalInfo& e3, const OrbitalInfo& e4) const;

protected:
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

            coeff = coeff * sqrt(double(e1.MaxNumElectrons() * e2.MaxNumElectrons() *
                                        e3.MaxNumElectrons() * e4.MaxNumElectrons()));

            double radial = integrals->GetTwoElectronIntegral(k, e1, e2, e3, e4);

            total += coeff * radial;
        }

        k = k+2;
    }

#ifdef INCLUDE_EXTRA_BOX_DIAGRAMS
    // Include the box diagrams with "wrong" parity.
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
                if(int(q - e1.M() - e2.M() + 1.)%2)
                    coeff = - coeff;

                coeff = coeff * sqrt(double(e1.MaxNumElectrons() * e2.MaxNumElectrons() *
                                            e3.MaxNumElectrons() * e4.MaxNumElectrons()));
                
                total += coeff * radial;
            }
        }
        
        k = k+2;
    }
#endif
    
    return total;
}

template <class TwoElectronIntegralType>
double TwoElectronCoulombOperator<TwoElectronIntegralType>::GetReducedMatrixElement(int k, const OrbitalInfo& e1, const OrbitalInfo& e2, const OrbitalInfo& e3, const OrbitalInfo& e4) const
{
    MathConstant* math = MathConstant::Instance();

    if(!math->sum_is_even(e1.L(), e3.L(), k) || !math->sum_is_even(e2.L(), e4.L(), k))
        return 0.0;

    double total = math->Electron3j(e1.TwoJ(), e3.TwoJ(), k) *
            math->Electron3j(e2.TwoJ(), e4.TwoJ(), k);

    if(total)
    {
        total *= sqrt(double(e1.MaxNumElectrons() * e2.MaxNumElectrons() *
                             e3.MaxNumElectrons() * e4.MaxNumElectrons()));
        total *= integrals->GetTwoElectronIntegral(k, e1, e2, e3, e4);
    }

    return total;
}

#endif
