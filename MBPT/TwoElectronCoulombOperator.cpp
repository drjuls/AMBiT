#include "TwoElectronCoulombOperator.h"
#include "HartreeFock/SpinorOperator.h"
#include "HartreeFock/HartreeY.h"
#include "Configuration/ElectronInfo.h"
#include "Universal/MathConstant.h"

namespace Ambit
{
double TwoElectronCoulombOperator::GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4) const
{
    if((e1.L() + e2.L() + e3.L() + e4.L())%2)
        return 0.;

    int two_q = e1.TwoM() - e3.TwoM();
    if(two_q != - e2.TwoM() + e4.TwoM())
        return 0.;

    int kmin = mmax(abs(e1.TwoJ() - e3.TwoJ()), abs(e2.TwoJ() - e4.TwoJ()))/2;
    int k = kmin;
    if((e1.L() + e3.L() + k)%2)
        k++;

    int kmax = mmin(e1.TwoJ() + e3.TwoJ(), e2.TwoJ() + e4.TwoJ())/2;
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

        k += 2;
    }

    // Include the box diagrams with "wrong" parity.
    if(include_off_parity)
    {
        k = kmin;
        if((e1.L() + e3.L() + k)%2 == 0)
            k++;

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

            k += 2;
        }
    }

    return total;
}

double TwoElectronCoulombOperator::GetReducedMatrixElement(int k, const OrbitalInfo& e1, const OrbitalInfo& e2, const OrbitalInfo& e3, const OrbitalInfo& e4) const
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
}
