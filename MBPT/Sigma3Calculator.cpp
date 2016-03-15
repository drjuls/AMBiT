#include "Sigma3Calculator.h"
#include "Include.h"
#include "Universal/PhysicalConstant.h"

Sigma3Calculator::Sigma3Calculator(pOrbitalManagerConst orbitals, pSlaterIntegrals two_body, bool include_valence):
    MBPTCalculator(orbitals), two_body(two_body), include_valence(include_valence),
    deep(orbitals->deep), high(orbitals->high)
{}

Sigma3Calculator::~Sigma3Calculator(void)
{
    two_body->clear();
}

unsigned int Sigma3Calculator::GetStorageSize()
{
    unsigned int total = two_body->CalculateTwoElectronIntegrals(deep, valence, valence, valence, true);
    if(include_valence)
        total += two_body->CalculateTwoElectronIntegrals(valence, valence, valence, high, true);

    return total;
}

void Sigma3Calculator::UpdateIntegrals()
{
    SetValenceEnergies();

    two_body->CalculateTwoElectronIntegrals(deep, valence, valence, valence);
    if(include_valence)
        two_body->CalculateTwoElectronIntegrals(valence, valence, valence, high);
}

double Sigma3Calculator::GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
                                          const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const
{
    double value = 0.;

    if((e1.L() + e2.L() + e3.L() + e4.L() + e5.L() + e6.L())%2 ||
       (e1.TwoM() + e2.TwoM() + e3.TwoM() != e4.TwoM() + e5.TwoM() + e6.TwoM()))
        return value;

    // There are no sign changes, since there are the same number
    // of permutations on both sides

    value = GetSecondOrderSigma3(e1, e2, e3, e4, e5, e6) +
            GetSecondOrderSigma3(e2, e3, e1, e5, e6, e4) +
            GetSecondOrderSigma3(e3, e1, e2, e6, e4, e5) +
            GetSecondOrderSigma3(e3, e2, e1, e6, e5, e4) +
            GetSecondOrderSigma3(e2, e1, e3, e5, e4, e6) +
            GetSecondOrderSigma3(e1, e3, e2, e4, e6, e5);

    return value;
}

double Sigma3Calculator::GetSecondOrderSigma3(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
           const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const
{
    // core state limits
    int two_mn = e1.TwoM() + e2.TwoM() - e4.TwoM();

    // k1 limits
    int q1 = (e1.TwoM() - e4.TwoM())/2;
    int k1min = kmin(e1, e4);
    while(k1min < abs(q1))
        k1min += 2;
    int k1max = kmax(e1, e4);

    // k2 limits
    int q2 = (e6.TwoM() - e3.TwoM())/2;
    int k2min = kmin(e3, e6);
    while(k2min < abs(q2))
        k2min += 2;
    int k2max = kmax(e3, e6);

    if((k1min > k1max) || (k2min > k2max))
        return 0.;

    int ln_parity = (k1min + e2.L())%2;
    const double ValenceEnergy = ValenceEnergies.find(e2.Kappa())->second;

    double total = 0.;
    int k1, k2;
    MathConstant* constants = MathConstant::Instance();

    // Summation over k1
    for(k1 = k1min; k1 <= k1max; k1+=2)
    {
        double coeff14 = constants->Electron3j(e1.TwoJ(), e4.TwoJ(), k1, -e1.TwoM(), e4.TwoM()) *
                         constants->Electron3j(e1.TwoJ(), e4.TwoJ(), k1);

        if(coeff14)
        {
            // Summation over k2
            for(k2 = k2min; k2 <= k2max; k2+=2)
            {
                double coeff36 = constants->Electron3j(e3.TwoJ(), e6.TwoJ(), k2, -e3.TwoM(), e6.TwoM()) *
                                 constants->Electron3j(e3.TwoJ(), e6.TwoJ(), k2);

                if(coeff36)
                {
                    // Summation over core state 'n'
                    auto it_n = deep->begin();
                    while(it_n != deep->end())
                    {
                        const OrbitalInfo& sn = it_n->first;
                        int two_Jn = sn.TwoJ();

                        if((two_Jn >= abs(two_mn)) && (sn.L()%2 == ln_parity))
                        {
                            const double En = it_n->second->Energy();

                            double coeff = constants->Electron3j(e2.TwoJ(), two_Jn, k1, -e2.TwoM(), two_mn) *
                                           constants->Electron3j(two_Jn, e5.TwoJ(), k2, -two_mn, e5.TwoM());

                            if(coeff)
                            {
                                coeff *= constants->Electron3j(e2.TwoJ(), two_Jn, k1) *
                                         constants->Electron3j(two_Jn, e5.TwoJ(), k2) *
                                         coeff14 * coeff36 /(En - ValenceEnergy + delta);

                                coeff *= sqrt(e1.MaxNumElectrons() * e2.MaxNumElectrons() * e3.MaxNumElectrons() *
                                              e4.MaxNumElectrons() * e5.MaxNumElectrons() * e6.MaxNumElectrons())
                                         * double(two_Jn + 1);

                                double R1 = two_body->GetTwoElectronIntegral(k1, sn, e4, e2, e1);
                                double R2 = two_body->GetTwoElectronIntegral(k2, sn, e3, e5, e6);

                                total += coeff * R1 * R2;
                            }
                        }
                        ++it_n;
                    }
                }
            }
        }
    }

    // phase = -1 * (-1)^(m2 + m3 + m4 + m5) = -1 * (-1)^(m1 - m6)
    if((abs(e1.TwoM() - e6.TwoM())/2)%2 == 0)
        total = - total;

    return total;
}
