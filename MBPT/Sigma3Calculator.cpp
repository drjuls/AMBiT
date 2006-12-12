#include "Sigma3Calculator.h"
#include "Include.h"
#include "Universal/Constant.h"
#include "Universal/CoulombIntegrator.h"
#include "HartreeFock/StateIntegrator.h"

Sigma3Calculator::Sigma3Calculator(Lattice* lattice, const Core* atom_core, const ExcitedStates* excited_states):
    MBPTCalculator(lattice, atom_core, excited_states)
{
    UseBrillouinWignerPT();
}

double Sigma3Calculator::GetSecondOrderSigma3(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
           const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const
{
    const double NuclearInverseMass = core->GetNuclearInverseMass();

    std::vector<double> density(MaxStateSize);
    std::vector<double> Pot14(MaxStateSize);
    std::vector<double> Pot36(MaxStateSize);
    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);
    const double* dR = lattice->dR();

    const DiscreteState& s1(*excited->GetState(e1));
    const DiscreteState& s2(*excited->GetState(e2));
    const DiscreteState& s3(*excited->GetState(e3));
    const DiscreteState& s4(*excited->GetState(e4));
    const DiscreteState& s5(*excited->GetState(e5));
    const DiscreteState& s6(*excited->GetState(e6));

    // k1 limits
    int q1 = (e1.TwoM() - e4.TwoM())/2;
    unsigned int k1min = abs(int(e1.L()) - int(e4.L()));
    if(fabs(e1.J() - e4.J()) > double(k1min))
        k1min += 2;
    while(k1min < abs(q1))
        k1min += 2;
    unsigned int k1max = (e1.TwoJ() + e4.TwoJ())/2;

    // k2 limits
    int q2 = (e6.TwoM() - e3.TwoM())/2;
    unsigned int k2min = abs(int(e3.L()) - int(e6.L()));
    if(fabs(e3.J() - e6.J()) > double(k2min))
        k2min += 2;
    while(k2min < abs(q2))
        k2min += 2;
    unsigned int k2max = (e3.TwoJ() + e6.TwoJ())/2;

    // core state limits
    int two_mn = e1.TwoM() + e2.TwoM() - e4.TwoM();
    double mn = double(two_mn)/2.;
    int ln_parity = (k1min + e2.L())%2;

    const double ValenceEnergy = ValenceEnergies.find(e2.Kappa())->second;

    double total = 0.;
    unsigned int k1, k2;
    unsigned int i;

    // Summation over k1
    for(k1 = k1min; k1 <= k1max; k1+=2)
    {
        double coeff14 =  Constant::Electron3j(e1.TwoJ(), e4.TwoJ(), k1, 1, -1)
                       * Constant::Electron3j(e1.TwoJ(), e4.TwoJ(), k1, -e1.TwoM(), e4.TwoM());

        if(coeff14)
        {   for(i=0; i<mmin(s1.Size(), s4.Size()); i++)
            {
                density[i] = s1.f[i] * s4.f[i] + Constant::AlphaSquared * s1.g[i] * s4.g[i];
            }
            I.FastCoulombIntegrate(density, Pot14, k1, mmin(s2.Size(), s4.Size()));

            double SMS_14 = 0.;
            if(NuclearInverseMass && (k1 == 1))
                SMS_14 = -SI.IsotopeShiftIntegral(s4, s1);

            // Summation over k2
            for(k2 = k2min; k2 <= k2max; k2+=2)
            {
                double coeff36 =  Constant::Electron3j(e3.TwoJ(), e6.TwoJ(), k2, 1, -1)
                                * Constant::Electron3j(e3.TwoJ(), e6.TwoJ(), k2, -e3.TwoM(), e6.TwoM());

                if(coeff36)
                {   for(i=0; i<mmin(s3.Size(), s6.Size()); i++)
                    {
                        density[i] = s3.f[i] * s6.f[i] + Constant::AlphaSquared * s3.g[i] * s6.g[i];
                    }
                    I.FastCoulombIntegrate(density, Pot36, k2, mmin(s3.Size(), s6.Size()));

                    double SMS_36 = 0.;
                    if(NuclearInverseMass && (k2 == 1))
                        SMS_36 = -SI.IsotopeShiftIntegral(s6, s3);

                    // Summation over core state 'n'
                    ConstStateIterator itn = core->GetConstStateIterator();
                    while(!itn.AtEnd())
                    {   const DiscreteState& sn = *itn.GetState();
                        unsigned int two_Jn = (unsigned int)(sn.J() * 2.);

                        if((two_Jn >= abs(two_mn)) && (sn.L()%2 == ln_parity))
                        {
                            ElectronInfo en(sn.RequiredPQN(), sn.Kappa(), two_mn);

                            double coeff = double(en.MaxNumElectrons())/(sn.Energy() - ValenceEnergy + delta)
                                * Constant::Electron3j(e2.TwoJ(), two_Jn, k1, 1, -1)
                                * Constant::Electron3j(e2.TwoJ(), two_Jn, k1, -e2.TwoM(), two_mn)
                                * Constant::Electron3j(two_Jn, e5.TwoJ(), k2, 1, -1)
                                * Constant::Electron3j(two_Jn, e5.TwoJ(), k2, -two_mn, e5.TwoM())
                                * coeff14 * coeff36;
                            
                            if(coeff)
                            {
                                // R1 = R_k1 (12, 4n) = R_k1 (21, n4)
                                double R1 = 0.;
                                for(i=0; i < mmin(s2.Size(), sn.Size()); i++)
                                    R1 = R1 + Pot14[i] * (s2.f[i] * sn.f[i] + Constant::AlphaSquared * s2.g[i] * sn.g[i]) * dR[i];

                                if(SMS_14)
                                {   R1 = R1 - NuclearInverseMass * SMS_14 * SI.IsotopeShiftIntegral(s2, sn);
                                }

                                // R2 = R_k2 (n3, 56)
                                double R2 = 0.;
                                for(i=0; i < mmin(sn.Size(), s5.Size()); i++)
                                    R2 = R2 + Pot36[i] * (sn.f[i] * s5.f[i] + Constant::AlphaSquared * sn.g[i] * s5.g[i]) * dR[i];

                                if(SMS_36)
                                {   R2 = R2 + NuclearInverseMass * SMS_36 * SI.IsotopeShiftIntegral(s5, sn);
                                }

                                total += coeff * R1 * R2;
                            }
                        }
                        itn.Next();
                    }
                }
            }
        }
    }

    total = total * sqrt(double(e1.MaxNumElectrons() * e2.MaxNumElectrons() *
                                e3.MaxNumElectrons() * e4.MaxNumElectrons() *
                                e5.MaxNumElectrons() * e6.MaxNumElectrons()));

    if(((e2.TwoM() + e3.TwoM() + e4.TwoM() + e5.TwoM())/2)%2 == 1)
        total = - total;

    return total;
}
