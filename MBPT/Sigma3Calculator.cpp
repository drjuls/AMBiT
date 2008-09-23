#include "Sigma3Calculator.h"
#include "Include.h"
#include "Universal/Constant.h"
#include "Universal/CoulombIntegrator.h"
#include "HartreeFock/StateIntegrator.h"

inline void swap(unsigned int& i1, unsigned int& i2)
{   unsigned int temp = i1;
    i1 = i2;
    i2 = temp;
}

Sigma3Calculator::Sigma3Calculator(Lattice* lattice, const Core* atom_core, const ExcitedStates* excited_states):
    MBPTCalculator(lattice, atom_core, excited_states)
{
    UpdateStateIndexes();
    Update();
}

/*
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
            I.FastCoulombIntegrate(density, Pot14, k1, mmin(s1.Size(), s4.Size()));

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
*/

double Sigma3Calculator::GetSecondOrderSigma3(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
           const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6)
{
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
    int ln_parity = (k1min + e2.L())%2;

    const double ValenceEnergy = ValenceEnergies.find(e2.Kappa())->second;

    double total = 0.;
    unsigned int k1, k2;

    // Summation over k1
    for(k1 = k1min; k1 <= k1max; k1+=2)
    {
        double coeff14 =  Constant::Electron3j(e1.TwoJ(), e4.TwoJ(), k1, 1, -1)
                       * Constant::Electron3j(e1.TwoJ(), e4.TwoJ(), k1, -e1.TwoM(), e4.TwoM());

        if(coeff14)
        {
            // Summation over k2
            for(k2 = k2min; k2 <= k2max; k2+=2)
            {
                double coeff36 =  Constant::Electron3j(e3.TwoJ(), e6.TwoJ(), k2, 1, -1)
                                * Constant::Electron3j(e3.TwoJ(), e6.TwoJ(), k2, -e3.TwoM(), e6.TwoM());

                if(coeff36)
                {
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
                                // R1 = R_k1 (12, 4n)
                                double R1 = GetTwoElectronIntegral(k1, e1, e2, e4, en);

                                // R2 = R_k2 (n3, 56)
                                double R2 = GetTwoElectronIntegral(k2, en, e3, e5, e6);

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

    // phase = -1 * (-1)^(m1 - m6)
//    if(((e1.TwoM() - e6.TwoM())/2 + 1)%2 == 1)
//        total = - total;
    if(((e2.TwoM() + e3.TwoM() + e4.TwoM() + e5.TwoM())/2)%2 == 1)
        total = - total;

    return -total;
}

void Sigma3Calculator::UpdateStateIndexes()
{
    NumStates = core->NumStates() + excited->NumStates();
    state_index.clear();
    reverse_state_index.clear();

    ConstStateIterator it_i = core->GetConstStateIterator();
    unsigned int i;

    // Iterate through states, assign in order
    it_i.First(); i = 0;
    while(!it_i.AtEnd())
    {
        state_index.insert(std::pair<StateInfo, unsigned int>(StateInfo(it_i.GetState()), i));
        reverse_state_index.insert(std::pair<unsigned int, StateInfo>(i, StateInfo(it_i.GetState())));

        it_i.Next(); i++;
    }

    it_i = excited->GetConstStateIterator();
    it_i.First();
    while(!it_i.AtEnd())
    {
        // Check not already present from open shell core
        if(state_index.find(StateInfo(it_i.GetState())) == state_index.end())
        {   state_index.insert(std::pair<StateInfo, unsigned int>(StateInfo(it_i.GetState()), i));
            reverse_state_index.insert(std::pair<unsigned int, StateInfo>(i, StateInfo(it_i.GetState())));
        }

        it_i.Next(); i++;
    }
}

void Sigma3Calculator::Update()
{
    SMSIntegrals.clear();
    TwoElectronIntegrals.clear();

    StateIntegrator SI(excited->GetLattice());
    unsigned int i, j;

    const DiscreteState *si, *sj;

    // Get SMS integrals
    i = 0;
    while(i < NumStates)
    {
        si = core->GetState(reverse_state_index.find(i)->second);
        if(!si)
            si = excited->GetState(reverse_state_index.find(i)->second);

        j = i;
        while(j < NumStates)
        {
            sj = core->GetState(reverse_state_index.find(j)->second);
            if(!sj)
                sj = excited->GetState(reverse_state_index.find(j)->second);

            // i.pqn <= j.pqn, so calculate using derivative of i instead of j
            SMSIntegrals.insert(std::pair<unsigned int, double>(i* NumStates + j,
                -SI.IsotopeShiftIntegral(*sj, *si)));

            j++;
        }
        i++;
    }
}

double Sigma3Calculator::GetTwoElectronIntegral(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4)
{
    const double NuclearInverseMass = core->GetNuclearInverseMass();

    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;
    unsigned int i3 = state_index.find(s3)->second;
    unsigned int i4 = state_index.find(s4)->second;

    bool sms_sign = TwoElectronIntegralOrdering(i1, i2, i3, i4);

    unsigned int key = k  * NumStates*NumStates*NumStates*NumStates +
                       i1 * NumStates*NumStates*NumStates +
                       i2 * NumStates*NumStates +
                       i3 * NumStates +
                       i4;

    double radial = 0.;
    if(TwoElectronIntegrals.find(key) != TwoElectronIntegrals.end())
    {
        radial = TwoElectronIntegrals.find(key)->second;
    }
    else
    {   // Check triangle and parity conditions on k
        if(((k + s1.L() + s3.L())%2 == 1) ||
             (double(k) < fabs(s1.J() - s3.J())) ||
             (double(k) > s1.J() + s3.J()) ||
           ((k + s2.L() + s4.L())%2 == 1) ||
             (double(k) < fabs(s2.J() - s4.J())) ||
             (double(k) > s2.J() + s4.J()))
            return 0.;

        const State* s_1 = core->GetState(reverse_state_index.find(i1)->second);
        const State* s_2 = excited->GetState(reverse_state_index.find(i2)->second);
        const State* s_3 = excited->GetState(reverse_state_index.find(i3)->second);
        const State* s_4 = excited->GetState(reverse_state_index.find(i4)->second);

        unsigned int p;
        CoulombIntegrator CI(excited->GetLattice());
        const double* dR = excited->GetLattice()->dR();

        // Get density24
        std::vector<double> density(mmin(s_2->Size(), s_4->Size()));
        for(p=0; p<density.size(); p++)
        {
            density[p] = s_2->f[p] * s_4->f[p] + Constant::AlphaSquared * s_2->g[p] * s_4->g[p];
        }
        density.resize(core->GetHFPotential().size());

        // Get Pot24
        std::vector<double> Pot24(density.size());
        CI.FastCoulombIntegrate(density, Pot24, k);

        unsigned int limit = mmin(s_1->Size(), s_3->Size());
        limit = mmin(limit, Pot24.size());
        for(p=0; p<limit; p++)
        {
            radial += (s_1->f[p] * s_3->f[p] + Constant::AlphaSquared * s_1->g[p] * s_3->g[p])
                        * Pot24[p] * dR[p];
        }

        TwoElectronIntegrals.insert(std::pair<unsigned int, double>(key, radial));
    }

    if(NuclearInverseMass && (k == 1))
    {
        double SMS = NuclearInverseMass * SMSIntegrals.find(i1*NumStates + i3)->second * SMSIntegrals.find(i2*NumStates + i4)->second;
        if(!sms_sign)
            SMS = -SMS;
        radial = radial - SMS;
    }

    return radial;
}

double Sigma3Calculator::GetSMSIntegral(const StateInfo& s1, const StateInfo& s2) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;

    if(i1 <= i2)
        return SMSIntegrals.find(i1 * NumStates + i2)->second;
    else
        return -SMSIntegrals.find(i2 * NumStates + i1)->second;
}

bool Sigma3Calculator::TwoElectronIntegralOrdering(unsigned int& i1, unsigned int& i2, unsigned int& i3, unsigned int& i4) const
{
    bool sms_sign = true;

    // Ordering of indices:
    // (i1 <= i3) && (i2 <= i4) && (i1 <= i2) && (if i1 == i2, then (i3 <= i4))
    // therefore (i1 <= i2 <= i4) and (i1 <= i3)
    if(i3 < i1)
    {   swap(i3, i1);
        sms_sign = !sms_sign;
    }
    if(i4 < i2)
    {   swap(i4, i2);
        sms_sign = !sms_sign;
    }
    if(i2 < i1)
    {   swap(i2, i1);
        swap(i3, i4);
    }
    if((i1 == i2) && (i4 < i3))
        swap(i3, i4);

    return sms_sign;
}
