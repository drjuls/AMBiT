#include "Sigma3Calculator.h"
#include "Include.h"
#include "Universal/Constant.h"
#include "Universal/CoulombIntegrator.h"
#include "HartreeFock/StateIntegrator.h"

unsigned int Sigma3Integrals::GetStorageSize(const ExcitedStates& valence)
{
    UpdateStateIndexes(valence);

    unsigned int sizesms = 0;
    unsigned int size2 = 0;

    // SMS integrals: need all < n | p | b > and < a | p | b >
    // where n is core; a and b are valence

    // k == 1
    // triangle(k, j1, j2) -> j1 + j2 >= 1 (always satisfied)
    //                     -> |j1 - j2| <= 1
    // (l1 + l2 + 1)%2 == 0

    std::set<unsigned int>::const_iterator it_n, it_a, it_b;

    it_b = valence_states.begin();
    while(it_b != valence_states.end())
    {
        StateInfo b = reverse_state_index.find(*it_b)->second;

        it_n = core_states.begin();
        while(it_n != core_states.end())
        {
            StateInfo n = reverse_state_index.find(*it_n)->second;

            if((abs(b.TwoJ() - n.TwoJ()) <= 2) &&
               ((n.L() + b.L() + 1)%2 == 0))
                sizesms++;

            it_n++;
        }

        it_a = valence_states.begin();
        while((it_a != valence_states.end()) && (*it_a <= *it_b))
        {
            StateInfo a = reverse_state_index.find(*it_a)->second;

            if((abs(b.TwoJ() - a.TwoJ()) <= 2) &&
               ((a.L() + b.L() + 1)%2 == 0))
                sizesms++;

            it_a++;
        }

        it_b++;
    }

    *logstream << "Num sms integrals: " << sizesms << std::endl;

     // Two-electron Slater integrals for Sigma3 diagrams are of the form:
    //   R^k(n b, a c), b <= c
    // where n is core; a, b and c are valence

    std::set<unsigned int>::const_iterator it_c;
    unsigned int k, kmax;

    it_n = core_states.begin();
    while(it_n != core_states.end())
    {
        StateInfo n = reverse_state_index.find(*it_n)->second;

        it_a = valence_states.begin();
        while(it_a != valence_states.end())
        {
            StateInfo a = reverse_state_index.find(*it_a)->second;

            // (n, a) -> limits on k
            k = abs(int(n.L()) - int(a.L()));
            if(abs(int(n.TwoJ()) - int(a.TwoJ())) > 2 * k)
                k += 2;

            kmax = n.TwoJ() + a.TwoJ();

            while(k <= kmax)
            {
                it_c = valence_states.begin();
                while(it_c != valence_states.end())
                {
                    StateInfo c = reverse_state_index.find(*it_c)->second;

                    unsigned int Lsum = c.L() + k;
                    unsigned int TwoJmin = abs(int(2*k) - int(c.TwoJ()));
                    unsigned int TwoJmax = 2*k + c.TwoJ();

                    it_b = valence_states.begin();
                    while((it_b != valence_states.end()) && (*it_b <= *it_c))
                    {
                        StateInfo b = reverse_state_index.find(*it_b)->second;

                        if((b.TwoJ() >= TwoJmin) &&
                           (b.TwoJ() <= TwoJmax) &&
                           ((Lsum + b.L())%2 == 0))
                            size2++;

                        it_b++;
                    }
                    it_c++;
                }
                k += 2;
            }
            it_a++;
        }
        it_n++;
    }

   *logstream << "Num two-electron integrals: " << size2 << std::endl;

    return (size2 + sizesms);
}

void Sigma3Integrals::UpdateOneElectronIntegrals()
{
    unsigned int key;
    double value;

    StateIntegrator SI(core.GetLattice());

    // SMS integrals: need all < n | p | b > and < a | p | b >
    // where n is core; a and b are valence

    // k == 1
    // triangle(k, j1, j2) -> j1 + j2 >= 1 (always satisfied)
    //                     -> |j1 - j2| <= 1
    // (l1 + l2 + 1)%2 == 0

    std::set<unsigned int>::const_iterator it_n, it_a, it_b;
    const DiscreteState* sn, *sa, *sb;

    it_b = valence_states.begin();
    while(it_b != valence_states.end())
    {
        StateInfo b = reverse_state_index.find(*it_b)->second;
        sb = excited.GetState(b);

        it_n = core_states.begin();
        while(it_n != core_states.end())
        {
            StateInfo n = reverse_state_index.find(*it_n)->second;

            if((abs(b.TwoJ() - n.TwoJ()) <= 2) &&
               ((n.L() + b.L() + 1)%2 == 0))
            {
                sn = core.GetState(n);

                key = *it_n * NumStates + *it_b;
                value = -SI.IsotopeShiftIntegral(*sb, *sn);
                if(SMSIntegrals.find(key) != SMSIntegrals.end())
                    *outstream << "sms overwrite" << std::endl;
                SMSIntegrals.insert(std::pair<unsigned int, double>(key, value));
            }

            it_n++;
        }

        it_a = valence_states.begin();
        while((it_a != valence_states.end()) && (*it_a <= *it_b))
        {
            StateInfo a = reverse_state_index.find(*it_a)->second;

            if((abs(b.TwoJ() - a.TwoJ()) <= 2) &&
               ((a.L() + b.L() + 1)%2 == 0))
            {
                sa = excited.GetState(a);

                key = *it_a * NumStates + *it_b;
                value = -SI.IsotopeShiftIntegral(*sb, *sa);

                if(SMSIntegrals.find(key) != SMSIntegrals.end())
                    *outstream << "sms overwrite" << std::endl;
                SMSIntegrals.insert(std::pair<unsigned int, double>(key, value));
            }

            it_a++;
        }

        it_b++;
    }
}

void Sigma3Integrals::UpdateTwoElectronIntegrals()
{
    LongKey key;
    double value;

    CoulombIntegrator CI(core.GetLattice());
    std::vector<double> density(core.GetHFPotential().size());
    std::vector<double> potential(core.GetHFPotential().size());
    const double* dR = core.GetLattice()->dR();
    unsigned int p;  // just a counter

    // Two-electron Slater integrals for Sigma3 diagrams are of the form:
    //   R^k(n b, a c), b <= c
    // where n is core; a, b and c are valence

    std::set<unsigned int>::const_iterator it_n, it_a, it_b, it_c;
    const DiscreteState *sn;
    const State *sa, *sb, *sc;
    unsigned int k, kmax;

    it_n = core_states.begin();
    while(it_n != core_states.end())
    {
        StateInfo n = reverse_state_index.find(*it_n)->second;
        sn = core.GetState(n);

        it_a = valence_states.begin();
        while(it_a != valence_states.end())
        {
            StateInfo a = reverse_state_index.find(*it_a)->second;
            sa = excited.GetState(a);

            // (n, a) -> limits on k
            k = abs(int(n.L()) - int(a.L()));
            if(abs(int(n.TwoJ()) - int(a.TwoJ())) > 2 * k)
                k += 2;

            kmax = n.TwoJ() + a.TwoJ();

            // Get density
            if(k <= kmax)
            {   for(p=0; p < mmin(sn->Size(), sa->Size()); p++)
                {   density[p] = sn->f[p] * sa->f[p] + Constant::AlphaSquared * sn->g[p] * sa->g[p];
                }
            }

            while(k <= kmax)
            {
                // n, a, k -> potential
                CI.FastCoulombIntegrate(density, potential, k, mmin(sn->Size(), sa->Size()));

                it_c = valence_states.begin();
                while(it_c != valence_states.end())
                {
                    StateInfo c = reverse_state_index.find(*it_c)->second;
                    sc = excited.GetState(c);

                    unsigned int Lsum = c.L() + k;
                    unsigned int TwoJmin = abs(int(2*k) - int(c.TwoJ()));
                    unsigned int TwoJmax = 2*k + c.TwoJ();

                    it_b = valence_states.begin();
                    while((it_b != valence_states.end()) && (*it_b <= *it_c))
                    {
                        StateInfo b = reverse_state_index.find(*it_b)->second;
                        sb = excited.GetState(b);

                        if((b.TwoJ() >= TwoJmin) &&
                           (b.TwoJ() <= TwoJmax) &&
                           ((Lsum + b.L())%2 == 0))
                        {
                            // Ordering should be correct
                            key = k     * NumStates*NumStates*NumStates*NumStates +
                                  *it_n * NumStates*NumStates*NumStates +
                                  *it_b * NumStates*NumStates +
                                  *it_a * NumStates +
                                  *it_c;

                            value = 0.;
                            unsigned int limit = mmin(sb->Size(), sc->Size());
                            limit = mmin(limit, potential.size());
                            for(p=0; p<limit; p++)
                            {
                                value += (sb->f[p] * sc->f[p] + Constant::AlphaSquared * sb->g[p] * sc->g[p])
                                         * potential[p] * dR[p];
                            }

                            TwoElectronIntegrals.insert(std::pair<LongKey, double>(key, value));
                        }
                        it_b++;
                    }
                    it_c++;
                }
                k += 2;
            }
            it_a++;
        }
        it_n++;
    }
}

Sigma3Calculator::Sigma3Calculator(Lattice* lattice, const Core* atom_core, const ExcitedStates* excited_states):
    MBPTCalculator(lattice, atom_core, excited_states), integrals(NULL)
{}

Sigma3Calculator::~Sigma3Calculator(void)
{
    if(integrals)
        delete integrals;
}

unsigned int Sigma3Calculator::GetStorageSize(const ExcitedStates* valence_states)
{
    if(!integrals)
        integrals = new Sigma3Integrals(excited);
 
    return integrals->GetStorageSize(*valence_states);
}

void Sigma3Calculator::UpdateIntegrals(const ExcitedStates* valence_states)
{
    SetValenceEnergies();
    if(integrals)
        delete integrals;
    
    integrals = new Sigma3Integrals(excited);
    integrals->Update(*valence_states);
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
    unsigned int k1min = kmin(e1, e4);
    while(k1min < abs(q1))
        k1min += 2;
    unsigned int k1max = kmax(e1, e4);

    // k2 limits
    int q2 = (e6.TwoM() - e3.TwoM())/2;
    unsigned int k2min = kmin(e3, e6);
    while(k2min < abs(q2))
        k2min += 2;
    unsigned int k2max = kmax(e3, e6);

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
                                double R1 = integrals->GetTwoElectronIntegral(k1, e1, e2, e4, en);

                                // R2 = R_k2 (n3, 56)
                                double R2 = integrals->GetTwoElectronIntegral(k2, en, e3, e5, e6);

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
