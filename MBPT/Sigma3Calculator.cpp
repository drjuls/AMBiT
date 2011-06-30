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
        OrbitalInfo b = reverse_state_index.find(*it_b)->second;

        it_n = core_states.begin();
        while(it_n != core_states.end())
        {
            OrbitalInfo n = reverse_state_index.find(*it_n)->second;

            if((abs(b.TwoJ() - n.TwoJ()) <= 2) &&
               ((n.L() + b.L() + 1)%2 == 0))
                sizesms++;

            it_n++;
        }

        it_a = valence_states.begin();
        while((it_a != valence_states.end()) && (*it_a <= *it_b))
        {
            OrbitalInfo a = reverse_state_index.find(*it_a)->second;

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
        OrbitalInfo n = reverse_state_index.find(*it_n)->second;

        it_a = valence_states.begin();
        while(it_a != valence_states.end())
        {
            OrbitalInfo a = reverse_state_index.find(*it_a)->second;

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
                    OrbitalInfo c = reverse_state_index.find(*it_c)->second;

                    unsigned int Lsum = c.L() + k;
                    unsigned int TwoJmin = abs(int(2*k) - int(c.TwoJ()));
                    unsigned int TwoJmax = 2*k + c.TwoJ();

                    it_b = valence_states.begin();
                    while((it_b != valence_states.end()) && (*it_b <= *it_c))
                    {
                        OrbitalInfo b = reverse_state_index.find(*it_b)->second;

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
    const Orbital* sn, *sa, *sb;

    it_b = valence_states.begin();
    while(it_b != valence_states.end())
    {
        OrbitalInfo b = reverse_state_index.find(*it_b)->second;
        sb = excited.GetState(b);

        it_n = core_states.begin();
        while(it_n != core_states.end())
        {
            OrbitalInfo n = reverse_state_index.find(*it_n)->second;

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
            OrbitalInfo a = reverse_state_index.find(*it_a)->second;

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
    const Orbital *sn;
    const SingleParticleWavefunction *sa, *sb, *sc;
    unsigned int k, kmax;

    it_n = core_states.begin();
    while(it_n != core_states.end())
    {
        OrbitalInfo n = reverse_state_index.find(*it_n)->second;
        sn = core.GetState(n);

        it_a = valence_states.begin();
        while(it_a != valence_states.end())
        {
            OrbitalInfo a = reverse_state_index.find(*it_a)->second;
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
                    OrbitalInfo c = reverse_state_index.find(*it_c)->second;
                    sc = excited.GetState(c);

                    unsigned int Lsum = c.L() + k;
                    unsigned int TwoJmin = abs(int(2*k) - int(c.TwoJ()));
                    unsigned int TwoJmax = 2*k + c.TwoJ();

                    it_b = valence_states.begin();
                    while((it_b != valence_states.end()) && (*it_b <= *it_c))
                    {
                        OrbitalInfo b = reverse_state_index.find(*it_b)->second;
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

                            value = value * Constant::Electron3j(n.TwoJ(), a.TwoJ(), k)
                                          * Constant::Electron3j(b.TwoJ(), c.TwoJ(), k)
                                    * sqrt(double(n.MaxNumElectrons() * a.MaxNumElectrons() *
                                                  b.MaxNumElectrons() * c.MaxNumElectrons()));
                                    

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

double Sigma3Integrals::GetTwoElectronIntegral(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;
    unsigned int i3 = state_index.find(s3)->second;
    unsigned int i4 = state_index.find(s4)->second;

    bool sms_sign = true;

    // Ordering of indices, assuming i1 is core, others are valence, and therefore i1 is smallest:
    // (i2 <= i4)
    // therefore (i1 <= i2 <= i4) and (i1 <= i3)
    if(i4 < i2)
    {   swap(i4, i2);
        sms_sign = false;
    }

    LongKey key = k  * NumStates*NumStates*NumStates*NumStates +
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
    {   *errstream << "SlaterIntegrals::GetTwoElectronIntegral() failed to find integral."
                   << "\n  key = " << key << "  num_states = " << NumStates 
                   << "\n  R^" << k << " ( " << s1.Name() << " " << s2.Name()
                   << ", " << s3.Name() << " " << s4.Name() << ") :" << std::endl;
    }

    if(include_valence_sms && (k == 1))
    {   double SMS = GetNuclearInverseMass();
        if(SMS)
        {   SMS = SMS * SMSIntegrals.find(i1*NumStates + i3)->second * SMSIntegrals.find(i2*NumStates + i4)->second;
            if(!sms_sign)
                SMS = -SMS;
            radial = radial - SMS;
        }
    }

    return radial;
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

    // Get maximum value of twoJ in core
    core_maxTwoJ = 0;
    ConstStateIterator it = core->GetConstStateIterator();
    while(!it.AtEnd())
    {
        if(core_maxTwoJ < it.GetState()->TwoJ())
            core_maxTwoJ = it.GetState()->TwoJ();
        it.Next();
    }
}

double Sigma3Calculator::GetSecondOrderSigma3(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
           const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6)
{
    // core state limits
    int two_mn = e1.TwoM() + e2.TwoM() - e4.TwoM();
    if(core_maxTwoJ < abs(two_mn))
        return 0.;

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

    if((k1min > k1max) || (k2min > k2max))
        return 0.;

    int ln_parity = (k1min + e2.L())%2;
    const double ValenceEnergy = ValenceEnergies.find(e2.Kappa())->second;

    double total = 0.;
    unsigned int k1, k2;

    // Summation over k1
    for(k1 = k1min; k1 <= k1max; k1+=2)
    {
        double coeff14 = Constant::Electron3j(e1.TwoJ(), e4.TwoJ(), k1, -e1.TwoM(), e4.TwoM());

        if(coeff14)
        {
            // Summation over k2
            for(k2 = k2min; k2 <= k2max; k2+=2)
            {
                double coeff36 = Constant::Electron3j(e3.TwoJ(), e6.TwoJ(), k2, -e3.TwoM(), e6.TwoM());

                if(coeff36)
                {
                    // Summation over core state 'n'
                    ConstStateIterator itn = core->GetConstStateIterator();
                    while(!itn.AtEnd())
                    {
                        const OrbitalInfo sn = itn.GetOrbitalInfo();
                        unsigned int two_Jn = sn.TwoJ();

                        if((two_Jn >= abs(two_mn)) && (sn.L()%2 == ln_parity))
                        {
                            const double En = itn.GetState()->Energy();

                            double coeff = Constant::Electron3j(e2.TwoJ(), two_Jn, k1, -e2.TwoM(), two_mn)
                                * Constant::Electron3j(two_Jn, e5.TwoJ(), k2, -two_mn, e5.TwoM())
                                * coeff14 * coeff36 /(En - ValenceEnergy + delta);

                            if(coeff)
                            {
                                // R1 = R_k1 (12, 4n) * (j1   j4  k1) * (j2   jn  k1) * sqrt[j1, j2, j4, jn]
                                //                      (1/2 -1/2  0) * (1/2 -1/2  0)
                                //    = R_k1 (n4, 21) * (jn   j2  k1) * (j4   j1  k1) * sqrt[jn, j4, j2, j1]
                                //                      (1/2 -1/2  0) * (1/2 -1/2  0)
                                double R1 = integrals->GetTwoElectronIntegral(k1, sn, e4, e2, e1);

                                // R2 = R_k2 (n3, 56) * (jn   j5  k2) * (j3   j6  k2) * sqrt[jn, j3, j5, j6]
                                //                      (1/2 -1/2  0) * (1/2 -1/2  0)
                                double R2 = integrals->GetTwoElectronIntegral(k2, sn, e3, e5, e6);

                                total += coeff * R1 * R2;
                            }
                        }
                        itn.Next();
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
