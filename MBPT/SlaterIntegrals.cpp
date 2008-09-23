#include "Include.h"
#include "SlaterIntegrals.h"

#include "HartreeFock/StateIntegrator.h"
#include "Universal/CoulombIntegrator.h"

void SlaterIntegrals::UpdateStateIndexes(const ExcitedStates& valence)
{
    NumStates = core.NumStates() + excited.NumStates();
    state_index.clear();
    reverse_state_index.clear();

    ConstStateIterator it_i = core.GetConstStateIterator();
    unsigned int i;

    // Iterate through states, assign in order
    it_i.First(); i = 0;
    while(!it_i.AtEnd())
    {
        core_states.insert(i);
        state_index.insert(std::pair<StateInfo, unsigned int>(StateInfo(it_i.GetState()), i));
        reverse_state_index.insert(std::pair<unsigned int, StateInfo>(i, StateInfo(it_i.GetState())));

        it_i.Next(); i++;
    }

    it_i = excited.GetConstStateIterator();
    it_i.First();
    while(!it_i.AtEnd())
    {
        StateInfo info(it_i.GetState());

        // Check not already present from open shell core
        if(state_index.find(info) == state_index.end())
        {   state_index.insert(std::pair<StateInfo, unsigned int>(info, i));
            reverse_state_index.insert(std::pair<unsigned int, StateInfo>(i, info));
        }

        excited_states.insert(i);
        if(valence.GetState(info))
            valence_states.insert(i);

        it_i.Next(); i++;
    }
}

unsigned int SlaterIntegrals::GetStorageSize(const ExcitedStates& valence)
{
    UpdateStateIndexes(valence);

    unsigned int size1 = 0;
    unsigned int size2 = 0;
    unsigned int sizesms = 0;

    std::set<unsigned int>::const_iterator it_n, it_alpha;

    *outstream << "Core states: " << core_states.size() << std::endl;
    *outstream << "Valence states: " << valence_states.size() << std::endl;
    *outstream << "Excited states: " << excited_states.size() << std::endl;
    *outstream << "NumStates: " << NumStates << std::endl;

    // One electron integrals: need all < n | H | alpha >
    // where n is core; alpha is excited (superset of valence)

    it_n = core_states.begin();
    while(it_n != core_states.end())
    {
        StateInfo n = reverse_state_index.find(*it_n)->second;
        it_alpha = excited_states.begin();
        while(it_alpha != excited_states.end())
        {
            StateInfo alpha = reverse_state_index.find(*it_alpha)->second;
            if(n.Kappa() == alpha.Kappa())
                size1++;

            it_alpha++;
        }

        it_n++;
    }

    *logstream << "Num one-electron integrals: " << size1 << std::endl;

    // SMS integrals: need all < n | p | alpha > and < a | p | alpha >
    // where n is core; a is valence; alpha is excited (superset of valence)

    // k == 1
    // triangle(k, j1, j2) -> j1 + j2 >= 1 (always satisfied)
    //                     -> |j1 - j2| <= 1
    // (l1 + l2 + 1)%2 == 0

    std::set<unsigned int>::const_iterator it_a;

    it_alpha = excited_states.begin();
    while(it_alpha != excited_states.end())
    {
        StateInfo alpha = reverse_state_index.find(*it_alpha)->second;
    
        it_n = core_states.begin();
        while(it_n != core_states.end())
        {
            StateInfo n = reverse_state_index.find(*it_n)->second;

            if((abs(alpha.TwoJ() - n.TwoJ()) <= 2) &&
               ((n.L() + alpha.L() + 1)%2 == 0))
                sizesms++;

            it_n++;
        }

        it_a = valence_states.begin();
        while((it_a != valence_states.end()) && (*it_a <= *it_alpha))
        {
            StateInfo a = reverse_state_index.find(*it_a)->second;

            if((abs(alpha.TwoJ() - a.TwoJ()) <= 2) &&
               ((a.L() + alpha.L() + 1)%2 == 0))
                sizesms++;

            it_a++;
        }

        it_alpha++;
    }

    *logstream << "Num sms integrals: " << sizesms << std::endl;

    // Two-electron Slater integrals for MBPT diagrams are of two types:
    //   R^k(n a, alpha beta), a <= beta
    //   R^k(n m, alpha beta), n <= m; if n == m then alpha <= beta
    // where n and m are core; a is valence; alpha and beta are excited (superset of valence)

    std::set<unsigned int>::const_iterator it_m, it_beta;
    unsigned int k, kmax;

    it_n = core_states.begin();
    while(it_n != core_states.end())
    {
        StateInfo n = reverse_state_index.find(*it_n)->second;
        
        it_alpha = excited_states.begin();
        while(it_alpha != excited_states.end())
        {
            StateInfo alpha = reverse_state_index.find(*it_alpha)->second;

            // (n, alpha) -> limits on k
            k = abs(int(n.L()) - int(alpha.L()));
            if(abs(int(n.TwoJ()) - int(alpha.TwoJ())) > 2 * k)
                k += 2;

            kmax = n.L() + alpha.L();
            if(n.TwoJ() + alpha.TwoJ() < 2 * kmax)
                kmax -= 2;

            while(k <= kmax)
            {
                it_beta = excited_states.begin();
                while(it_beta != excited_states.end())
                {
                    StateInfo beta = reverse_state_index.find(*it_beta)->second;

                    unsigned int Lsum = beta.L() + k;
                    unsigned int TwoJmin = abs(2*k - int(beta.TwoJ()));
                    unsigned int TwoJmax = 2*k + beta.TwoJ();

                    it_a = valence_states.begin();

                    // a <= beta, unless beta is not a valence state
                    unsigned int alimit;
                    if(valence_states.find(*it_beta) != valence_states.end())
                        alimit = *it_beta;
                    else
                        alimit = NumStates;

                    while((it_a != valence_states.end()) && (*it_a <= alimit))
                    {
                        StateInfo a = reverse_state_index.find(*it_a)->second;

                        if((a.TwoJ() >= TwoJmin) &&
                           (a.TwoJ() <= TwoJmax) &&
                           ((Lsum + a.L())%2 == 0))
                            size2++;

                        it_a++;
                    }

                    // n <= m; if n == m then alpha <= beta
                    it_m = it_n;
                    if(*it_alpha > *it_beta)
                        it_m++;

                    while(it_m != core_states.end())
                    {
                        StateInfo m = reverse_state_index.find(*it_m)->second;

                        if((m.TwoJ() >= TwoJmin) &&
                           (m.TwoJ() <= TwoJmax) &&
                           ((Lsum + m.L())%2 == 0))
                            size2++;

                        it_m++;
                    }
                    
                    it_beta++;
                }

                k += 2;
            }

            it_alpha++;
        }

        it_n++;
    }

    *logstream << "Num two-electron integrals: " << size2 << std::endl;

    return (size1 + size2 + sizesms);
}

void SlaterIntegrals::Clear()
{
    state_index.clear();
    reverse_state_index.clear();
    OneElectronIntegrals.clear();
    SMSIntegrals.clear();
    TwoElectronIntegrals.clear();
}

void SlaterIntegrals::Update(const ExcitedStates& valence)
{
    Clear();
    UpdateStateIndexes(valence);

    UpdateOneElectronIntegrals();
    UpdateTwoElectronIntegrals();
}

void SlaterIntegrals::UpdateOneElectronIntegrals()
{
    unsigned int key;
    double value;

    std::set<unsigned int>::const_iterator it_n, it_alpha;

    const DiscreteState* sn;
    const State* salpha;

    StateIntegrator SI(core.GetLattice());

    // One electron integrals: need all < n | H | alpha >
    // where n is core; alpha is excited (superset of valence)

    it_n = core_states.begin();
    while(it_n != core_states.end())
    {
        StateInfo n = reverse_state_index.find(*it_n)->second;
        sn = core.GetState(n);
        
        it_alpha = excited_states.begin();
        while(it_alpha != excited_states.end())
        {
            StateInfo alpha = reverse_state_index.find(*it_alpha)->second;

            if(n.Kappa() == alpha.Kappa())
            {
                salpha = excited.GetState(alpha);

                key = *it_n * NumStates + *it_alpha;
                value = SI.HamiltonianMatrixElement(*sn, *salpha, core);
                if(OneElectronIntegrals.find(key) != OneElectronIntegrals.end())
                    *outstream << "one electron overwrite" << std::endl;
                OneElectronIntegrals.insert(std::pair<unsigned int, double>(key, value));
            }

            it_alpha++;
        }

        it_n++;
    }

    // SMS integrals: need all < n | p | alpha > and < a | p | alpha >
    // where n is core; a is valence; alpha is excited (superset of valence)

    // k == 1
    // triangle(k, j1, j2) -> j1 + j2 >= 1 (always satisfied)
    //                     -> |j1 - j2| <= 1
    // (l1 + l2 + 1)%2 == 0

    std::set<unsigned int>::const_iterator it_a;
    const State* sa;

    it_alpha = excited_states.begin();
    while(it_alpha != excited_states.end())
    {
        StateInfo alpha = reverse_state_index.find(*it_alpha)->second;
        salpha = excited.GetState(alpha);

        it_n = core_states.begin();
        while(it_n != core_states.end())
        {
            StateInfo n = reverse_state_index.find(*it_n)->second;

            if((abs(alpha.TwoJ() - n.TwoJ()) <= 2) &&
               ((n.L() + alpha.L() + 1)%2 == 0))
            {
                sn = core.GetState(n);

                key = *it_n * NumStates + *it_alpha;
                value = -SI.IsotopeShiftIntegral(*salpha, *sn);
                if(SMSIntegrals.find(key) != SMSIntegrals.end())
                    *outstream << "sms overwrite" << std::endl;
                SMSIntegrals.insert(std::pair<unsigned int, double>(key, value));
            }

            it_n++;
        }

        it_a = valence_states.begin();
        while((it_a != valence_states.end()) && (*it_a <= *it_alpha))
        {
            StateInfo a = reverse_state_index.find(*it_a)->second;

            if((abs(alpha.TwoJ() - a.TwoJ()) <= 2) &&
               ((a.L() + alpha.L() + 1)%2 == 0))
            {
                sa = excited.GetState(a);

                if(*it_a < *it_alpha)
                {   key = *it_a * NumStates + *it_alpha;
                    value = -SI.IsotopeShiftIntegral(*salpha, *sa);
                }
                else
                {   key = *it_alpha * NumStates + *it_a;
                    value = -SI.IsotopeShiftIntegral(*sa, *salpha);
                }

                if(SMSIntegrals.find(key) != SMSIntegrals.end())
                    *outstream << "sms overwrite" << std::endl;
                SMSIntegrals.insert(std::pair<unsigned int, double>(key, value));
            }

            it_a++;
        }

        it_alpha++;
    }
}

void SlaterIntegrals::UpdateTwoElectronIntegrals()
{
    unsigned long long int key;
    double value;

    CoulombIntegrator CI(core.GetLattice());
    std::vector<double> density(core.GetHFPotential().size());
    std::vector<double> potential(core.GetHFPotential().size());
    const double* dR = core.GetLattice()->dR();
    unsigned int p;  // just a counter

    // Two-electron Slater integrals for MBPT diagrams are of two types:
    //   R^k(n a, alpha beta), a <= beta
    //   R^k(n m, alpha beta), n <= m; if n == m then alpha <= beta
    // where n and m are core; a is valence; alpha and beta are excited (superset of valence)

    std::set<unsigned int>::const_iterator it_n, it_m, it_a, it_alpha, it_beta;
    const DiscreteState *sn, *sm;
    const State *sa, *salpha, *sbeta;
    unsigned int k, kmax;

    it_n = core_states.begin();
    while(it_n != core_states.end())
    {
        StateInfo n = reverse_state_index.find(*it_n)->second;
        sn = core.GetState(n);

        it_alpha = excited_states.begin();
        while(it_alpha != excited_states.end())
        {
            StateInfo alpha = reverse_state_index.find(*it_alpha)->second;
            salpha = excited.GetState(alpha);

            // (n, alpha) -> limits on k
            k = abs(int(n.L()) - int(alpha.L()));
            if(abs(int(n.TwoJ()) - int(alpha.TwoJ())) > 2 * k)
                k += 2;

            kmax = n.TwoJ() + alpha.TwoJ();

            // Get density
            if(k <= kmax)
            {   for(p=0; p < mmin(sn->Size(), salpha->Size()); p++)
                {   density[p] = sn->f[p] * salpha->f[p] + Constant::AlphaSquared * sn->g[p] * salpha->g[p];
                }
            }

            while(k <= kmax)
            {
                // n, alpha, k -> potential
                CI.FastCoulombIntegrate(density, potential, k, mmin(sn->Size(), salpha->Size()));

                it_beta = excited_states.begin();
                while(it_beta != excited_states.end())
                {
                    StateInfo beta = reverse_state_index.find(*it_beta)->second;
                    sbeta = excited.GetState(beta);

                    unsigned int Lsum = beta.L() + k;
                    unsigned int TwoJmin = abs(int(2*k) - int(beta.TwoJ()));
                    unsigned int TwoJmax = 2*k + beta.TwoJ();

                    it_a = valence_states.begin();

                    // a <= beta, unless beta is not a valence state
                    unsigned int alimit;
                    if(valence_states.find(*it_beta) != valence_states.end())
                        alimit = *it_beta;
                    else
                        alimit = NumStates;

                    while((it_a != valence_states.end()) && (*it_a <= alimit))
                    {
                        StateInfo a = reverse_state_index.find(*it_a)->second;
                        sa = excited.GetState(a);

                        if((a.TwoJ() >= TwoJmin) &&
                           (a.TwoJ() <= TwoJmax) &&
                           ((Lsum + a.L())%2 == 0))
                        {
                            // Fix ordering
                            unsigned int i1, i2, i3, i4;
                            i1 = *it_n;
                            i2 = *it_a;
                            i3 = *it_alpha;
                            i4 = *it_beta;
                            TwoElectronIntegralOrdering(i1, i2, i3, i4);
                            key = k  * NumStates*NumStates*NumStates*NumStates +
                                  i1 * NumStates*NumStates*NumStates +
                                  i2 * NumStates*NumStates +
                                  i3 * NumStates +
                                  i4;

                            value = 0.;
                            unsigned int limit = mmin(sa->Size(), sbeta->Size());
                            limit = mmin(limit, potential.size());
                            for(p=0; p<limit; p++)
                            {
                                value += (sa->f[p] * sbeta->f[p] + Constant::AlphaSquared * sa->g[p] * sbeta->g[p])
                                         * potential[p] * dR[p];
                            }

                            TwoElectronIntegrals.insert(std::pair<unsigned long long int, double>(key, value));
                        }

                        it_a++;
                    }

                    // n <= m; if n == m then alpha <= beta
                    it_m = it_n;
                    if(*it_alpha > *it_beta)
                        it_m++;

                    while(it_m != core_states.end())
                    {
                        StateInfo m = reverse_state_index.find(*it_m)->second;
                        sm = core.GetState(m);

                        if((m.TwoJ() >= TwoJmin) &&
                           (m.TwoJ() <= TwoJmax) &&
                           ((Lsum + m.L())%2 == 0))
                        {
                            // Ordering should be correct
                            key = k         * NumStates*NumStates*NumStates*NumStates +
                                  *it_n     * NumStates*NumStates*NumStates +
                                  *it_m     * NumStates*NumStates +
                                  *it_alpha * NumStates +
                                  *it_beta;

                            value = 0.;
                            unsigned int limit = mmin(sm->Size(), sbeta->Size());
                            limit = mmin(limit, potential.size());
                            for(p=0; p<limit; p++)
                            {
                                value += (sm->f[p] * sbeta->f[p] + Constant::AlphaSquared * sm->g[p] * sbeta->g[p])
                                         * potential[p] * dR[p];
                            }

                            TwoElectronIntegrals.insert(std::pair<unsigned long long int, double>(key, value));
                        }

                        it_m++;
                    }
                    
                    it_beta++;
                }

                k += 2;
            }

            it_alpha++;
        }

        it_n++;
    }
}

double SlaterIntegrals::GetOneElectronIntegral(const StateInfo& s1, const StateInfo& s2) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;

    if(i1 <= i2)
        return OneElectronIntegrals.find(i1 * NumStates + i2)->second;
    else
        return OneElectronIntegrals.find(i2 * NumStates + i1)->second;
}

double SlaterIntegrals::GetSMSIntegral(const StateInfo& s1, const StateInfo& s2) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;

    if(i1 <= i2)
        return SMSIntegrals.find(i1 * NumStates + i2)->second;
    else
        return -SMSIntegrals.find(i2 * NumStates + i1)->second;
}

bool SlaterIntegrals::TwoElectronIntegralOrdering(unsigned int& i1, unsigned int& i2, unsigned int& i3, unsigned int& i4) const
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

double SlaterIntegrals::GetTwoElectronIntegral(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;
    unsigned int i3 = state_index.find(s3)->second;
    unsigned int i4 = state_index.find(s4)->second;

    bool sms_sign = TwoElectronIntegralOrdering(i1, i2, i3, i4);

    unsigned long long int key = k  * NumStates*NumStates*NumStates*NumStates +
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
