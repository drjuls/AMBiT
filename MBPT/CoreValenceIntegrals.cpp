#include "Include.h"
#include "CoreValenceIntegrals.h"
#include "HartreeFock/StateIntegrator.h"
#include "Universal/CoulombIntegrator.h"

unsigned int CoreValenceIntegrals::GetStorageSize(const ExcitedStates& valence)
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

            kmax = (n.TwoJ() + alpha.TwoJ())/2;

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

void CoreValenceIntegrals::UpdateOneElectronIntegrals()
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

        // a <= alpha, unless alpha is not a valence state
        unsigned int alimit;
        if(valence_states.find(*it_alpha) != valence_states.end())
            alimit = *it_alpha;
        else
            alimit = NumStates;

        it_a = valence_states.begin();
        while((it_a != valence_states.end()) && (*it_a <= alimit))
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

void CoreValenceIntegrals::UpdateTwoElectronIntegrals()
{
    LongKey key;
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

            kmax = (n.TwoJ() + alpha.TwoJ())/2;

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

                            TwoElectronIntegrals.insert(std::pair<LongKey, double>(key, value));
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

                            TwoElectronIntegrals.insert(std::pair<LongKey, double>(key, value));
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
