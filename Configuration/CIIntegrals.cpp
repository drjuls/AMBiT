#include "Include.h"
#include "CIIntegrals.h"
#include "HartreeFock/StateIntegrator.h"
#include "Universal/CoulombIntegrator.h"
#include "MBPT/MBPTCalculator.h"

CIIntegrals::~CIIntegrals()
{
    std::map<int, SigmaPotential*>::iterator it = Sigma1.begin();
    while(it != Sigma1.end())
    {   delete it->second;
        it++;
    }
    Sigma1.clear();
}

unsigned int CIIntegrals::GetStorageSize() const
{
    unsigned int num_states = states.NumStates();
    unsigned int size = 0;

    unsigned int i1, i2, i3, i4;
    unsigned int k, kmax;

    ConstStateIterator it_1 = states.GetConstStateIterator();
    ConstStateIterator it_2 = states.GetConstStateIterator();
    ConstStateIterator it_3 = states.GetConstStateIterator();
    ConstStateIterator it_4 = states.GetConstStateIterator();

    it_2.First(); i2 = 0;
    while(!it_2.AtEnd())
    {
        const DiscreteState* s_2 = it_2.GetState();
        it_4 = it_2; i4 = i2;
        while(!it_4.AtEnd())
        {
            const DiscreteState* s_4 = it_4.GetState();

            // Limits on k
            k = abs(int(s_2->L()) - int(s_4->L()));
            if(fabs(s_2->J() - s_4->J()) > double(k))
                k += 2;

            kmax = s_2->L() + s_4->L();
            if(s_2->J() + s_4->J() < double(kmax))
                kmax -= 2;

            while(k <= kmax)
            {
                // s1 is the smallest
                it_1.First(); i1 = 0;
                while((i1 <= i2) && (it_1.GetState()->RequiredPQN() <= max_pqn_1))
                {
                    const DiscreteState* s_1 = it_1.GetState();

                    it_3 = it_1; i3 = i1;
                    unsigned int i3_limit;
                    if(i1 == i2)
                        i3_limit = i4;
                    else
                        i3_limit = num_states;
                    while((i3 <= i3_limit) && !it_3.AtEnd())
                    {
                        const DiscreteState* s_3 = it_3.GetState();

                        // Check max_pqn conditions and k conditions
                        if(((s_2->RequiredPQN() <= max_pqn_2) || (s_3->RequiredPQN() <= max_pqn_2)) &&
                           (s_2->RequiredPQN() <= max_pqn_3) && (s_3->RequiredPQN() <= max_pqn_3) &&
                           (int(k) >= abs(int(s_1->L()) - int(s_3->L()))) &&
                           (double(k) >= fabs(s_1->J() - s_3->J())) &&
                           (k <= s_1->L() + s_3->L()) &&
                           (double(k) <= s_1->J() + s_3->J()))
                        {
                                size++;
                        }

                        it_3.Next(); i3++;
                    }
                    it_1.Next(); i1++;
                }
                k+=2;
            }
            it_4.Next(); i4++;
        }
        it_2.Next(); i2++;
    }

    *logstream << "Num Integrals: " << size << std::endl;
    return size;
}

void CIIntegrals::Update(const std::string& id)
{
    StateIntegrator SI(*states.GetLattice());
    CoulombIntegrator CI(*states.GetLattice());

    NumStates = states.NumStates();
    state_index.clear();
    OneElectronIntegrals.clear();
    SMSIntegrals.clear();
    TwoElectronIntegrals.clear();

    unsigned int i, j;

    ConstStateIterator it_i = states.GetConstStateIterator();
    ConstStateIterator it_j = states.GetConstStateIterator();

    // Calculate sigma1 potentials
    if(include_sigma1)
    {
        std::map<int, SigmaPotential*>::iterator it_sigma = Sigma1.begin();
        while(it_sigma != Sigma1.end())
        {   delete it_sigma->second;
            it_sigma++;
        }
        Sigma1.clear();

        unsigned int max_l = 0;
        it_i.First();
        while(!it_i.AtEnd())
        {   max_l = mmax(it_i.GetState()->L(), max_l);
            it_i.Next();
        }

        for(int kappa = - (int)max_l - 1; kappa <= (int)max_l; kappa++)
            if(kappa != 0)
            {
                std::string sigma_id = id + ".";
                if(kappa > 0)
                    sigma_id = sigma_id + Constant::SpectroscopicNotation[kappa];
                else
                    sigma_id = sigma_id + Constant::SpectroscopicNotation[-kappa-1];

                if(kappa < -1)
                    sigma_id = sigma_id + '+';

                sigma_id = sigma_id + ".matrix";

                SigmaPotential* pot = new SigmaPotential(states.GetLattice(), sigma_id);

                if(!pot->Size() && PT)
                {   double valence_energy = 0.;
                    unsigned int pqn = 10;

                    // Get leading state (for energy denominator)
                    it_i.First();
                    while(!it_i.AtEnd())
                    {   const DiscreteState* ds = it_i.GetState();
                        if((ds->Kappa() == kappa) && (ds->RequiredPQN() < pqn))
                        {   pqn = ds->RequiredPQN();
                            valence_energy = ds->Energy();
                        }
                        it_i.Next();
                    }

                    PT->UseBrillouinWignerPT(valence_energy);
                    PT->GetSecondOrderSigma(kappa, pot);

                    if(!id.empty())
                        pot->Store();
                }

                Sigma1.insert(std::pair<int, SigmaPotential*>(kappa, pot));
            }
    }

    // Get single particle integrals
    it_i.First(); i = 0;
    while(!it_i.AtEnd())
    {
        // Use it_i to build index of states
        state_index.insert(std::pair<StateInfo, unsigned int>(StateInfo(it_i.GetState()), i));
        reverse_state_index.insert(std::pair<unsigned int, StateInfo>(i, StateInfo(it_i.GetState())));

        it_j = it_i; j = i;
        while(!it_j.AtEnd())
        {
            const DiscreteState* si = it_i.GetState();
            const DiscreteState* sj = it_j.GetState();

            if(si->Kappa() == sj->Kappa())
            {   double integral = SI.HamiltonianMatrixElement(*it_i.GetState(), *it_j.GetState(), *states.GetCore());
                if(include_sigma1)
                {   SigmaPotential* pot = Sigma1[si->Kappa()];
                    integral += pot->GetMatrixElement(si->f, sj->f);
                }
                OneElectronIntegrals.insert(std::pair<unsigned int, double>(i* NumStates + j, integral));
            }

            // i.pqn <= j.pqn, so calculate using derivative of i instead of j
            SMSIntegrals.insert(std::pair<unsigned int, double>(i* NumStates + j,
                -SI.IsotopeShiftIntegral(*it_j.GetState(), *it_i.GetState())));

            // Overlap integrals
            const State& p1 = *it_i.GetState();
            const State& p2 = *it_j.GetState();
            if(p1.L() == p2.L())
            {
                double overlap = 0.;
                const double* dR = states.GetLattice()->dR();
                for(unsigned int x=0; x<mmin(p1.Size(), p2.Size()); x++)
                    overlap += (p1.f[x] * p2.f[x] + Constant::AlphaSquared * p1.g[x] * p2.g[x]) * dR[x];

                OverlapIntegrals.insert(std::pair<unsigned int, double>(i * NumStates + j, overlap));
            }

            it_j.Next(); j++;
        }
        it_i.Next(); i++;
    }

    // Two electron integrals
    const double* dR = states.GetLattice()->dR();
    const double* R = states.GetLattice()->R();
    const double core_pol = states.GetCore()->GetPolarisability();
    const double core_rad = states.GetCore()->GetClosedShellRadius();

    unsigned int i1, i2, i3, i4;
    unsigned int k, kmax;
    unsigned int p;  // just a counter

    ConstStateIterator it_1 = states.GetConstStateIterator();
    ConstStateIterator it_2 = states.GetConstStateIterator();
    ConstStateIterator it_3 = states.GetConstStateIterator();
    ConstStateIterator it_4 = states.GetConstStateIterator();

    // Get 2 -> 4
    it_2.First(); i2 = 0;
    while(!it_2.AtEnd())
    {
        const DiscreteState* s_2 = it_2.GetState();
        it_4 = it_2; i4 = i2;
        while(!it_4.AtEnd())
        {
            const DiscreteState* s_4 = it_4.GetState();

            // Limits on k
            k = abs(int(s_2->L()) - int(s_4->L()));
            if(fabs(s_2->J() - s_4->J()) > double(k))
                k += 2;

            kmax = s_2->L() + s_4->L();
            if(s_2->J() + s_4->J() < double(kmax))
                kmax -= 2;

            // Get density24
            std::vector<double> density(mmin(s_2->Size(), s_4->Size()));
            if(k <= kmax)
            {   
                for(p=0; p<density.size(); p++)
                {
                    density[p] = s_2->f[p] * s_4->f[p] + Constant::AlphaSquared * s_2->g[p] * s_4->g[p];
                }
                density.resize(states.GetCore()->GetHFPotential().size());
            }

            while(k <= kmax)
            {
                // Get Pot24
                std::vector<double> Pot24(density.size());
                CI.FastCoulombIntegrate(density, Pot24, k);

                // s1 is the smallest
                it_1.First(); i1 = 0;
                while((i1 <= i2) && (it_1.GetState()->RequiredPQN() <= max_pqn_1))
                {
                    const DiscreteState* s_1 = it_1.GetState();

                    it_3 = it_1; i3 = i1;
                    unsigned int i3_limit;
                    if(i1 == i2)
                        i3_limit = i4;
                    else
                        i3_limit = NumStates;
                    while((i3 <= i3_limit) && !it_3.AtEnd())
                    {
                        const DiscreteState* s_3 = it_3.GetState();

                        // Check max_pqn conditions and k conditions
                        if(((s_2->RequiredPQN() <= max_pqn_2) || (s_3->RequiredPQN() <= max_pqn_2)) &&
                           (s_2->RequiredPQN() <= max_pqn_3) && (s_3->RequiredPQN() <= max_pqn_3) &&
                           (int(k) >= abs(int(s_1->L()) - int(s_3->L()))) &&
                           (double(k) >= fabs(s_1->J() - s_3->J())) &&
                           (k <= s_1->L() + s_3->L()) &&
                           (double(k) <= s_1->J() + s_3->J()))
                        {
                            unsigned int key = k  * NumStates*NumStates*NumStates*NumStates +
                                               i1 * NumStates*NumStates*NumStates +
                                               i2 * NumStates*NumStates + 
                                               i3 * NumStates +
                                               i4;

                            double radial = 0.;
                            unsigned int limit = mmin(s_1->Size(), s_3->Size());
                            limit = mmin(limit, Pot24.size());
                            for(p=0; p<limit; p++)
                            {
                                radial += (s_1->f[p] * s_3->f[p] + Constant::AlphaSquared * s_1->g[p] * s_3->g[p])
                                          * Pot24[p] * dR[p];
                            }

                            if(core_pol && k == 1)
                            {   
                                double R1 = 0.;
                                double R2 = 0.;
                                for(p=0; p<limit; p++)
                                {
                                    double r2 = R[p]*R[p] + core_rad*core_rad;
                                    R1 += (s_1->f[p] * s_3->f[p] + Constant::AlphaSquared * s_1->g[p] * s_3->g[p])/r2 * dR[p];
                                    R2 += density[p]/r2 * dR[p];
                                }

                                radial -= core_pol * R1 * R2;
                            }

                            TwoElectronIntegrals.insert(std::pair<unsigned int, double>(key, radial));
                        }
                        it_3.Next(); i3++;
                    }
                    it_1.Next(); i1++;
                }
                k+=2;
            }
            it_4.Next(); i4++;
        }
        it_2.Next(); i2++;
    }
}

inline void swap(unsigned int& i1, unsigned int& i2)
{   unsigned int temp = i1;
    i1 = i2;
    i2 = temp;
}

double CIIntegrals::GetOneElectronIntegral(const StateInfo& s1, const StateInfo& s2) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;

    if(i1 <= i2)
        return OneElectronIntegrals.find(i1 * NumStates + i2)->second;
    else
        return OneElectronIntegrals.find(i2 * NumStates + i1)->second;
}

double CIIntegrals::GetSMSIntegral(const StateInfo& s1, const StateInfo& s2) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;

    if(i1 <= i2)
        return SMSIntegrals.find(i1 * NumStates + i2)->second;
    else
        return -SMSIntegrals.find(i2 * NumStates + i1)->second;
}

double CIIntegrals::GetOverlapIntegral(const StateInfo& s1, const StateInfo& s2) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;

    if(i1 <= i2)
        return OverlapIntegrals.find(i1 * NumStates + i2)->second;
    else
        return OverlapIntegrals.find(i2 * NumStates + i1)->second;
}

double CIIntegrals::GetTwoElectronIntegral(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;
    unsigned int i3 = state_index.find(s3)->second;
    unsigned int i4 = state_index.find(s4)->second;

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
    {   // Check triangle conditions on k
        if((int(k) < abs(int(s1.L()) - int(s3.L()))) ||
             (double(k) < fabs(s1.J() - s3.J())) ||
           (k > s1.L() + s3.L()) ||
             (double(k) > s1.J() + s3.J()) ||
           (int(k) < abs(int(s2.L()) - int(s4.L()))) ||
             (double(k) < fabs(s2.J() - s4.J())) ||
           (k > s2.L() + s4.L()) ||
             (double(k) > s2.J() + s4.J()))
            return 0.;

        const State* s_1 = states.GetState(reverse_state_index.find(i1)->second);
        const State* s_2 = states.GetState(reverse_state_index.find(i2)->second);
        const State* s_3 = states.GetState(reverse_state_index.find(i3)->second);
        const State* s_4 = states.GetState(reverse_state_index.find(i4)->second);

        unsigned int p;
        CoulombIntegrator CI(*states.GetLattice());
        const double* R = states.GetLattice()->R();
        const double* dR = states.GetLattice()->dR();
        const double core_pol = states.GetCore()->GetPolarisability();
        const double core_rad = states.GetCore()->GetClosedShellRadius();

        // Get density24
        std::vector<double> density(mmin(s_2->Size(), s_4->Size()));
        for(p=0; p<density.size(); p++)
        {
            density[p] = s_2->f[p] * s_4->f[p] + Constant::AlphaSquared * s_2->g[p] * s_4->g[p];
        }
        density.resize(states.GetCore()->GetHFPotential().size());

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

        if(core_pol && k == 1)
        {
            double R1 = 0.;
            double R2 = 0.;
            for(p=0; p<limit; p++)
            {
                double r2 = R[p]*R[p] + core_rad*core_rad;
                R1 += (s_1->f[p] * s_3->f[p] + Constant::AlphaSquared * s_1->g[p] * s_3->g[p])/r2 * dR[p];
                R2 += density[p]/r2 * dR[p];
            }

            radial -= core_pol * R1 * R2;
        }
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
