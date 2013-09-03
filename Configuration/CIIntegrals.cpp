#include "Include.h"
#include "CIIntegrals.h"
#include "HartreeFock/StateIntegrator.h"
#include "Universal/CoulombIntegrator.h"
#include "MBPT/CoreMBPTCalculator.h"
#include "Universal/PhysicalConstant.h"

inline void swap(unsigned int& i1, unsigned int& i2)
{   unsigned int temp = i1;
    i1 = i2;
    i2 = temp;
}

unsigned int CIIntegrals::GetStorageSize() const
{
    unsigned int num_states = states.NumStates();
    unsigned int size1 = 0;
    unsigned int size2 = 0;

    unsigned int i1, i2, i3, i4;
    unsigned int k, kmax;

    ConstStateIterator it_1 = states.GetConstStateIterator();
    ConstStateIterator it_2 = states.GetConstStateIterator();
    ConstStateIterator it_3 = states.GetConstStateIterator();
    ConstStateIterator it_4 = states.GetConstStateIterator();

    // One electron integrals
    it_1.First();
    while(!it_1.AtEnd())
    {
        it_2 = it_1;
        while(!it_2.AtEnd())
        {
            const Orbital* si = it_1.GetState();
            const Orbital* sj = it_2.GetState();

            // Calculate any remaining one electron integrals
            if(si->Kappa() == sj->Kappa())
                size1++;

            it_2.Next();
        }
        it_1.Next();
    }
    *logstream << "Num one-electron integrals: " << size1 << std::endl;

    // Two-electron integrals
    it_2.First(); i2 = 0;
    while(!it_2.AtEnd())
    {
        const Orbital* s_2 = it_2.GetState();
        it_4 = it_2; i4 = i2;
        while(!it_4.AtEnd())
        {
            const Orbital* s_4 = it_4.GetState();

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
                while((i1 <= i2) && (it_1.GetState()->GetPQN() <= max_pqn_1))
                {
                    const Orbital* s_1 = it_1.GetState();

                    it_3 = it_1; i3 = i1;
                    unsigned int i3_limit;
                    if(i1 == i2)
                        i3_limit = i4;
                    else
                        i3_limit = num_states;
                    while((i3 <= i3_limit) && !it_3.AtEnd())
                    {
                        const Orbital* s_3 = it_3.GetState();

                        // Check max_pqn conditions and k conditions
                        if(((s_2->GetPQN() <= max_pqn_2) || (s_3->GetPQN() <= max_pqn_2)) &&
                           (s_2->GetPQN() <= max_pqn_3) && (s_3->GetPQN() <= max_pqn_3) &&
                           ((s_1->L() + s_3->L() + k)%2 == 0) &&
                           (k >= (unsigned int)abs(int(s_1->L()) - int(s_3->L()))) &&
                           (double(k) >= fabs(s_1->J() - s_3->J())) &&
                           (k <= s_1->L() + s_3->L()) &&
                           (double(k) <= s_1->J() + s_3->J()))
                        {
                                size2++;
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

    *logstream << "Num two-electron integrals: " << size2 << std::endl;
    return (size1 + size2);
}

void CIIntegrals::Clear()
{
    state_index.clear();
    reverse_state_index.clear();
    OneElectronIntegrals.clear();
    SMSIntegrals.clear();
    TwoElectronIntegrals.clear();
    OverlapIntegrals.clear();

    UpdateStateIndexes();
}

void CIIntegrals::Update()
{
    Clear();

    UpdateOneElectronIntegrals();
    UpdateTwoElectronIntegrals();
}

void CIIntegrals::UpdateStateIndexes()
{
    NumStates = states.NumStates();
    state_index.clear();
    reverse_state_index.clear();

    ConstStateIterator it_i = states.GetConstStateIterator();
    unsigned int i;

    // Iterate through states, assign in order
    it_i.First(); i = 0;
    while(!it_i.AtEnd())
    {
        state_index.insert(std::pair<OrbitalInfo, unsigned int>(OrbitalInfo(it_i.GetState()), i));
        reverse_state_index.insert(std::pair<unsigned int, OrbitalInfo>(i, OrbitalInfo(it_i.GetState())));

        it_i.Next(); i++;
    }
}

void CIIntegrals::UpdateOneElectronIntegrals()
{
    StateIntegrator SI(states.GetLattice());
    unsigned int i, j;

    ConstStateIterator it_i = states.GetConstStateIterator();
    ConstStateIterator it_j = states.GetConstStateIterator();

    // If fp is NULL, calculate one electron integrals, otherwise read them in.
    std::string file1 = read_id + ".one.int";
    FILE* fp = NULL;

    if(!read_id.empty())
    {   fp = fopen(file1.c_str(), "rb");
        if(fp)
        {   ReadOneElectronIntegrals(fp);
            fclose(fp);
        }
    }

    // Get single particle integrals
    it_i.First(); i = 0;
    while(!it_i.AtEnd())
    {
        it_j = it_i; j = i;
        while(!it_j.AtEnd())
        {
            const Orbital* si = it_i.GetState();
            const Orbital* sj = it_j.GetState();

            // Calculate any remaining one electron integrals
            if(si->Kappa() == sj->Kappa())
            {
                unsigned int key = i * NumStates + j;

                // Check that this integral doesn't already exist
                // (it may have been read in)
                if(OneElectronIntegrals.find(key) == OneElectronIntegrals.end())
                {
                    double integral = SI.HamiltonianMatrixElement(*si, *sj, *states.GetCore());
                    OneElectronIntegrals.insert(std::pair<unsigned int, double>(key, integral));
                }
            }

            // i.pqn <= j.pqn, so calculate using derivative of i instead of j
            SMSIntegrals.insert(std::pair<unsigned int, double>(i* NumStates + j,
                -SI.IsotopeShiftIntegral(*it_j.GetState(), *it_i.GetState())));

            // Overlap integrals
            const SingleParticleWavefunction& p1 = *it_i.GetState();
            const SingleParticleWavefunction& p2 = *it_j.GetState();
            if(p1.L() == p2.L())
            {
                double overlap = 0.;
                const double* dR = states.GetLattice()->dR();
                for(unsigned int x=0; x<mmin(p1.Size(), p2.Size()); x++)
                    overlap += (p1.f[x] * p2.f[x] + p1.g[x] * p2.g[x]) * dR[x];

                OverlapIntegrals.insert(std::pair<unsigned int, double>(i * NumStates + j, overlap));
            }

            it_j.Next(); j++;
        }
        it_i.Next(); i++;
    }
}

void CIIntegrals::UpdateTwoElectronIntegrals()
{
    // Read stored two electron integrals.
    std::string file1 = read_id + ".two.int";
    FILE* fp = NULL;

    if(!read_id.empty())
    {   fp = fopen(file1.c_str(), "rb");
        if(fp)
        {   ReadTwoElectronIntegrals(fp);
            fclose(fp);
        }
    }

    // Calculate any remaining two electron integrals.
    CoulombIntegrator CI(states.GetLattice());
    std::vector<double> density(states.GetCore()->GetHFPotential().size());
    std::vector<double> Pot24(states.GetCore()->GetHFPotential().size());
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
        const Orbital* s_2 = it_2.GetState();
        it_4 = it_2; i4 = i2;
        while(!it_4.AtEnd())
        {
            const Orbital* s_4 = it_4.GetState();

            // Limits on k
            k = abs(int(s_2->L()) - int(s_4->L()));
            if(fabs(s_2->J() - s_4->J()) > double(k))
                k += 2;

            kmax = s_2->L() + s_4->L();
            if(s_2->J() + s_4->J() < double(kmax))
                kmax -= 2;

            // Get density24
            if(k <= kmax)
            {   for(p=0; p<mmin(s_2->Size(), s_4->Size()); p++)
                {   density[p] = s_2->f[p] * s_4->f[p] + s_2->g[p] * s_4->g[p];
                }
            }

            while(k <= kmax)
            {
                // Get Pot24
                CI.FastCoulombIntegrate(density, Pot24, k, mmin(s_2->Size(), s_4->Size()));

                // s1 is the smallest
                it_1.First(); i1 = 0;
                while((i1 <= i2) && (it_1.GetState()->GetPQN() <= max_pqn_1))
                {
                    const Orbital* s_1 = it_1.GetState();

                    it_3 = it_1; i3 = i1;
                    unsigned int i3_limit;
                    if(i1 == i2)
                        i3_limit = i4;
                    else
                        i3_limit = NumStates;
                    while((i3 <= i3_limit) && !it_3.AtEnd())
                    {
                        const Orbital* s_3 = it_3.GetState();

                        // Check max_pqn conditions and k conditions
                        if(((s_2->GetPQN() <= max_pqn_2) || (s_3->GetPQN() <= max_pqn_2)) &&
                           (s_2->GetPQN() <= max_pqn_3) && (s_3->GetPQN() <= max_pqn_3) &&
                           ((s_1->L() + s_3->L() + k)%2 == 0) &&
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

                            // Check that this integral doesn't already exist
                            // (it may have been read in)
                            if(TwoElectronIntegrals.find(key) == TwoElectronIntegrals.end())
                            {
                                double radial = 0.;
                                unsigned int limit = mmin(s_1->Size(), s_3->Size());
                                limit = mmin(limit, Pot24.size());
                                for(p=0; p<limit; p++)
                                {
                                    radial += (s_1->f[p] * s_3->f[p] + s_1->g[p] * s_3->g[p])
                                            * Pot24[p] * dR[p];
                                }

                                if(core_pol && k == 1)
                                {
                                    double R1 = 0.;
                                    double R2 = 0.;
                                    for(p=0; p<limit; p++)
                                    {
                                        double r2 = R[p]*R[p] + core_rad*core_rad;
                                        R1 += (s_1->f[p] * s_3->f[p] + s_1->g[p] * s_3->g[p])/r2 * dR[p];
                                        R2 += density[p]/r2 * dR[p];
                                    }

                                    radial -= core_pol * R1 * R2;
                                }

                                TwoElectronIntegrals.insert(std::pair<unsigned int, double>(key, radial));
                            }
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

double CIIntegrals::GetOneElectronIntegral(const OrbitalInfo& s1, const OrbitalInfo& s2) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;

    if(i1 <= i2)
        return OneElectronIntegrals.find(i1 * NumStates + i2)->second;
    else
        return OneElectronIntegrals.find(i2 * NumStates + i1)->second;
}

double CIIntegrals::GetSMSIntegral(const OrbitalInfo& s1, const OrbitalInfo& s2) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;

    if(i1 <= i2)
        return SMSIntegrals.find(i1 * NumStates + i2)->second;
    else
        return -SMSIntegrals.find(i2 * NumStates + i1)->second;
}

double CIIntegrals::GetOverlapIntegral(const OrbitalInfo& s1, const OrbitalInfo& s2) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;

    if(i1 <= i2)
        return OverlapIntegrals.find(i1 * NumStates + i2)->second;
    else
        return OverlapIntegrals.find(i2 * NumStates + i1)->second;
}

double CIIntegrals::GetTwoElectronIntegral(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4) const
{
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

        const SingleParticleWavefunction* s_1 = states.GetState(reverse_state_index.find(i1)->second);
        const SingleParticleWavefunction* s_2 = states.GetState(reverse_state_index.find(i2)->second);
        const SingleParticleWavefunction* s_3 = states.GetState(reverse_state_index.find(i3)->second);
        const SingleParticleWavefunction* s_4 = states.GetState(reverse_state_index.find(i4)->second);

        unsigned int p;
        CoulombIntegrator CI(states.GetLattice());
        const double* R = states.GetLattice()->R();
        const double* dR = states.GetLattice()->dR();
        const double core_pol = states.GetCore()->GetPolarisability();
        const double core_rad = states.GetCore()->GetClosedShellRadius();

        // Get density24
        std::vector<double> density(mmin(s_2->Size(), s_4->Size()));
        for(p=0; p<density.size(); p++)
        {
            density[p] = s_2->f[p] * s_4->f[p] + s_2->g[p] * s_4->g[p];
        }
        density.resize(states.GetCore()->GetHFPotential().size());

        // Get Pot24
        std::vector<double> Pot24(density.size());
        CI.FastCoulombIntegrate(density, Pot24, k);

        unsigned int limit = mmin(s_1->Size(), s_3->Size());
        limit = mmin(limit, Pot24.size());
        for(p=0; p<limit; p++)
        {
            radial += (s_1->f[p] * s_3->f[p] + s_1->g[p] * s_3->g[p])
                        * Pot24[p] * dR[p];
        }

        if(core_pol && k == 1)
        {
            double R1 = 0.;
            double R2 = 0.;
            for(p=0; p<limit; p++)
            {
                double r2 = R[p]*R[p] + core_rad*core_rad;
                R1 += (s_1->f[p] * s_3->f[p] + s_1->g[p] * s_3->g[p])/r2 * dR[p];
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

bool CIIntegrals::TwoElectronIntegralOrdering(unsigned int& i1, unsigned int& i2, unsigned int& i3, unsigned int& i4) const
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

void CIIntegrals::SetIdentifier(const std::string& storage_id)
{
    read_id = storage_id;

    if((NumProcessors > 1) && !storage_id.empty())
    {   std::stringstream proc;
        proc << ProcessorRank;
        write_id = storage_id + '_' + proc.str();
    }
    else
        write_id = storage_id;
}

void CIIntegrals::ReadOneElectronIntegrals(FILE* fp)
{
    // Make state index for stored integrals
    // Get start and end principal quantum numbers for each L
    std::map<unsigned int, OrbitalInfo> stored_state_index;
    unsigned int max_l, L;
    fread(&max_l, sizeof(unsigned int), 1, fp);
    
    unsigned int* start_pqn = new unsigned int[max_l+1];
    unsigned int* end_pqn = new unsigned int[max_l+1];
    unsigned int min_pqn = 100, max_pqn = 0;

    for(L = 0; L <= max_l; L++)
    {   fread(start_pqn+L, sizeof(unsigned int), 1, fp);
        fread(end_pqn+L, sizeof(unsigned int), 1, fp);

        if(start_pqn[L] < min_pqn)
            min_pqn = start_pqn[L];
        if(end_pqn[L] > max_pqn)
            max_pqn = end_pqn[L];
    }

    // States are ordered first by pqn, then L, then Kappa
    unsigned int index = 0;
    for(unsigned int pqn = min_pqn; pqn <= max_pqn; pqn++)
    {
        for(L = 0; L <= max_l; L++)
        {   if((pqn >= start_pqn[L]) && (pqn <= end_pqn[L]))
            {
                if(L != 0)
                {   stored_state_index.insert(std::pair<unsigned int, OrbitalInfo>(index, OrbitalInfo(pqn, int(L))));
                    index++;
                }
                stored_state_index.insert(std::pair<unsigned int, OrbitalInfo>(index, OrbitalInfo(pqn, -int(L+1))));
                index++;
            }
        }
    }

    delete[] start_pqn;
    delete[] end_pqn;

    unsigned int NumStoredStates = stored_state_index.size();

    // Get integrals
    unsigned int size;
    fread(&size, sizeof(unsigned int), 1, fp);

    unsigned int i, stored_key, new_key;
    unsigned int stored_s1, stored_s2;
    unsigned int i1, i2;
    std::map<OrbitalInfo, unsigned int>::const_iterator it;
    double value;

    for(i=0; i<size; i++)
    {
        fread(&stored_key, sizeof(unsigned int), 1, fp);
        fread(&value, sizeof(double), 1, fp);

        // Get stored states
        stored_s2 = stored_key%NumStoredStates;
            stored_key = stored_key/NumStoredStates;
        stored_s1 = stored_key%NumStoredStates;

        // Get new state_indexes. Check each state exists in new basis.
        it = state_index.find(stored_state_index.find(stored_s1)->second);
        if(it != state_index.end())
        {   i1 = it->second;

            it = state_index.find(stored_state_index.find(stored_s2)->second);
            if(it != state_index.end())
            {   i2 = it->second;

                // Check ordering
                if(i1 > i2)
                    swap(i1, i2);

                new_key = i1 * NumStates + i2;

                OneElectronIntegrals.insert(std::pair<unsigned int, double>(new_key, value));
            }
        }
    }
}

void CIIntegrals::ReadTwoElectronIntegrals(FILE* fp)
{
    // Make state index for stored integrals
    // Get start and end principal quantum numbers for each L
    std::map<unsigned int, OrbitalInfo> stored_state_index;
    unsigned int max_l, L;
    fread(&max_l, sizeof(unsigned int), 1, fp);
    
    unsigned int* start_pqn = new unsigned int[max_l+1];
    unsigned int* end_pqn = new unsigned int[max_l+1];
    unsigned int min_pqn = 100, max_pqn = 0;

    for(L = 0; L <= max_l; L++)
    {   fread(start_pqn+L, sizeof(unsigned int), 1, fp);
        fread(end_pqn+L, sizeof(unsigned int), 1, fp);

        if(start_pqn[L] < min_pqn)
            min_pqn = start_pqn[L];
        if(end_pqn[L] > max_pqn)
            max_pqn = end_pqn[L];
    }

    // States are ordered first by pqn, then L, then Kappa
    unsigned int index = 0;
    for(unsigned int pqn = min_pqn; pqn <= max_pqn; pqn++)
    {
        for(L = 0; L <= max_l; L++)
        {   if((pqn >= start_pqn[L]) && (pqn <= end_pqn[L]))
            {
                if(L != 0)
                {   stored_state_index.insert(std::pair<unsigned int, OrbitalInfo>(index, OrbitalInfo(pqn, int(L))));
                    index++;
                }
                stored_state_index.insert(std::pair<unsigned int, OrbitalInfo>(index, OrbitalInfo(pqn, -int(L+1))));
                index++;
            }
        }
    }

    delete[] start_pqn;
    delete[] end_pqn;

    unsigned int NumStoredStates = stored_state_index.size();

    // Get integrals
    unsigned int size;
    fread(&size, sizeof(unsigned int), 1, fp);

    unsigned int i, stored_key, new_key;
    unsigned int k, stored_s1, stored_s2, stored_s3, stored_s4;
    unsigned int i1, i2, i3, i4;
    std::map<OrbitalInfo, unsigned int>::const_iterator it;
    double value;

    for(i=0; i<size; i++)
    {
        fread(&stored_key, sizeof(unsigned int), 1, fp);
        fread(&value, sizeof(double), 1, fp);

        // Get k and stored states
        stored_s4 = stored_key%NumStoredStates;
            stored_key = stored_key/NumStoredStates;
        stored_s3 = stored_key%NumStoredStates;
            stored_key = stored_key/NumStoredStates;
        stored_s2 = stored_key%NumStoredStates;
            stored_key = stored_key/NumStoredStates;
        stored_s1 = stored_key%NumStoredStates;
            stored_key = stored_key/NumStoredStates;
        k = stored_key%NumStoredStates;

        // Get new state_indexes. Check each state exists in new basis.
        it = state_index.find(stored_state_index.find(stored_s1)->second);
        if(it != state_index.end())
        {   i1 = it->second;
            
            it = state_index.find(stored_state_index.find(stored_s2)->second);
            if(it != state_index.end())
            {   i2 = it->second;
                
                it = state_index.find(stored_state_index.find(stored_s3)->second);
                if(it != state_index.end())
                {   i3 = it->second;
                    
                    it = state_index.find(stored_state_index.find(stored_s4)->second);
                    if(it != state_index.end())
                    {   i4 = it->second;
                        TwoElectronIntegralOrdering(i1, i2, i3, i4);
                    
                        new_key = k  * NumStates*NumStates*NumStates*NumStates +
                                  i1 * NumStates*NumStates*NumStates +
                                  i2 * NumStates*NumStates +
                                  i3 * NumStates +
                                  i4;                                

                        TwoElectronIntegrals.insert(std::pair<unsigned int, double>(new_key, value));
                    }
                }
            }
        }
    }
}

void CIIntegrals::WriteOneElectronIntegrals(bool use_read_id) const
{
    std::string filename;
    if(use_read_id)
        filename = read_id + ".one.int";
    else
        filename = write_id + ".one.int";
    FILE* fp = fopen(filename.c_str(), "wb");

    if(fp)
    {   // Store information about the basis
        // Get number of states in each wave
        unsigned int max_l = 0;
        std::map<OrbitalInfo, unsigned int>::const_iterator it = state_index.begin();
        while(it != state_index.end())
        {   if(it->first.L() > max_l)
                max_l = it->first.L();
            it++;
        }

        unsigned int L;
        unsigned int* start_pqn = new unsigned int[max_l+1];
        unsigned int* end_pqn = new unsigned int[max_l+1];
        for(L = 0; L <= max_l; L++)
        {   start_pqn[L] = 100;
            end_pqn[L] = 0;
        }

        it = state_index.begin();
        while(it != state_index.end())
        {   if(it->first.PQN() < start_pqn[it->first.L()])
                start_pqn[it->first.L()] = it->first.PQN();
            if(it->first.PQN() > end_pqn[it->first.L()])
                end_pqn[it->first.L()] = it->first.PQN();
            it++;
        }

        fwrite(&max_l, sizeof(unsigned int), 1, fp);
        for(L = 0; L <= max_l; L++)
        {   fwrite(start_pqn+L, sizeof(unsigned int), 1, fp);
            fwrite(end_pqn+L, sizeof(unsigned int), 1, fp);
        }

        delete[] start_pqn;
        delete[] end_pqn;

        // Store integrals
        unsigned int size = OneElectronIntegrals.size();
        fwrite(&size, sizeof(unsigned int), 1, fp);

        std::map<unsigned int, double>::const_iterator it2 = OneElectronIntegrals.begin();
        while(it2 != OneElectronIntegrals.end())
        {
            fwrite(&(it2->first), sizeof(unsigned int), 1, fp);
            fwrite(&(it2->second), sizeof(double), 1, fp);

            it2++;
        }

        fclose(fp);
    }
}

void CIIntegrals::WriteTwoElectronIntegrals(bool use_read_id) const
{
    std::string filename;
    if(use_read_id)
        filename = read_id + ".two.int";
    else
        filename = write_id + ".two.int";
    FILE* fp = fopen(filename.c_str(), "wb");

    if(fp)
    {   // Store information about the basis
        // Get number of states in each wave
        unsigned int max_l = 0;
        std::map<OrbitalInfo, unsigned int>::const_iterator it = state_index.begin();
        while(it != state_index.end())
        {   if(it->first.L() > max_l)
                max_l = it->first.L();
            it++;
        }

        unsigned int L;
        unsigned int* start_pqn = new unsigned int[max_l+1];
        unsigned int* end_pqn = new unsigned int[max_l+1];
        for(L = 0; L <= max_l; L++)
        {   start_pqn[L] = 100;
            end_pqn[L] = 0;
        }

        it = state_index.begin();
        while(it != state_index.end())
        {   if(it->first.PQN() < start_pqn[it->first.L()])
                start_pqn[it->first.L()] = it->first.PQN();
            if(it->first.PQN() > end_pqn[it->first.L()])
                end_pqn[it->first.L()] = it->first.PQN();
            it++;
        }

        fwrite(&max_l, sizeof(unsigned int), 1, fp);
        for(L = 0; L <= max_l; L++)
        {   fwrite(start_pqn+L, sizeof(unsigned int), 1, fp);
            fwrite(end_pqn+L, sizeof(unsigned int), 1, fp);
        }

        delete[] start_pqn;
        delete[] end_pqn;

        // Store integrals
        unsigned int size = TwoElectronIntegrals.size();
        fwrite(&size, sizeof(unsigned int), 1, fp);

        std::map<unsigned int, double>::const_iterator it2 = TwoElectronIntegrals.begin();
        while(it2 != TwoElectronIntegrals.end())
        {
            fwrite(&(it2->first), sizeof(unsigned int), 1, fp);
            fwrite(&(it2->second), sizeof(double), 1, fp);

            it2++;
        }

        fclose(fp);
    }
}
