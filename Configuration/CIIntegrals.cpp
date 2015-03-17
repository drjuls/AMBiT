#include "Include.h"
#include "CIIntegrals.h"
//#include "MBPT/CoreMBPTCalculator.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/MathConstant.h"

inline void swap(unsigned int& i1, unsigned int& i2)
{   unsigned int temp = i1;
    i1 = i2;
    i2 = temp;
}

unsigned int CIIntegrals::GetStorageSize() const
{
    unsigned int num_states = states->size();
    unsigned int size1 = 0;
    unsigned int size2 = 0;

    unsigned int i1, i2, i3, i4;
    unsigned int k, kmax;

    OrbitalMap::const_iterator it_1 = states->begin();
    OrbitalMap::const_iterator it_2 = states->begin();
    OrbitalMap::const_iterator it_3 = states->begin();
    OrbitalMap::const_iterator it_4 = states->begin();
    OrbitalMap::const_iterator it_end = states->end();

    // One electron integrals
    while(it_1 != states->end())
    {
        it_2 = it_1;
        while(it_2 != states->end())
        {
            pOrbitalConst si = it_1->second;
            pOrbitalConst sj = it_2->second;

            // Calculate any remaining one electron integrals
            if(si->Kappa() == sj->Kappa())
                size1++;

            it_2++;
        }
        it_1++;
    }
    *logstream << "Num one-electron integrals: " << size1 << std::endl;

    // Two-electron integrals
    it_2 = states->begin(); i2 = 0;
    while(it_2 != it_end)
    {
        pOrbitalConst s_2 = it_2->second;

        if(two_body_reverse_symmetry)
        {   // i2 <= i4
            it_4 = it_2; i4 = i2;
        }
        else
        {   it_4 = states->begin(); i4 = 0;
        }

        while(it_4 != it_end)
        {
            pOrbitalConst s_4 = it_4->second;

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
                it_1 = states->begin(); i1 = 0;
                while((i1 <= i2) && (i1 <= i4) && (it_1->second->PQN() <= max_pqn_1))
                {
                    pOrbitalConst s_1 = it_1->second;

                    // i1 is smallest && (if i1 == i2, then (i3 <= i4))
                    //                && (if i1 == i3, then (i2 <= i4))
                    //                && (if i1 == i4, then (i2 <= i3))

                    it_3 = it_1; i3 = i1;

                    // if i1 == i3, then (i2 <= i4)
                    if(!two_body_reverse_symmetry && (i2 > i4))
                    {   it_3++; i3++;
                    }
                    // if i1 == i4, then (i2 <= i3)
                    if(!two_body_reverse_symmetry && (i1 == i4))
                    {   it_3 = it_2; i3 = i2;
                    }

                    unsigned int i3_limit;
                    if(i1 == i2)
                        i3_limit = i4;
                    else
                        i3_limit = num_states;

                    while((i3 <= i3_limit) && it_3 != it_end)
                    {
                        pOrbitalConst s_3 = it_3->second;

                        // Check max_pqn conditions and k conditions
                        if(((s_2->PQN() <= max_pqn_2) || (s_3->PQN() <= max_pqn_2)) &&
                           (s_2->PQN() <= max_pqn_3) && (s_3->PQN() <= max_pqn_3) &&
                           ((s_1->L() + s_3->L() + k)%2 == 0) &&
                           (k >= (unsigned int)abs(int(s_1->L()) - int(s_3->L()))) &&
                           (double(k) >= fabs(s_1->J() - s_3->J())) &&
                           (k <= s_1->L() + s_3->L()) &&
                           (double(k) <= s_1->J() + s_3->J()))
                        {
                                size2++;
                        }

                        it_3++; i3++;
                    }
                    it_1++; i1++;
                }
                k+=2;
            }
            it_4++; i4++;
        }
        it_2++; i2++;
    }

    *logstream << "Num two-electron integrals: " << size2 << std::endl;
    return (size1 + size2);
}

void CIIntegrals::clear()
{
    state_index.clear();
    reverse_state_index.clear();
    OneElectronIntegrals.clear();
    TwoElectronIntegrals.clear();
    OverlapIntegrals.clear();

    UpdateStateIndexes();
}

unsigned int CIIntegrals::size() const
{
    return OneElectronIntegrals.size() + TwoElectronIntegrals.size();
}

void CIIntegrals::Update()
{
    clear();

    UpdateOneElectronIntegrals();
    UpdateTwoElectronIntegrals();
}

void CIIntegrals::UpdateStateIndexes()
{
    NumStates = states->size();
    state_index.clear();
    reverse_state_index.clear();

    // Iterate through states, assign in order
    OrbitalMap::const_iterator it_i = states->begin();
    unsigned int i = 0;
    while(it_i != states->end())
    {
        state_index.insert(std::make_pair(it_i->first, i));
        reverse_state_index.insert(std::make_pair(i, it_i->first));

        it_i++; i++;
    }
}

void CIIntegrals::UpdateOneElectronIntegrals()
{
    unsigned int i, j;

    OrbitalMap::const_iterator it_i = states->begin();
    OrbitalMap::const_iterator it_j = states->begin();

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

    pOPIntegrator integrator = one_body_operator->GetOPIntegrator();

    // Get single particle integrals
    it_i = states->begin(); i = 0;
    while(it_i != states->end())
    {
        it_j = it_i; j = i;
        while(it_j != states->end())
        {
            pOrbitalConst si = it_i->second;
            pOrbitalConst sj = it_j->second;

            // Calculate any remaining one electron integrals
            if(si->Kappa() == sj->Kappa())
            {
                unsigned int key = i * NumStates + j;

                // Check that this integral doesn't already exist
                // (it may have been read in)
                if(OneElectronIntegrals.find(key) == OneElectronIntegrals.end())
                {
                    double integral = one_body_operator->GetMatrixElement(*si, *sj);
                    OneElectronIntegrals.insert(std::pair<unsigned int, double>(key, integral));
                }
            }

            // Overlap integrals
            if(si->L() == sj->L())
            {
                double overlap = integrator->GetInnerProduct(*si, *sj);
                OverlapIntegrals.insert(std::pair<unsigned int, double>(i * NumStates + j, overlap));
            }

            it_j++; j++;
        }
        it_i++; i++;
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
    unsigned int i1, i2, i3, i4;
    unsigned int k, kmax;

    OrbitalMap::const_iterator it_1 = states->begin();
    OrbitalMap::const_iterator it_2 = states->begin();
    OrbitalMap::const_iterator it_3 = states->begin();
    OrbitalMap::const_iterator it_4 = states->begin();
    OrbitalMap::const_iterator it_end = states->end();

    // Get 2 -> 4
    it_2 = states->begin(); i2 = 0;
    while(it_2 != it_end)
    {
        pOrbitalConst s_2 = it_2->second;
        it_4 = it_2; i4 = i2;
        while(it_4 != it_end)
        {
            pOrbitalConst s_4 = it_4->second;

            // Limits on k
            k = abs(int(s_2->L()) - int(s_4->L()));
            if(fabs(s_2->J() - s_4->J()) > double(k))
                k += 2;

            kmax = s_2->L() + s_4->L();
            if(s_2->J() + s_4->J() < double(kmax))
                kmax -= 2;

            while(k <= kmax)
            {
                // Get Pot24
                hartreeY_operator->SetParameters(k, *s_2, *s_4);

                // s1 is the smallest
                it_1 = states->begin(); i1 = 0;
                while((i1 <= i2) && (it_1->second->PQN() <= max_pqn_1))
                {
                    pOrbitalConst s_1 = it_1->second;

                    it_3 = it_1; i3 = i1;
                    unsigned int i3_limit;
                    if(i1 == i2)
                        i3_limit = i4;
                    else
                        i3_limit = NumStates;

                    while((i3 <= i3_limit) && it_3 != it_end)
                    {
                        pOrbitalConst s_3 = it_3->second;

                        // Check max_pqn conditions and k conditions
                        if(((s_2->PQN() <= max_pqn_2) || (s_3->PQN() <= max_pqn_2)) &&
                           (s_2->PQN() <= max_pqn_3) && (s_3->PQN() <= max_pqn_3) &&
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
                                double radial = hartreeY_operator->GetMatrixElement(*s_1, *s_3);
                                TwoElectronIntegrals.insert(std::pair<unsigned int, double>(key, radial));
                            }
                        }
                        it_3++; i3++;
                    }

                    // Add case where i2 > i4 (swap i2 and i4)
                    if(!two_body_reverse_symmetry && (i2 != i4))
                    {
                        it_3 = it_1; i3 = i1;

                        // if (i1 == i3) then (i2 <= i4) so no swapped version necessary
                        it_3++; i3++;

                        // i2 and i4 are swapped, therefore with new labels
                        // if (i1 == i2) then (i4 <= i3)
                        if(i1 == i2)
                        {   it_3 = it_4; i3 = i4;
                        }

                        i3_limit = NumStates;

                        while((i3 <= i3_limit) && it_3 != it_end)
                        {
                            pOrbitalConst s_3 = it_3->second;

                            // Check max_pqn conditions and k conditions
                            if(((s_2->PQN() <= max_pqn_2) || (s_3->PQN() <= max_pqn_2)) &&
                               (s_2->PQN() <= max_pqn_3) && (s_3->PQN() <= max_pqn_3) &&
                               ((s_1->L() + s_3->L() + k)%2 == 0) &&
                               (int(k) >= abs(int(s_1->L()) - int(s_3->L()))) &&
                               (double(k) >= fabs(s_1->J() - s_3->J())) &&
                               (k <= s_1->L() + s_3->L()) &&
                               (double(k) <= s_1->J() + s_3->J()))
                            {
                                unsigned int key = k  * NumStates*NumStates*NumStates*NumStates +
                                i1 * NumStates*NumStates*NumStates +
                                i4 * NumStates*NumStates +
                                i3 * NumStates +
                                i2;

                                // Check that this integral doesn't already exist
                                // (it may have been read in)
                                if(TwoElectronIntegrals.find(key) == TwoElectronIntegrals.end())
                                {
                                    double radial = hartreeY_operator->GetMatrixElement(*s_1, *s_3, true);
                                    TwoElectronIntegrals.insert(std::pair<unsigned int, double>(key, radial));
                                }
                            }
                            it_3++; i3++;
                        }
                    }

                    it_1++; i1++;
                }
                k+=2;
            }
            it_4++; i4++;
        }
        it_2++; i2++;
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

    TwoElectronIntegralOrdering(i1, i2, i3, i4);

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
//        if(((k + s1.L() + s3.L())%2 == 1) ||
//             (double(k) < fabs(s1.J() - s3.J())) ||
//             (double(k) > s1.J() + s3.J()) ||
//           ((k + s2.L() + s4.L())%2 == 1) ||
//             (double(k) < fabs(s2.J() - s4.J())) ||
//             (double(k) > s2.J() + s4.J()))
//            return 0.;
//
//        pOrbitalConst s_1 = states->GetState(reverse_state_index.find(i1)->second);
//        pOrbitalConst s_2 = states->GetState(reverse_state_index.find(i2)->second);
//        pOrbitalConst s_3 = states->GetState(reverse_state_index.find(i3)->second);
//        pOrbitalConst s_4 = states->GetState(reverse_state_index.find(i4)->second);
//
//        hartreeY_operator->SetParameters(k, *s_2, *s_4);
//        radial = hartreeY_operator->GetMatrixElement(*s_1, *s_3, false);
        *errstream << "Cannot find integral: k = " << k << "; "
                   << s1.Name() << " " << s2.Name() << " " << s3.Name() << " " << s4.Name() << std::endl;
        exit(2);
    }

    return radial;
}

void CIIntegrals::TwoElectronIntegralOrdering(unsigned int& i1, unsigned int& i2, unsigned int& i3, unsigned int& i4) const
{
    if(two_body_reverse_symmetry)
    {   // Ordering of indices:
        // (i1 <= i3) && (i2 <= i4) && (i1 <= i2) && (if i1 == i2, then (i3 <= i4))
        // therefore (i1 <= i2 <= i4) and (i1 <= i3)
        if(i3 < i1)
            swap(i3, i1);
        if(i4 < i2)
            swap(i4, i2);
        if(i2 < i1)
        {   swap(i2, i1);
            swap(i3, i4);
        }
        if((i1 == i2) && (i4 < i3))
            swap(i3, i4);
    }
    else
    {   // Ordering of indices:
        // i1 is smallest && (if i1 == i2, then (i3 <= i4))
        //                && (if i1 == i3, then (i2 <= i4))
        //                && (if i1 == i4, then (i2 <= i3))

        // Assert one of i1, i3 is smallest
        if(mmin(i1, i3) > mmin(i2, i4))
        {   swap(i1, i2);
            swap(i3, i4);
        }
        // Assert i1 <= i3
        if(i1 > i3)
        {   swap(i1, i3);
            swap(i2, i4);
        }

        if((i1 == i2) && (i4 < i3))
            swap(i3, i4);
        if((i1 == i3) && (i4 < i2))
            swap(i2, i4);
        if((i1 == i4) && (i3 < i2))
            swap(i2, i3);
    }
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
