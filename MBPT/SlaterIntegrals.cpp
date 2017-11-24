#include "Include.h"

#ifdef AMBIT_USE_OPENMP
    #include <omp.h>
#endif

// Below purposely not included: this file is for a template class and should be included in the header.
// #include "SlaterIntegrals.h"

template <class MapType>
SlaterIntegrals<MapType>::SlaterIntegrals(pOrbitalManagerConst orbitals, bool two_body_reverse_symmetry_exists):
    SlaterIntegralsInterface(orbitals, two_body_reverse_symmetry_exists)
{
    NumStates = orbitals->size();
    SetUpMap();
}

template <class MapType>
SlaterIntegrals<MapType>::SlaterIntegrals(pOrbitalManagerConst orbitals, pHartreeY hartreeY_op, bool two_body_reverse_symmetry_exists):
    SlaterIntegralsInterface(orbitals, two_body_reverse_symmetry_exists), hartreeY_operator(hartreeY_op)
{
    NumStates = orbitals->size();
    SetUpMap();
}

template <class MapType>
SlaterIntegrals<MapType>::SlaterIntegrals(pOrbitalManagerConst orbitals, pHartreeY hartreeY_op):
    SlaterIntegralsInterface(orbitals, hartreeY_op->ReverseSymmetryExists()), hartreeY_operator(hartreeY_op)
{
    NumStates = orbitals->size();
    SetUpMap();
}

template <class MapType>
unsigned int SlaterIntegrals<MapType>::CalculateTwoElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, pOrbitalMapConst orbital_map_3, pOrbitalMapConst orbital_map_4, bool check_size_only)
{
    // NOTE: For each set of orbitals, we actually calculate
    //       R^k(12,34) = < 4 | Y^k_{31} | 2 >
    // on the assumption that i1 and i2 are smaller.

    unsigned int i1, i2, i3, i4;
    int k;
    pOrbitalConst s1, s2, s3, s4;

    std::set<KeyType> found_keys;   // For check_size_only
    if(check_size_only)
        hartreeY_operator->SetLightWeightMode(true);

    // Get Y^k_{31}
    auto it_1 = orbital_map_1->begin();
#ifdef AMBIT_USE_OPENMP
    // The HartreeY operator is not thread-safe, so make a separate clone for each thread
    std::vector<pHartreeY> hartreeY_operators;
    for(int ii = 0; ii < omp_get_max_threads(); ++ii){
        hartreeY_operators.emplace_back(hartreeY_operator->Clone());
    }
    #pragma omp parallel
    {
    #pragma omp single nowait
#endif
    while(it_1 != orbital_map_1->end())
    {
        i1 = orbitals->state_index.at(it_1->first);
        s1 = it_1->second;

        auto it_3 = orbital_map_3->begin();
        if(two_body_reverse_symmetry && orbital_map_1 == orbital_map_3)
        {   it_3 = it_1;
            i3 = i1;
        }

        while(it_3 != orbital_map_3->end())
        {
            i3 = orbitals->state_index.at(it_3->first);
            s3 = it_3->second;

            // Limits on k. This is the expensive part to calculate
#ifdef AMBIT_USE_OPENMP
            #pragma omp task firstprivate(i1, i3, s1, s3) private(k, i2, i4, s2, s4)
            {
            k = hartreeY_operators[omp_get_thread_num()]->SetOrbitals(s3, s1);
#else
            k = hartreeY_operator->SetOrbitals(s3, s1);
#endif
            while(k != -1)
            {
                auto it_2 = orbital_map_2->begin();
                while(it_2 != orbital_map_2->end())
                {
                    i2 = orbitals->state_index.at(it_2->first);
                    s2 = it_2->second;

                    auto it_4 = orbital_map_4->begin();
                    while(it_4 != orbital_map_4->end())
                    {
                        i4 = orbitals->state_index.at(it_4->first);
                        s4 = it_4->second;

                        // Check max_pqn conditions and k conditions
                        if(((s1->L() + s2->L() + s3->L() + s4->L())%2 == 0) &&
                           (2 * k >= abs(s2->TwoJ() - s4->TwoJ())) &&
                           (2 * k <= s2->TwoJ() + s4->TwoJ()))
                        {
                            KeyType key = GetKey(k, i1, i2, i3, i4);

                            if(check_size_only)
                            {
#ifdef AMBIT_USE_OPENMP
                                #pragma omp critical(TWO_ELECTRON_SLATER)
#endif
                                found_keys.insert(key);
                            }
                            else
                            {   // Check that this integral doesn't already exist
                                if(TwoElectronIntegrals.find(key) == TwoElectronIntegrals.end())
                                {
#ifdef AMBIT_USE_OPENMP
                                    double radial = hartreeY_operators[omp_get_thread_num()]->GetMatrixElement(*s4, *s2);
                                    #pragma omp critical(TWO_ELECTRON_SLATER)
#else
                                    double radial = hartreeY_operator->GetMatrixElement(*s4, *s2);
#endif
                                    TwoElectronIntegrals.insert(std::pair<KeyType, double>(key, radial));
                                }
                            }
                        }
                        it_4++;
                    }
                    it_2++;
                }
#ifdef AMBIT_USE_OPENMP
                k = hartreeY_operators[omp_get_thread_num()]->NextK();
            } // K loop
            } // OpenMP task
#else
                k = hartreeY_operator->NextK();
            } // K loop
#endif
            it_3++;
        }
        it_1++;
    }
#ifdef AMBIT_USE_OPENMP
    #pragma omp taskwait
    } // OpenMP parallel region
#endif

    if(check_size_only)
    {   hartreeY_operator->SetLightWeightMode(false);
        return found_keys.size();
    }
    else
        return TwoElectronIntegrals.size();
}

template <class MapType>
auto SlaterIntegrals<MapType>::GetKey(unsigned int k, unsigned int i1, unsigned int i2, unsigned int i3, unsigned int i4) const -> KeyType
{
    if(two_body_reverse_symmetry)
    {   // Ordering of indices:
        // (i1 <= i3) && (i2 <= i4) && (i1 <= i2) && (if i1 == i2, then (i3 <= i4))
        // therefore (i1 <= i2 <= i4) and (i1 <= i3)
        if(i3 < i1)
            std::swap(i3, i1);
        if(i4 < i2)
            std::swap(i4, i2);
        if(i2 < i1)
        {   std::swap(i2, i1);
            std::swap(i3, i4);
        }
        if((i1 == i2) && (i4 < i3))
            std::swap(i3, i4);
    }
    else
    {   // Ordering of indices:
        // i1 is smallest && (if i1 == i2, then (i3 <= i4))
        //                && (if i1 == i3, then (i2 <= i4))
        //                && (if i1 == i4, then (i2 <= i3))

        // Assert one of i1, i3 is smallest
        if(mmin(i1, i3) > mmin(i2, i4))
        {   std::swap(i1, i2);
            std::swap(i3, i4);
        }
        // Assert i1 <= i3
        if(i1 > i3)
        {   std::swap(i1, i3);
            std::swap(i2, i4);
        }

        if((i1 == i2) && (i4 < i3))
            std::swap(i3, i4);
        if((i1 == i3) && (i4 < i2))
            std::swap(i2, i4);
        if((i1 == i4) && (i3 < i2))
            std::swap(i2, i3);
    }

    KeyType key = k  * NumStates*NumStates*NumStates*NumStates +
                  i1 * NumStates*NumStates*NumStates +
                  i2 * NumStates*NumStates +
                  i3 * NumStates +
                  i4;
    return key;
}

template <class MapType>
auto SlaterIntegrals<MapType>::GetKey(ExpandedKeyType expanded_key) const -> KeyType
{
    return GetKey(std::get<0>(expanded_key), std::get<1>(expanded_key), std::get<2>(expanded_key), std::get<3>(expanded_key), std::get<4>(expanded_key));
}

template <class MapType>
double SlaterIntegrals<MapType>::GetTwoElectronIntegral(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4) const
{
    unsigned int i1 = orbitals->state_index.at(s1);
    unsigned int i2 = orbitals->state_index.at(s2);
    unsigned int i3 = orbitals->state_index.at(s3);
    unsigned int i4 = orbitals->state_index.at(s4);

    KeyType key = GetKey(k, i1, i2, i3, i4);
    double radial = 0.;

    auto it = TwoElectronIntegrals.find(key);
    if(it != TwoElectronIntegrals.end())
    {
        radial = it->second;
    }
    else if((s1.L() + s3.L() + k) == 0)
    {   // Only print error if requested integral has correct parity rules
#ifdef AMBIT_USE_OPENMP
        #pragma omp critical(ERRSTREAM)
#endif
        *errstream << "SlaterIntegrals::GetTwoElectronIntegral() failed to find integral."
                   << "\n  R^" << k << " ( " << s1.Name() << " " << s2.Name()
                   << ", " << s3.Name() << " " << s4.Name() << "):  key = "
                   << key << "  num_states = " << NumStates << "\n";
    }

    return radial;
}

template <class MapType>
void SlaterIntegrals<MapType>::Read(const std::string& filename)
{
    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
    {   *errstream << "SlaterIntegrals::Read: file " << filename << " not found." << std::endl;
        return;
    }

    OrbitalIndex old_state_index;
    ReadOrbitalIndexes(old_state_index, fp);

    unsigned int old_key_size;
    fread(&old_key_size, sizeof(unsigned int), 1, fp);

    switch(old_key_size)
    {
        case sizeof(unsigned long long int):
        {
            unsigned int num_integrals;
            unsigned long long int old_key;
            double value;

            fread(&num_integrals, sizeof(unsigned int), 1, fp);

            unsigned long long int old_num_states = old_state_index.size();

            for(unsigned int i = 0; i < num_integrals; i++)
            {
                fread(&old_key, sizeof(unsigned long long int), 1, fp);
                fread(&value, sizeof(double), 1, fp);

                ExpandedKeyType temp_expanded = ReverseKey(old_num_states, old_key);
                KeyType new_key = GetKey(temp_expanded);

                auto it = TwoElectronIntegrals.find(new_key);
                if(it == TwoElectronIntegrals.end())
                    TwoElectronIntegrals[new_key] = value;
                else
                    it->second += value;
            }
            break;
        }

        case sizeof(unsigned int):
        {
            unsigned int num_integrals;
            unsigned int old_key;
            double value;

            fread(&num_integrals, sizeof(unsigned int), 1, fp);

            unsigned long long int old_num_states = old_state_index.size();

            for(unsigned int i = 0; i < num_integrals; i++)
            {
                fread(&old_key, sizeof(unsigned int), 1, fp);
                fread(&value, sizeof(double), 1, fp);

                ExpandedKeyType temp_expanded = ReverseKey(old_num_states, old_key);
                KeyType new_key = GetKey(temp_expanded);

                auto it = TwoElectronIntegrals.find(new_key);
                if(it == TwoElectronIntegrals.end())
                    TwoElectronIntegrals[new_key] = value;
                else
                    it->second += value;
            }
            break;
        }
    }

    fclose(fp);
}

template <class MapType>
auto SlaterIntegrals<MapType>::ReverseKey(KeyType num_states, KeyType key) -> ExpandedKeyType
{
    KeyType running_power = num_states * num_states * num_states * num_states;
    KeyType remainder = key;
    ExpandedKeyType expanded_key;

    std::get<0>(expanded_key) = remainder/running_power;
    remainder -= std::get<0>(expanded_key) * running_power;
    running_power = running_power/num_states;

    std::get<1>(expanded_key) = remainder/running_power;
    remainder -= std::get<1>(expanded_key) * running_power;
    running_power = running_power/num_states;

    std::get<2>(expanded_key) = remainder/running_power;
    remainder -= std::get<2>(expanded_key) * running_power;
    running_power = running_power/num_states;

    std::get<3>(expanded_key) = remainder/running_power;
    remainder -= std::get<3>(expanded_key) * running_power;
    running_power = running_power/num_states;

    std::get<4>(expanded_key) = remainder/running_power;
    remainder -= std::get<4>(expanded_key) * running_power;

    return expanded_key;
}

template <class MapType>
void SlaterIntegrals<MapType>::Write(const std::string& filename) const
{
    if(ProcessorRank == 0)
    {
        FILE* fp = fopen(filename.c_str(), "wb");

        // Write state index
        WriteOrbitalIndexes(orbitals->state_index, fp);

        unsigned int KeyType_size = sizeof(KeyType);
        fwrite(&KeyType_size, sizeof(unsigned int), 1, fp);

        unsigned int num_integrals = size();
        fwrite(&num_integrals, sizeof(unsigned int), 1, fp);

        for(auto& pair: TwoElectronIntegrals)
        {
            const double value = pair.second;   // Convert to double
            fwrite(&pair.first, sizeof(KeyType), 1, fp);
            fwrite(&value, sizeof(double), 1, fp);
        }

        fclose(fp);
    }
}
