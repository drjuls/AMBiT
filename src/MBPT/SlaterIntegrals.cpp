#include "Basis/OrbitalManager.h"
#include "HartreeFock/Orbital.h"
#include "HartreeFock/OrbitalInfo.h"
#include "Include.h"
#include <absl/container/btree_map.h>
#include <tuple>
#include <utility>

#ifdef AMBIT_USE_OPENMP
    #include <omp.h>
#endif

// Below purposely not included: this file is for a template class and should be included in the header.
// #include "SlaterIntegrals.h"

namespace Ambit
{
template <class MapType>
SlaterIntegrals<MapType>::SlaterIntegrals(pOrbitalManagerConst orbitals, bool two_body_reverse_symmetry_exists):
    SlaterIntegralsInterface(orbitals, two_body_reverse_symmetry_exists)
{
    NumStates = orbitals->size();
}

template <class MapType>
SlaterIntegrals<MapType>::SlaterIntegrals(pOrbitalManagerConst orbitals, pHartreeY hartreeY_op, bool two_body_reverse_symmetry_exists):
    SlaterIntegralsInterface(orbitals, two_body_reverse_symmetry_exists), hartreeY_operator(hartreeY_op)
{
    NumStates = orbitals->size();
}

template <class MapType>
SlaterIntegrals<MapType>::SlaterIntegrals(pOrbitalManagerConst orbitals, pHartreeY hartreeY_op):
    SlaterIntegralsInterface(orbitals, hartreeY_op->ReverseSymmetryExists()), hartreeY_operator(hartreeY_op)
{
    NumStates = orbitals->size();
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

    /* First, run over the orbital indices and calculate some de-duplicated set of valid
     * integrals (as well as eliminating integrals which are zero by symmetry
     * requirments), then store their "expanded keys" (a.k.a. a tuple of orbital indices
     * plus k). Do not calculate any integrals at this point - that comes later. The
     * reason for splitting this step out and doing it in serial is because we store the
     * calculated integrals in a shared hash table and we want to avoid lock contention
     * between lots of OpenMP threads. 
     *
     * We therefore calculate which integrals need to be calculated (represented by
     * orbital indices in a tuple (k, i1, i2, i3, i4)) and essentially pre-bake a set
     * of tasks, which can later be executed in parallel (via OpenMP task parallelism).
     *
     * The only wrinkle is that we have a weak set of ordering constraints, due to the
     * fact that we can re-use some of the results between integrals with shared
     * (k, i1, i3), as implemented in the HartreeY (Hartree Screening) operator. As a
     * result, we want to group all integrals with common (k, i1, i3) and ensure they
     * are all calculated in a block to reuse the common HartreeY components. Currently,
     * this we do this via an ordered map (using Abseil's B-Tree implementation) where
     * we map (k, i1, i3) tuples to a list of all the (i2, i4) orbital indices such that
     * (k, i1, i2, i3, i4) is a nonzero and unique Slater integral. This can then be
     * iterated over in parallel, with each (k, i1, i3) key corresponding to a single
     * OpenMP task.
     */
    absl::btree_map<std::tuple<int, unsigned, unsigned>, std::list<std::pair<unsigned, unsigned>>> iteration_tree;
    std::vector<double> values;

    absl::flat_hash_set<KeyType> found_keys;   // For check_size_only
    hartreeY_operator->SetLightWeightMode(true);

    for(auto it_1 = orbital_map_1->begin(); it_1 < orbital_map_1->end(); it_1++)
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

            // Limits on k. Only doing this in "lightweight" mode for now
            k = hartreeY_operator->SetOrbitals(s3, s1);
            while(k != -1)
            {
                auto key_tuple = std::make_tuple(k, i1, i3);
                iteration_tree[key_tuple] = std::list<std::pair<unsigned, unsigned> >{0};

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

                            // Check if we actually inserted a new key into the set
                            bool inserted = found_keys.insert(key).second;
                            if(inserted)
                            {
                                // Store the expanded key *in-order* so we can amortise
                                // shared integral elements
                                iteration_tree[key_tuple].push_back(std::make_pair(i2, i4));
                                // Also insert a placeholder entry in the main hashtable
                                TwoElectronIntegrals.insert(std::pair<KeyType, double>(key, 0.0));
                            }
                        }
                        it_4++;
                    }
                    it_2++;
                }
                k = hartreeY_operator->NextK();
            } // K loop
            it_3++;
        }
    }

    hartreeY_operator->SetLightWeightMode(false);
    size_t num_integrals = found_keys.size();

    // If we're only checking the number of integrals, then we're done so return early
    if(check_size_only)
        return(num_integrals);

    // Now actually calculate the integrals in parallel, via OpenMP tasks
#ifdef AMBIT_USE_OPENMP
    // The HartreeY operator is not thread-safe, so make a separate clone for each thread
    std::vector<pHartreeY> hartreeY_operators;
    for(int ii = 0; ii < omp_get_max_threads(); ++ii){
        hartreeY_operators.emplace_back(hartreeY_operator->Clone());
    }
#pragma omp parallel
{
#pragma omp single nowait
    {
#endif
    for(auto it = iteration_tree.begin(); it != iteration_tree.end(); it++)
    {
#ifdef AMBIT_USE_OPENMP
#pragma omp task default(none) \
                 firstprivate(it) \
                 shared(hartreeY_operators, orbital_map_1, orbital_map_2, \
                        orbital_map_3, orbital_map_4, orbitals, TwoElectronIntegrals)
        {
#endif
        auto key = it->first;
        int k = std::get<0>(key);
        unsigned i1 = std::get<1>(key);
        unsigned i3 = std::get<2>(key);

        OrbitalInfo s1 = orbitals->reverse_state_index.at(i1);
        OrbitalInfo s3 = orbitals->reverse_state_index.at(i3);

        pOrbitalConst orb1 = orbital_map_1->GetState(s1);
        pOrbitalConst orb3 = orbital_map_3->GetState(s3);

#ifdef AMBIT_USE_OPENMP
        hartreeY_operators[omp_get_thread_num()]->SetParameters(k, orb3, orb1);
#else
        hartreeY_operator->SetParameters(k, orb3, orb1);
#endif

        auto list = it->second;
        for(auto item : list)
        {
            unsigned i2 = item.first;
            unsigned i4 = item.second;
            
            OrbitalInfo s2 = orbitals->reverse_state_index.at(i2);
            OrbitalInfo s4 = orbitals->reverse_state_index.at(i4);

            pOrbitalConst orb2 = orbital_map_2->GetState(s2);
            pOrbitalConst orb4 = orbital_map_4->GetState(s4);

#ifdef AMBIT_USE_OPENMP
            double radial = hartreeY_operators[omp_get_thread_num()]->GetMatrixElement(*orb4, *orb2);
#else
            double radial = hartreeY_operator->GetMatrixElement(*orb4, *orb2);
#endif
            KeyType key = GetKey(k, i1, i2, i3, i4);
            TwoElectronIntegrals[key] = radial;
        }
    }
#ifdef AMBIT_USE_OPENMP
        } // omp task
    } // omp single
} // omp parallel
#endif
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
    else if((s1.L() + s3.L() + k)%2 == 0 && (s2.L() + s4.L() + k)%2 == 0)
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
bool SlaterIntegrals<MapType>::Read(const std::string& filename)
{
    FILE* fp = file_err_handler->fopen(filename.c_str(), "rb");
    if(!fp)
    {
        return false;
    }

    OrbitalIndex old_state_index;
    ReadOrbitalIndexes(old_state_index, fp);
    unsigned long long int old_num_states = old_state_index.size();
    ReverseOrbitalIndex rev_old_state_index = GetReverseIndex(old_state_index);

    unsigned int old_key_size;
    unsigned int num_integrals;
    double value;
    file_err_handler->fread(&old_key_size, sizeof(unsigned int), 1, fp);
    file_err_handler->fread(&num_integrals, sizeof(unsigned int), 1, fp);

    unsigned int old_key_4;
    unsigned long long int old_key_8;
    ExpandedKeyType temp_expanded;

    OrbitalInfo s1, s2, s3, s4;
    OrbitalIndex::const_iterator i1, i2, i3, i4;

    for(unsigned int i = 0; i < num_integrals; i++)
    {
        if(old_key_size == sizeof(unsigned long long int))
        {
            file_err_handler->fread(&old_key_8, sizeof(unsigned long long int), 1, fp);
            temp_expanded = ReverseKey(old_num_states, old_key_8);
        }
        else
        {
            file_err_handler->fread(&old_key_4, sizeof(unsigned int), 1, fp);
            temp_expanded = ReverseKey(old_num_states, old_key_4);
        }

        file_err_handler->fread(&value, sizeof(double), 1, fp);

        s1 = rev_old_state_index.at(std::get<1>(temp_expanded));
        s2 = rev_old_state_index.at(std::get<2>(temp_expanded));
        s3 = rev_old_state_index.at(std::get<3>(temp_expanded));
        s4 = rev_old_state_index.at(std::get<4>(temp_expanded));

        // Old valence states no longer being included is not an error
        i1 = orbitals->state_index.find(s1);
        i2 = orbitals->state_index.find(s2);
        i3 = orbitals->state_index.find(s3);
        i4 = orbitals->state_index.find(s4);

        if(i1 != orbitals->state_index.end() && i2 != orbitals->state_index.end() &&
           i3 != orbitals->state_index.end() && i4 != orbitals->state_index.end())
        {
            KeyType new_key = GetKey(std::get<0>(temp_expanded), i1->second, i2->second, i3->second, i4->second);

            auto it = TwoElectronIntegrals.find(new_key);
            if(it == TwoElectronIntegrals.end())
                TwoElectronIntegrals[new_key] = value;
            else
                it->second += value;
        }
    }

    file_err_handler->fclose(fp);
    return true;
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

    std::get<4>(expanded_key) = remainder;

    return expanded_key;
}

template <class MapType>
void SlaterIntegrals<MapType>::Write(const std::string& filename) const
{
    if(ProcessorRank == 0)
    {
        FILE* fp = file_err_handler->fopen(filename.c_str(), "wb");

        // Write state index
        WriteOrbitalIndexes(orbitals->state_index, fp);

        unsigned int KeyType_size = sizeof(KeyType);
        file_err_handler->fwrite(&KeyType_size, sizeof(unsigned int), 1, fp);

        unsigned int num_integrals = size();
        file_err_handler->fwrite(&num_integrals, sizeof(unsigned int), 1, fp);

        for(auto& pair: TwoElectronIntegrals)
        {
            const double value = pair.second;   // Convert to double
            file_err_handler->fwrite(&pair.first, sizeof(KeyType), 1, fp);
            file_err_handler->fwrite(&value, sizeof(double), 1, fp);
        }

        file_err_handler->fclose(fp);
    }
}
}
