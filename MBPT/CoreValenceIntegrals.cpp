#include "Include.h"

// Below purposely not included: this file is for a template class and should be included in the header.
// #include "CoreValenceIntegrals.h"

template <class MapType>
CoreValenceIntegrals<MapType>::CoreValenceIntegrals(pOrbitalManagerConst orbitals, pHFElectronOperatorConst one_body, pHartreeY hartreeY_op):
    CoreValenceIntegrals(one_body, pSlaterIntegrals(new SlaterIntegrals<MapType>(orbitals, hartreeY_op)))
{}

template <class MapType>
CoreValenceIntegrals<MapType>::CoreValenceIntegrals(pOrbitalManagerConst orbitals, pHFElectronOperatorConst one_body, pSlaterIntegrals bare_integrals):
    SlaterIntegrals<MapType>(orbitals, false), core_PT(nullptr),
    include_core(true), include_core_subtraction(true), include_core_extra_box(true),
    include_valence(false), include_valence_subtraction(false), include_valence_extra_box(false)
{
    core_PT.reset(new CoreMBPTCalculator(this->orbitals, one_body, bare_integrals));
}

template <class MapType>
CoreValenceIntegrals<MapType>::CoreValenceIntegrals(pOrbitalManagerConst orbitals, pCoreMBPTCalculator core_mbpt_calculator):
    SlaterIntegrals<MapType>(orbitals, false), core_PT(core_mbpt_calculator),
    include_core(true), include_core_subtraction(true), include_core_extra_box(true),
    include_valence(false), include_valence_subtraction(false), include_valence_extra_box(false)
{}

template <class MapType>
CoreValenceIntegrals<MapType>::~CoreValenceIntegrals()
{}

template <class MapType>
unsigned int CoreValenceIntegrals<MapType>::CalculateTwoElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, pOrbitalMapConst orbital_map_3, pOrbitalMapConst orbital_map_4, bool check_size_only)
{
    unsigned int i1, i2, i3, i4;
    int k, kmax;

    std::set<KeyType> found_keys;   // For check_size_only

    if(!check_size_only)
        core_PT->UpdateIntegrals();

    auto it_1 = orbital_map_1->begin();
    while(it_1 != orbital_map_1->end())
    {
        i1 = this->orbitals->state_index.at(it_1->first);
        const auto& s1 = it_1->first;

        auto it_3 = orbital_map_3->begin();
        while(it_3 != orbital_map_3->end())
        {
            i3 = this->orbitals->state_index.at(it_3->first);
            const auto& s3 = it_3->first;

            auto it_2 = orbital_map_2->begin();
            while(it_2 != orbital_map_2->end())
            {
                i2 = this->orbitals->state_index.at(it_2->first);
                const auto& s2 = it_2->first;

                auto it_4 = orbital_map_4->begin();
                while(it_4 != orbital_map_4->end())
                {
                    i4 = this->orbitals->state_index.at(it_4->first);
                    const auto& s4 = it_4->first;

                    // Check parity conservation
                    if((s1.L() + s2.L() + s3.L() + s4.L())%2 == 0)
                    {
                        // Limits on k
                        k = mmax(abs(s1.L() - s3.L()), abs(s2.L() - s4.L()));
                        if((abs(s1.TwoJ() - s3.TwoJ()) > 2 * k) || abs(s2.TwoJ() - s4.TwoJ()) > 2 * k)
                        {
                            if(include_core_extra_box || include_valence_extra_box)
                                k++;
                            else
                                k+=2;
                        }

                        kmax = mmin(s1.L() + s3.L(), s2.L() + s4.L());
                        if((s1.TwoJ() + s3.TwoJ() <  2 * kmax) || (s2.TwoJ() + s4.TwoJ() <  2 * kmax))
                        {
                            if(include_core_extra_box || include_valence_extra_box)
                                kmax--;
                            else
                                kmax-=2;
                        }
                            
                        while(k <= kmax)
                        {
                            // Usual multipolarity rules
                            if((s2.L() + s4.L() + k)%2 == 0)
                            {
                                KeyType key = this->GetKey(k, i1, i2, i3, i4);

                                if(check_size_only)
                                {
                                    if(include_core || include_core_subtraction || include_valence || include_valence_subtraction)
                                        found_keys.insert(key);
                                }
                                else
                                {   // Check that this integral doesn't already exist
                                    if(this->TwoElectronIntegrals.find(key) == this->TwoElectronIntegrals.end())
                                    {
                                        double radial = 0;
                                        if(include_core)
                                            radial += core_PT->GetTwoElectronDiagrams(k, s1, s2, s3, s4);
                                        if(include_core_subtraction)
                                            radial += core_PT->GetTwoElectronSubtraction(k, s1, s2, s3, s4);

                                        this->TwoElectronIntegrals.insert(std::pair<KeyType, double>(key, radial));
                                    }
                                }
                            }
                            // Wrong multipolarity
                            else if(include_core_extra_box || include_valence_extra_box)
                            {
                                KeyType key = this->GetKey(k, i1, i2, i3, i4);

                                if(check_size_only)
                                {
                                    found_keys.insert(key);
                                }
                                else
                                {   // Check that this integral doesn't already exist
                                    if(this->TwoElectronIntegrals.find(key) == this->TwoElectronIntegrals.end())
                                    {
                                        double radial = 0;
                                        if(include_core_extra_box)
                                            radial += core_PT->GetTwoElectronBoxDiagrams(k, s1, s2, s3, s4);

                                        this->TwoElectronIntegrals.insert(std::pair<KeyType, double>(key, radial));
                                    }
                                }
                            }

                            if(include_core_extra_box || include_valence_extra_box)
                                k++;
                            else
                                k+=2;
                        }
                    }
                    it_4++;
                }
                it_2++;
            }
            it_3++;
        }
        it_1++;
    }
    
    if(check_size_only)
        return found_keys.size();
    else
        return this->TwoElectronIntegrals.size();
}

template <class MapType>
void CoreValenceIntegrals<MapType>::IncludeCore(bool include_mbpt, bool include_subtraction, bool include_wrong_parity_box_diagrams)
{
    include_core = include_mbpt;
    include_core_subtraction = include_subtraction;
    include_core_extra_box = include_wrong_parity_box_diagrams;
}

template <class MapType>
void CoreValenceIntegrals<MapType>::IncludeValence(bool include_mbpt, bool include_subtraction, bool include_wrong_parity_box_diagrams)
{
    include_valence = include_mbpt;
    include_valence_subtraction = include_subtraction;
    include_valence_extra_box = include_wrong_parity_box_diagrams;
}
