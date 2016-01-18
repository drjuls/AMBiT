#include "TransitionIntegrals.h"

unsigned int TransitionIntegrals::CalculateOneElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, bool check_size_only)
{
    unsigned int i1, i2;
    pOrbitalConst s1, s2;

    std::set<unsigned int> found_keys;   // For check_size_only

    auto it_1 = orbital_map_1->begin();
    while(it_1 != orbital_map_1->end())
    {
        i1 = orbitals->state_index.at(it_1->first);
        s1 = it_1->second;

        auto it_2 = orbital_map_2->begin();
        while(it_2 != orbital_map_2->end())
        {
            if(abs(it_1->first.TwoJ() - it_2->first.TwoJ()) <= 2 * op->GetK())
            {
                i2 = orbitals->state_index.at(it_2->first);
                s2 = it_2->second;

                unsigned int key = GetKey(i1, i2);
                if(check_size_only)
                    found_keys.insert(key);
                else if(!integrals.count(key))
                    integrals[key] = op->GetReducedMatrixElement(*s1, *s2);
            }
            it_2++;
        }
        it_1++;
    }

    if(check_size_only)
        return found_keys.size();
    else
        return integrals.size();
}
