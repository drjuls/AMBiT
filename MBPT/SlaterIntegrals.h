#ifndef SLATER_INTEGRALS_H
#define SLATER_INTEGRALS_H

#include "Basis/OrbitalManager.h"
#include "HartreeFock/HartreeY.h"
#include <map>
#include <unordered_map>

/** Class to hold Slater integrals \f$ R^k(12,34) \f$.
    storage_id is used to store and retrieve integrals in files.
    The existence of reversal symmetry
        \f$ R^k(12,34) = R^k(12,43) = R^k(21,34) \f$
    is specified by two_body_reverse_symmetry, which defaults to false.
    The symmetry is present in the usual Coulomb operator but is broken by
    operators such as specific mass shift as well as MBPT.

    If the size of the map becomes very large you may like to use a
    hash map to store the Slater integrals.
        std::unordered_map is smaller than map, and is fast
        google::sparse_hash_map is very space-efficient, but can take a long time to build
 */
template <class MapType>
class SlaterIntegrals
{
    typedef typename MapType::key_type KeyType;

public:
    SlaterIntegrals(pHartreeY hartreeY_op, pOrbitalManagerConst orbitals, bool two_body_reverse_symmetry_exists = false);
    virtual ~SlaterIntegrals() {}

    /** Calculate two-electron Slater integrals, \f$ R^k(12,34) \f$, and return number of integrals that will be stored.
        The first four arguments are the types of orbitals to store for each limb of the Slater integral
            1 ------- 3
                 |
            2 ------- 4
        Only unique integrals are stored in any case.
        If check_size_only is true, the integrals are not calculated, but the storage size is returned.
     */
    virtual unsigned int CalculateTwoElectronIntegrals(OrbitalClassification s1_type, OrbitalClassification s2_type, OrbitalClassification s3_type, OrbitalClassification s4_type, bool check_size_only = false);

    /** Clear all integrals. */
    void clear() { TwoElectronIntegrals.clear(); }

    /** Number of stored integrals. */
    unsigned int size() const { return TwoElectronIntegrals.size(); }

    /** GetTwoElectronIntegral(k, 1, 2, 3, 4) = R_k(12, 34): 1->3, 2->4
        PRE: s1, s2, s3, and s4 are all in orbital manager and conform to a orbital type pattern that has been
             calculated using CalculateTwoElectronIntegrals.
     */
    virtual double GetTwoElectronIntegral(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4) const;

protected:
    /** Change ordering of states so that it corresponds to a stored integral and return key. */
    virtual KeyType GetKey(unsigned int k, unsigned int i1, unsigned int i2, unsigned int i3, unsigned int i4) const;

    inline void swap(unsigned int& i1, unsigned int& i2) const
    {   unsigned int temp = i1;
        i1 = i2;
        i2 = temp;
    }

protected:
    pHartreeY hartreeY_operator;
    bool two_body_reverse_symmetry;

    pOrbitalManagerConst orbitals;
    KeyType NumStates;

    // TwoElectronIntegrals(k, i, j, l, m) = R_k(ij, lm): i->l, j->m
    MapType TwoElectronIntegrals;
};

typedef SlaterIntegrals<std::map<unsigned long long int, double>> SlaterIntegralsMap;
typedef boost::shared_ptr<SlaterIntegralsMap> pSlaterIntegralsMap;

typedef SlaterIntegrals<std::unordered_map<unsigned long long int, double>> SlaterIntegralsHash;

#if USE_GOOGLE_SPARSEHASH
    #include <google/sparse_hash_map>
    typedef SlaterIntegrals<LongKey, google::sparse_hash_map<LongKey, double>> SlaterIntegralsSparseHash;
#endif

#include "SlaterIntegrals.cpp"

#endif
