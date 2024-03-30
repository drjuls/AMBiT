#ifndef SLATER_INTEGRALS_H
#define SLATER_INTEGRALS_H

#include "Basis/OrbitalManager.h"
#include "HartreeFock/HartreeY.h"
#include <map>
#include <unordered_map>
#include <sparsehash/dense_hash_map>
#include <sparsehash/sparse_hash_map>
#include <absl/container/flat_hash_map.h>

namespace Ambit
{
/** Class to hold Slater integrals \f$ R^k(12,34) \f$.
    storage_id is used to store and retrieve integrals in files.
    The existence of reversal symmetry
        \f$ R^k(12,34) = R^k(12,43) = R^k(21,34) \f$
    is specified by two_body_reverse_symmetry, which defaults to false.
    The symmetry is present in the usual Coulomb operator but is broken by
    operators such as specific mass shift as well as MBPT.
*/
class SlaterIntegralsInterface
{
public:
    SlaterIntegralsInterface(pOrbitalManagerConst orbitals, bool two_body_reverse_symmetry_exists):
        orbitals(orbitals), two_body_reverse_symmetry(two_body_reverse_symmetry_exists)
    {}
    virtual ~SlaterIntegralsInterface() {}

    /** Calculate two-electron Slater integrals, \f$ R^k(12,34) \f$, and return number of integrals that will be stored.
        The first four arguments are the sets of orbitals to store for each limb of the Slater integral
            1 ------- 3
                 |
            2 ------- 4
        All orbitals must be present in orbital manager.
        Only unique integrals are stored in any case.
        If check_size_only is true, the integrals are not calculated, but the storage size is returned.
     */
    virtual unsigned int CalculateTwoElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, pOrbitalMapConst orbital_map_3, pOrbitalMapConst orbital_map_4, bool check_size_only = false) = 0;

    /** Clear all integrals. */
    virtual void clear() = 0;

    /** Number of stored integrals. */
    virtual unsigned int size() const = 0;

    /** Whether this is "reversible" in the sense \f$ R_k(12, 34) = R_k(14, 32) \f$.
     */
    virtual bool ReverseSymmetryExists() const { return two_body_reverse_symmetry; }

    /** Whether any off-parity radial integrals are non-zero, i.e. \f$ (l_1 + l_3 + k) \f$ is odd.
        The total parity must still hold: \f$(l_1 + l_2 + l_3 + l_4)\f$ is even.
     */
    virtual bool OffParityExists() const { return false; }

    /** GetTwoElectronIntegral(k, 1, 2, 3, 4) = R_k(12, 34): 1->3, 2->4
        PRE: s1, s2, s3, and s4 are all in orbital manager and conform to a orbital type pattern that has been
             calculated using CalculateTwoElectronIntegrals.
     */
    virtual double GetTwoElectronIntegral(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4) const = 0;

    /** Structure of files:
        store state_index for all orbitals, and then integrals
        -------------------------------------------------------------
        | size  | index |    value      | index |    value      | ...
        |       |       |   (double)    |       |   (double)    |
        -------------------------------------------------------------
     */
    virtual void Read(const std::string& filename) = 0;
    virtual void Write(const std::string& filename) const = 0;

protected:
    bool two_body_reverse_symmetry;
    pOrbitalManagerConst orbitals;
};

typedef std::shared_ptr<SlaterIntegralsInterface> pSlaterIntegrals;
typedef std::shared_ptr<const SlaterIntegralsInterface> pSlaterIntegralsConst;

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
class SlaterIntegrals : public SlaterIntegralsInterface
{
protected:
    typedef typename MapType::key_type KeyType;
    typedef std::tuple<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int> ExpandedKeyType;   // k, i1, i2, i3, i4

    /** Create without hartreeY operator; CalculateTwoElectronIntegrals() must be overwritten. */
    SlaterIntegrals(pOrbitalManagerConst orbitals, bool two_body_reverse_symmetry_exists);

public:
    SlaterIntegrals(pOrbitalManagerConst orbitals, pHartreeY hartreeY_op, bool two_body_reverse_symmetry_exists);
    SlaterIntegrals(pOrbitalManagerConst orbitals, pHartreeY hartreeY_op);  //!< Reversal symmetry is specified by hartreeY_op.
    virtual ~SlaterIntegrals() {}

    /** Calculate two-electron Slater integrals, \f$ R^k(12,34) \f$, and return number of integrals that will be stored.
        The first four arguments are the types of orbitals to store for each limb of the Slater integral
            1 ------- 3
                 |
            2 ------- 4
        Only unique integrals are stored in any case.
        If check_size_only is true, the integrals are not calculated, but the storage size is returned.
     */
    virtual unsigned int CalculateTwoElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, pOrbitalMapConst orbital_map_3, pOrbitalMapConst orbital_map_4, bool check_size_only = false) override;

    /** Clear all integrals. */
    virtual void clear() override { TwoElectronIntegrals.clear(); }

    /** Number of stored integrals. */
    virtual unsigned int size() const override { return TwoElectronIntegrals.size(); }

    /** Whether any off-parity radial integrals are non-zero. */
    virtual bool OffParityExists() const override { return hartreeY_operator->OffParityExists(); }

    /** GetTwoElectronIntegral(k, 1, 2, 3, 4) = R_k(12, 34): 1->3, 2->4
        PRE: s1, s2, s3, and s4 are all in orbital manager and conform to a orbital type pattern that has been
             calculated using CalculateTwoElectronIntegrals.
     */
    virtual double GetTwoElectronIntegral(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4) const override;

    virtual pHartreeY GetHartreeY() { return hartreeY_operator; }

    /** Read integrals, adding to existing keys or creating new ones. */
    virtual void Read(const std::string& filename) override;
    virtual void Write(const std::string& filename) const override;

protected:
    /** Change ordering of states so that it corresponds to a stored integral and return key. */
    virtual KeyType GetKey(unsigned int k, unsigned int i1, unsigned int i2, unsigned int i3, unsigned int i4) const;
    virtual KeyType GetKey(ExpandedKeyType expanded_key) const;

    ExpandedKeyType ReverseKey(KeyType num_states, KeyType key);

    template<typename U = MapType>
    typename boost::enable_if_c<std::is_same<U, google::dense_hash_map<KeyType, double>>::value, void>::type SetUpMap()
    {   TwoElectronIntegrals.set_empty_key(-1LL);
    }

    template<typename U = MapType>
    typename boost::enable_if_c<!std::is_same<U, google::dense_hash_map<KeyType, double>>::value, void>::type SetUpMap()
    {}

protected:
    pHartreeY hartreeY_operator;
    KeyType NumStates;

    // TwoElectronIntegrals(k, i, j, l, m) = R_k(ij, lm): i->l, j->m
    MapType TwoElectronIntegrals;
};

typedef SlaterIntegrals<std::map<unsigned long long int, double>> SlaterIntegralsMap;
typedef SlaterIntegrals<google::dense_hash_map<unsigned long long int, double>> SlaterIntegralsDenseHash;
typedef SlaterIntegrals<google::sparse_hash_map<unsigned long long int, double>> SlaterIntegralsSparseHash;
typedef SlaterIntegrals<absl::flat_hash_map<unsigned long long int, double>> SlaterIntegralsFlatHash;
}

#include "SlaterIntegrals.cpp"

#endif
