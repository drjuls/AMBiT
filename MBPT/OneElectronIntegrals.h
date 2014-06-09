#ifndef ONE_ELECTRON_OPERATOR_H
#define ONE_ELECTRON_OPERATOR_H

#include "HartreeFock/SpinorOperator.h"
#include "Basis/OrbitalManager.h"
#include "ElectronInfo.h"

/** OneElectronIntegrals takes a SpinorOperator and maps it to an electron operator for
    ManyBodyOperator, i.e. it implements GetMatrixElement(const ElectronInfo& e1, ...),
    adding the trivial angular dependence delta(kappa1, kappa2).delta(m1, m2).
    Radial integrals are stored using CalculateOneElectronIntegrals.
 */
class OneElectronIntegrals
{
public:
    OneElectronIntegrals(pOrbitalManagerConst pOrbitals, pSpinorOperatorConst pOperator):
        orbitals(pOrbitals), op(pOperator)
    {
        num_orbitals = orbitals->size();
    }

    /** Calculate one-electron integrals and return number of integrals stored.
        If check_size_only is true, the integrals are not calculated, but the storage size is returned.
     */
    virtual unsigned int CalculateOneElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, bool check_size_only = false);

    /** Clear all integrals. */
    virtual void clear() { integrals.clear(); }

    /** Number of stored integrals. */
    virtual unsigned int size() const { return integrals.size(); }

    inline double GetMatrixElement(const OrbitalInfo& e1, const OrbitalInfo& e2) const
    {
        if(e1.Kappa() == e2.Kappa())
        {
            unsigned int i1 = orbitals->state_index.at(e1);
            unsigned int i2 = orbitals->state_index.at(e2);

            return integrals.at(GetKey(i1, i2));
        }
        else
            return 0.0;
    }

    inline double GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2) const
    {
        if(e1.Kappa() == e2.Kappa() && e1.TwoM() == e2.TwoM())
        {
            unsigned int i1 = orbitals->state_index.at(e1);
            unsigned int i2 = orbitals->state_index.at(e2);

            return integrals.at(GetKey(i1, i2));
        }
        else
            return 0.0;
    }

    /** Read integrals, adding to existing keys or creating new ones. */
    virtual void Read(const std::string& filename);
    virtual void Write(const std::string& filename) const;

    /** Get the spinor operator. */
    pSpinorOperatorConst GetOperator() const { return op; }

protected:
    inline unsigned int GetKey(unsigned int i1, unsigned int i2) const
    {
        if(i1 < i2)
            return i1 * num_orbitals + i2;
        else
            return i2 * num_orbitals + i1;
    }

    pSpinorOperatorConst op;
    pOrbitalManagerConst orbitals;
    unsigned int num_orbitals;
    std::map<unsigned int, double> integrals;
};

typedef boost::shared_ptr<OneElectronIntegrals> pOneElectronIntegrals;
typedef boost::shared_ptr<const OneElectronIntegrals> pOneElectronIntegralsConst;

#endif
