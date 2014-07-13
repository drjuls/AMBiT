#ifndef TRANSITION_INTEGRALS_H
#define TRANSITION_INTEGRALS_H

#include "HartreeFock/SpinorOperator.h"
#include "Basis/OrbitalManager.h"
#include "Configuration/ElectronInfo.h"

/** TransitionIntegrals takes a SpinorMatrixElement and maps it to an electron operator for
    ManyBodyOperator, i.e. it implements GetMatrixElement(const ElectronInfo& e1, ...),
    adding angular momentum projection dependence.
    Radial integrals are stored using CalculateOneElectronIntegrals.
 */
class TransitionIntegrals
{
public:
    TransitionIntegrals(pOrbitalManagerConst pOrbitals, pSpinorMatrixElementConst pOperator):
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
        if(abs(e1.TwoJ() - e2.TwoJ()) < 2 * op->GetMaxK())
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
        double matrix_element = MathConstant::Instance()->Electron3j(e2.TwoJ(), e1.TwoJ(), op->GetMaxK(), e2.TwoM(), -e1.TwoM());

        if(matrix_element)
        {
            if((abs(e1.TwoJ() - e1.TwoM())/2)%2 == 1)
                matrix_element = -matrix_element;

            unsigned int i1 = orbitals->state_index.at(e1);
            unsigned int i2 = orbitals->state_index.at(e2);

            matrix_element *= integrals.at(GetKey(i1, i2));
        }

        return matrix_element;
    }

    /** Get the spinor operator. */
    pSpinorMatrixElementConst GetOperator() const { return op; }

protected:
    inline unsigned int GetKey(unsigned int i1, unsigned int i2) const
    {
        return i1 * num_orbitals + i2;
    }

    pSpinorMatrixElementConst op;
    pOrbitalManagerConst orbitals;
    unsigned int num_orbitals;
    std::map<unsigned int, double> integrals;
};

typedef boost::shared_ptr<TransitionIntegrals> pTransitionIntegrals;
typedef boost::shared_ptr<const TransitionIntegrals> pTransitionIntegralsConst;

#endif
