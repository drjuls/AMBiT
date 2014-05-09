#ifndef ONE_ELECTRON_OPERATOR_H
#define ONE_ELECTRON_OPERATOR_H

#include "HartreeFock/SpinorOperator.h"
#include "Basis/OrbitalManager.h"
#include "ElectronInfo.h"

/** OneElectronOperator takes a SpinorOperator and maps it to an electron operator for
    ManyBodyOperator, i.e. it implements GetMatrixElement(const ElectronInfo& e1, ...),
    adding the trivial angular dependence delta(kappa1, kappa2).delta(m1, m2).
    Initialise with an OrbitalManager to store radial integrals.
 */
template<class pSpinorOperatorType>
class OneElectronOperator
{
public:
    OneElectronOperator(pSpinorOperatorType pOperator, pOrbitalManagerConst pOrbitals):
        op(pOperator), orbitals(pOrbitals)
    {
        num_orbitals = orbitals->size();
        auto it = orbitals->state_index.begin();
        auto end = orbitals->state_index.end();

        const OrbitalMap& all_orbitals = *orbitals->all;
        while(it != end)
        {
            auto jt = it;
            while(jt != end)
            {
                if(it->first.Kappa() == jt->first.Kappa())
                {
                    const Orbital& left(*all_orbitals.GetState(it->first));
                    const Orbital& right(*all_orbitals.GetState(jt->first));
                    integrals[it->second * num_orbitals + jt->second] = op->GetMatrixElement(left, right);
                }
                jt++;
            }
            it++;
        }
    }

    inline double GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2) const
    {
        if(e1.Kappa() == e2.Kappa() && e1.TwoM() == e2.TwoM())
        {
            unsigned int i1 = orbitals->state_index.at(e1);
            unsigned int i2 = orbitals->state_index.at(e2);

            IntegralOrdering(i1, i2);
            return integrals.at(i1 * num_orbitals + i2);
        }
        else
            return 0.0;
    }

    inline pSpinorOperatorType GetSpinorOperator() const { return op; }

protected:
    inline void IntegralOrdering(unsigned int& i1, unsigned int& i2) const
    {
        if(i1 > i2)
            std::swap(i1, i2);
    }

    pSpinorOperatorType op;
    pOrbitalManagerConst orbitals;
    unsigned int num_orbitals;
    std::map<unsigned int, double> integrals;
};

#endif
