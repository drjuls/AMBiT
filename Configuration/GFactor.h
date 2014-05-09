#ifndef G_FACTOR_H
#define G_FACTOR_H

#include "ManyBodyOperator.h"
#include "HartreeFock/OpIntegrator.h"
#include "Basis/OrbitalManager.h"

/** Sz operator is a one body electron operator that calculates
        < a | Sz | b >
    It stores overlap integrals < a | b > where a.L() == b.L() and
    a, b are both in valence states.
 */
class SzOperator
{
public:
    SzOperator(pOPIntegrator integrator, pOrbitalManagerConst pOrbitals):
        orbitals(pOrbitals)
    {
        num_orbitals = orbitals->size();

        const OrbitalMap& valence_orbitals = *orbitals->GetOrbitalMap(OrbitalClassification::valence);
        auto it = valence_orbitals.begin();
        auto end = valence_orbitals.end();

        while(it != end)
        {
            unsigned int index_it = orbitals->state_index.at(it->first);

            auto jt = it;
            while(jt != end)
            {
                if(it->first.L() == jt->first.L())
                {
                    unsigned int key =
                         index_it * num_orbitals + orbitals->state_index.at(jt->first);

                    integrals[key] = integrator->GetInnerProduct(*it->second, *jt->second);
                }
                jt++;
            }
            it++;
        }
    }

    inline double GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2) const
    {
        double val = 0.;

        if((e1.L() == e2.L()) && (e1.TwoM() == e2.TwoM()))
        {
            unsigned int i1 = orbitals->state_index.at(e1);
            unsigned int i2 = orbitals->state_index.at(e2);

            IntegralOrdering(i1, i2);
            double overlap =  integrals.at(i1 * num_orbitals + i2);

            if(e1.Kappa() == e2.Kappa())
            {   val = e1.M()/(2*e1.L() + 1.);

                if(e1.Kappa() > 0)
                    val = -val;
            }
            else
            {   double Lplushalf = double(e1.L()) + 0.5;
                double M = e1.M();
                val = -sqrt((Lplushalf + M)*(Lplushalf - M))/(2.*Lplushalf);
            }
            val = val * overlap;
        }
        
        return val;
    }

protected:
    inline void IntegralOrdering(unsigned int& i1, unsigned int& i2) const
    {
        if(i1 > i2)
            std::swap(i1, i2);
    }

    pOrbitalManagerConst orbitals;
    unsigned int num_orbitals;
    std::map<unsigned int, double> integrals;
};

typedef boost::shared_ptr<const SzOperator> pSzOperatorConst;

/** GFactorCalculator uses the SzOperator to calculate Lande g-factors
    for all levels of a given symmetry in a LevelMap.
 */
class GFactorCalculator
{
public:
    GFactorCalculator(pOPIntegrator integrator, pOrbitalManagerConst orbitals):
        Sz(pSzOperatorConst(new SzOperator(integrator, orbitals)))
    {}

    inline void CalculateGFactors(LevelMap& levels, const Symmetry& sym) const
    {
        unsigned int NumSolutions = levels.size(sym);
        if(NumSolutions == 0 || sym.GetTwoJ() == 0)
            return;

        std::vector<double> total_Sz = Sz.GetMatrixElement(levels.begin(sym), levels.end(sym));

        auto Sz_it = total_Sz.begin();
        auto it = levels.begin(sym);
        while(it != levels.end(sym) && Sz_it != total_Sz.end())
        {
            it->second->SetgFactor(*Sz_it/sym.GetJ() + 1.);
            it++;
            Sz_it++;
        }
    }

protected:
    ManyBodyOperator<pSzOperatorConst> Sz;
};

#endif
