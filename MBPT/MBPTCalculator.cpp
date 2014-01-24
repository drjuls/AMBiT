#include "MBPTCalculator.h"
#include "Universal/MathConstant.h"
#include "Universal/CoulombIntegrator.h"
#include "HartreeFock/StateIntegrator.h"

MBPTCalculator::MBPTCalculator(pLattice lat, const Core* atom_core, const ExcitedStates* excited_states):
    lattice(lat), core(atom_core), excited(excited_states), delta(0.)
{
    SetValenceEnergies();
}

MBPTCalculator::~MBPTCalculator(void)
{}

void MBPTCalculator::SetValenceEnergies()
{
    ConstStateIterator it_i = excited->GetConstStateIterator();
    ValenceEnergies.clear();

    // Get maximum angular momentum in excited states
    unsigned int max_l = 0;
    it_i.First();
    while(!it_i.AtEnd())
    {   max_l = mmax(it_i.GetState()->L(), max_l);
        it_i.Next();
    }

    for(int kappa = - (int)max_l - 1; kappa <= (int)max_l; kappa++)
        if(kappa != 0)
        {
            double valence_energy = 0.;
            unsigned int pqn = 10;

            // Get leading state (for energy denominator)
            it_i.First();
            while(!it_i.AtEnd())
            {   pOrbitalConst ds = it_i.GetState();
                if((ds->Kappa() == kappa) && (ds->GetPQN() < pqn))
                {   pqn = ds->GetPQN();
                    valence_energy = ds->GetEnergy();
                }
                it_i.Next();
            }

            ValenceEnergies.insert(std::pair<int, double>(kappa, valence_energy));
        }
}
