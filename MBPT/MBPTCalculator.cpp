#include "MBPTCalculator.h"

MBPTCalculator::MBPTCalculator(pOrbitalManagerConst pOrbitals):
    orbitals(pOrbitals)
{
    SetValenceEnergies();
}

MBPTCalculator::~MBPTCalculator(void)
{}

void MBPTCalculator::SetValenceEnergies()
{
    pOrbitalMapConst valence = orbitals->GetOrbitalMap(OrbitalClassification::valence);
    ValenceEnergies.clear();

    // Get maximum angular momentum in excited states
    unsigned int max_l = 0;
    auto it_i = valence->begin();
    while(it_i != valence->end())
    {   max_l = mmax(it_i->first.L(), max_l);
        it_i++;
    }

    for(int kappa = - (int)max_l - 1; kappa <= (int)max_l; kappa++)
    {
        if(kappa != 0)
        {
            double valence_energy = 0.;
            unsigned int pqn = 10;

            // Get leading state (for energy denominator)
            it_i = valence->begin();
            while(it_i != valence->end())
            {   pOrbitalConst ds = it_i->second;
                if((ds->Kappa() == kappa) && (ds->PQN() < pqn))
                {   pqn = ds->PQN();
                    valence_energy = ds->Energy();
                }
                it_i++;
            }

            ValenceEnergies.insert(std::pair<int, double>(kappa, valence_energy));
        }
    }
}
