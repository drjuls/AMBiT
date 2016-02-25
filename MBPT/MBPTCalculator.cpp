#include "MBPTCalculator.h"

MBPTCalculator::MBPTCalculator(pOrbitalManagerConst pOrbitals):
    orbitals(pOrbitals), valence(pOrbitals->valence)
{
    SetValenceEnergies();
}

MBPTCalculator::~MBPTCalculator(void)
{}

void MBPTCalculator::SetValenceEnergies()
{
    pOrbitalMapConst valence = orbitals->GetOrbitalMap(OrbitalClassification::valence);
    pOrbitalMapConst particle = orbitals->GetOrbitalMap(OrbitalClassification::particle);
    pOrbitalMapConst hole = orbitals->GetOrbitalMap(OrbitalClassification::hole);
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
            double valence_energy = std::nan("");
            int pqn = 10;

            // Get leading excited state (for energy denominator)
            it_i = particle->begin();
            while(it_i != particle->end())
            {   pOrbitalConst ds = it_i->second;
                if((ds->Kappa() == kappa) && (ds->PQN() < pqn))
                {   pqn = ds->PQN();
                    valence_energy = ds->Energy();
                }
                it_i++;
            }

            if(std::isnan(valence_energy))
            {   // valence state not found in particle set, need to search holes
                pqn = 0;
                it_i = hole->begin();
                while(it_i != hole->end())
                {   pOrbitalConst ds = it_i->second;
                    if((ds->Kappa() == kappa) && (ds->PQN() > pqn))
                    {   pqn = ds->PQN();
                        valence_energy = ds->Energy();
                    }
                    it_i++;
                }
            }

            ValenceEnergies.insert(std::pair<int, double>(kappa, valence_energy));
        }
    }
}
