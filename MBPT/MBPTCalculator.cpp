#include "MBPTCalculator.h"
#include "HartreeFock/ConfigurationParser.h"

MBPTCalculator::MBPTCalculator(pOrbitalManagerConst pOrbitals, const std::string& fermi_orbitals):
    orbitals(pOrbitals), valence(pOrbitals->valence), fermi_orbitals(fermi_orbitals)
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

    // Get orbitals to use if specified
    std::vector<int> valence_orbitals = ConfigurationParser::ParseBasisSize(fermi_orbitals);

    for(int kappa = - (int)max_l - 1; kappa <= (int)max_l; kappa++)
    {
        if(kappa != 0)
        {
            int l = (kappa > 0)? kappa: -kappa-1;
            double valence_energy = std::nan("");

            // Check if found in valence_orbitals
            if(l < valence_orbitals.size() && valence_orbitals[l])
            {
                int pqn = valence_orbitals[l];
                pOrbitalConst ds = valence->GetState(OrbitalInfo(pqn, kappa));
                if(ds)
                {   valence_energy = ds->Energy();
                }
                else
                {   *errstream << "MBPTCalculator::SetValenceEnergies: MBPT/EnergyDenomOrbitals "
                               << OrbitalInfo(pqn, kappa).Name() << " not found." << std::endl;
                }
            }

            if(std::isnan(valence_energy))
            {
                // Get leading excited state
                int pqn = 10;
                it_i = particle->begin();
                while(it_i != particle->end())
                {   pOrbitalConst ds = it_i->second;
                    if((ds->Kappa() == kappa) && (ds->PQN() < pqn))
                    {   pqn = ds->PQN();
                        valence_energy = ds->Energy();
                    }
                    it_i++;
                }
            }

            if(std::isnan(valence_energy))
            {
                // valence state not found in particle set, need to search holes
                int pqn = 0;
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
