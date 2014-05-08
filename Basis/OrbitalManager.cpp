#include "OrbitalManager.h"
#include "Include.h"
#include "IndexedIterator.h"

void OrbitalManager::UpdateStateIndexes()
{
    state_index.clear();
    reverse_state_index.clear();

    // Iterate through states, assign in order
    for(indexed_iterator<OrbitalMap::const_iterator> it(all->begin()); it.base() != all->end(); it++)
    {
        state_index.insert(std::make_pair(it->first, it.index()));
        reverse_state_index.insert(std::make_pair(it.index(), it->first));
    }
}

void OrbitalManager::Read(const std::string& filename)
{
    FILE* fp = fopen(filename.c_str(), "rb");

    all->Read(fp);

    // Read OrbitalInfo for all other maps and get corresponding Orbital from all
    ReadInfo(fp, deep);
    ReadInfo(fp, hole);
    ReadInfo(fp, particle);
    ReadInfo(fp, high);
    ReadInfo(fp, core);
    ReadInfo(fp, excited);
    ReadInfo(fp, valence);
}

void OrbitalManager::Write(const std::string& filename) const
{
    FILE* fp = fopen(filename.c_str(), "wb");

    all->Write(fp);

    // Write OrbitalInfo for all other maps.
    WriteInfo(fp, deep);
    WriteInfo(fp, hole);
    WriteInfo(fp, particle);
    WriteInfo(fp, high);
    WriteInfo(fp, core);
    WriteInfo(fp, excited);
    WriteInfo(fp, valence);
}

void OrbitalManager::ReadInfo(FILE* fp, pOrbitalMap& orbitals)
{
    unsigned int size;
    int pqn, kappa;

    fread(&size, sizeof(unsigned int), 1, fp);

    for(unsigned int i = 0; i < size; i++)
    {
        fread(&pqn, sizeof(int), 1, fp);
        fread(&kappa, sizeof(int), 1, fp);

        pOrbital orb = all->GetState(OrbitalInfo(pqn, kappa));
        orbitals->AddState(orb);
    }
}

void OrbitalManager::WriteInfo(FILE* fp, const pOrbitalMap& orbitals) const
{
    unsigned int size = orbitals->size();
    fwrite(&size, sizeof(unsigned int), 1, fp);

    for(auto pair: *orbitals)
    {
        int pqn = pair.first.PQN();
        int kappa = pair.first.Kappa();

        fwrite(&pqn, sizeof(int), 1, fp);
        fwrite(&kappa, sizeof(int), 1, fp);
    }
}
