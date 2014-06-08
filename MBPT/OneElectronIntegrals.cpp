#include "OneElectronIntegrals.h"

unsigned int OneElectronIntegrals::CalculateOneElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, bool check_size_only)
{
    unsigned int i1, i2;
    pOrbitalConst s1, s2;

    std::set<unsigned int> found_keys;   // For check_size_only

    auto it_1 = orbital_map_1->begin();
    while(it_1 != orbital_map_1->end())
    {
        i1 = orbitals->state_index.at(it_1->first);
        s1 = it_1->second;

        auto it_2 = orbital_map_2->begin();
        while(it_2 != orbital_map_2->end())
        {
            if(it_1->first.Kappa() == it_2->first.Kappa())
            {
                i2 = orbitals->state_index.at(it_2->first);
                s2 = it_2->second;

                if(check_size_only)
                    found_keys.insert(GetKey(i1, i2));
                else
                    integrals[GetKey(i1, i2)] = op->GetMatrixElement(*s1, *s2);
            }
            it_2++;
        }
        it_1++;
    }

    if(check_size_only)
        return found_keys.size();
    else
        return integrals.size();
}

void OneElectronIntegrals::Read(const std::string& filename)
{
    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
    {   *errstream << "SlaterIntegrals::Read: file " << filename << " not found." << std::endl;
        exit(1);
    }

    OrbitalIndex old_state_index;
    ReadOrbitalIndexes(old_state_index, fp);
    unsigned int old_num_states = old_state_index.size();

    unsigned int num_integrals;
    unsigned int old_key;
    double value;

    fread(&num_integrals, sizeof(unsigned int), 1, fp);

    for(unsigned int i = 0; i < num_integrals; i++)
    {
        fread(&old_key, sizeof(unsigned int), 1, fp);
        fread(&value, sizeof(double), 1, fp);

        unsigned int i1 = old_key/old_num_states;
        unsigned int i2 = old_key - i1 * old_num_states;
        unsigned int new_key = GetKey(i1, i2);

        auto it = integrals.find(new_key);
        if(it == integrals.end())
            integrals[new_key] = value;
        else
            it->second += value;
    }

    fclose(fp);
}

void OneElectronIntegrals::Write(const std::string& filename) const
{
    FILE* fp = fopen(filename.c_str(), "wb");

    // Write state index
    WriteOrbitalIndexes(orbitals->state_index, fp);

    unsigned int num_integrals = size();
    fwrite(&num_integrals, sizeof(unsigned int), 1, fp);

    for(auto& pair: integrals)
    {
        fwrite(&pair.first, sizeof(unsigned int), 1, fp);
        fwrite(&pair.second, sizeof(double), 1, fp);
    }
    
    fclose(fp);
}
