#include "Include.h"
#include "OrbitalMap.h"
#include "Universal/Interpolator.h"

OrbitalMap* OrbitalMap::Clone(pLattice new_lattice) const
{
    OrbitalMap* ret;

    if(new_lattice && new_lattice != lattice)
    {   // Interpolate on to new lattice
        ret = new OrbitalMap(new_lattice);

        Interpolator interp(lattice);
        unsigned int order = 6;

        const double* R_old = lattice->R();
        const double* R = new_lattice->R();

        for(auto pair: *this)
        {
            pOrbitalConst s_old = pair.second;

            // Copy kappa, pqn, etc.
            pOrbital s = s_old->Clone();

            double real_orbital_size = R_old[s_old->size() - 1];
            unsigned int new_size = lattice->real_to_lattice(real_orbital_size);
            s->resize(new_size);

            for(unsigned int i = 0; i < new_size; i++)
            {   interp.Interpolate(s_old->f, R[i], s->f[i], s->dfdr[i], order);
                interp.Interpolate(s_old->g, R[i], s->g[i], s->dgdr[i], order);
            }
            
            ret->AddState(s);
        }
    }
    else
    {   //Simply copy orbitals
        ret = new OrbitalMap(lattice);

        auto it = begin();
        while(it != end())
        {
            pOrbital s = it->second->Clone();
            ret->AddState(s);

            it++;
        }
    }

    return ret;
}

pOrbitalConst OrbitalMap::GetState(const OrbitalInfo& info) const
{
    const_iterator it;

    if((it = m_orbitals.find(info)) != m_orbitals.end())
        return it->second;
    else
        return pOrbitalConst();
}

pOrbital OrbitalMap::GetState(const OrbitalInfo& info)
{
    iterator it;

    if((it = m_orbitals.find(info)) != m_orbitals.end())
        return it->second;
    else
        return pOrbital();
}

void OrbitalMap::AddState(pOrbital s)
{
    OrbitalInfo info(s);
    m_orbitals[info] = s;
}

void OrbitalMap::Write(FILE* fp) const
{
    unsigned int num_states = size();
    fwrite(&num_states, sizeof(unsigned int), 1, fp);

    for(auto it = begin(); it != end(); it++)
        it->second->Write(fp);
}

void OrbitalMap::Read(FILE* fp)
{
    clear();
    unsigned int num_core, i;
    unsigned int max_size = 0;

    // Read states
    fread(&num_core, sizeof(unsigned int), 1, fp);
    for(i = 0; i<num_core; i++)
    {
        pOrbital ds(new Orbital(-1));
        ds->Read(fp);
        AddState(ds);
        max_size = mmax(max_size, ds->size());
    }

    // Ensure lattice is large enough for all states
    if(max_size > lattice->size())
        lattice->resize(max_size);

}

void OrbitalMap::AddStates(const OrbitalMap& other)
{
    m_orbitals.insert(other.begin(), other.end());
}

unsigned int OrbitalMap::LargestOrbitalSize() const
{
    unsigned int biggest = 0;
    for(auto it: m_orbitals)
        biggest = mmax(biggest, it.second->size());

    return biggest;
}

void OrbitalMap::Print() const
{
    for(auto it: m_orbitals)
    {
        *outstream << std::setw(6) << it.first.Name()
                   << "  E = " << std::setprecision(12) << it.second->Energy()
                   << "  size: (" << it.second->size()
                   << ") " << std::setprecision(4) << lattice->R(it.second->size()) << std::endl;
    }
}
