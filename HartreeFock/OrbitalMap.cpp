#include "Include.h"
#include "OrbitalMap.h"
#include "Universal/Interpolator.h"

OrbitalMap::OrbitalMap(const OrbitalMap& other)
{
    (*this) = other;
}

OrbitalMap::~OrbitalMap(void)
{
    clear();
}

const OrbitalMap& OrbitalMap::operator=(const OrbitalMap& other)
{
    lattice = other.lattice;
    m_orbitals = other.m_orbitals;

    return *this;
}

const OrbitalMap& OrbitalMap::Copy(const OrbitalMap& other, pLattice new_lattice)
{
    bool interpolate = false;
    if(new_lattice && new_lattice != other.lattice)
    {   lattice = new_lattice;
        interpolate = true;
    }
    else
        lattice = other.lattice;

    clear();

    if(interpolate)
    {
        Interpolator interp(other.lattice);
        unsigned int order = 6;

        const double* R_old = other.lattice->R();
        const double* R = lattice->R();

        for(auto pair: other)
        {
            pOrbitalConst old_orbital = pair.second;

            // Copy kappa, pqn, etc.
            pOrbital s(new Orbital(old_orbital->Kappa(), old_orbital->PQN(), old_orbital->Energy()));

            double real_orbital_size = R_old[old_orbital->size() - 1];
            unsigned int new_size = lattice->real_to_lattice(real_orbital_size);
            s->resize(new_size);

            for(unsigned int i = 0; i < new_size; i++)
            {   interp.Interpolate(old_orbital->f, R[i], s->f[i], s->dfdr[i], order);
                interp.Interpolate(old_orbital->g, R[i], s->g[i], s->dgdr[i], order);
            }
            
            AddState(s);
        }

    }
    else
    {   //Simply copy orbitals
        auto it = other.begin();
        while(it != other.end())
        {
            pOrbitalConst old_orbital = it->second;

            pOrbital s(new Orbital(*old_orbital));
            AddState(s);

            it++;
        }
    }

    return *this;
}

OrbitalMap OrbitalMap::Copy(pLattice new_lattice) const
{
    OrbitalMap ret(new_lattice);
    ret.Copy(*this, new_lattice);
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
        it->second->SpinorFunction::Write(fp);
}

void OrbitalMap::Read(FILE* fp)
{
    unsigned int num_core, i;
    Orbital ds(-1);
    unsigned int max_size = 0;

    // Read states
    fread(&num_core, sizeof(unsigned int), 1, fp);
    for(i = 0; i<num_core; i++)
    {
        ds.Read(fp);
        pOrbital current = GetState(OrbitalInfo(&ds));
        if(current)
        {   *current = ds;
            max_size = mmax(max_size, ds.size());
        }
    }

    // Ensure lattice is large enough for new states
    lattice->R(max_size);
}

void OrbitalMap::AddStates(OrbitalMap& other)
{
    for(auto it = other.begin(); it != other.end(); it++)
        AddState(it->second);
}
