#include "Include.h"
#include "StateManager.h"
#include "StateIterator.h"
#include "Universal/Interpolator.h"

StateManager::StateManager(Lattice* lat):
    lattice(lat)
{}

StateManager::StateManager(const StateManager& other, Lattice* new_lattice)
{
    if(new_lattice)
        lattice = new_lattice;
    else
        lattice = other.lattice;

    bool interpolate = !(*lattice == *other.lattice);
    Interpolator interp(other.lattice);
    unsigned int order = 6;

    const double* R_old = other.lattice->R();
    const double* R = lattice->R();

    ConstStateIterator it = other.GetConstStateIterator();
    it.First();
    while(!it.AtEnd())
    {
        pOrbitalConst ds_old = it.GetState();

        // Copy kappa, pqn, etc.
        pOrbital ds(new Orbital(*ds_old));

        if(interpolate)
        {
            unsigned int new_size = lattice->real_to_lattice(R_old[ds_old->Size() - 1]);
            double dfdr, dgdr;

            ds->ReSize(new_size);
            for(unsigned int i = 0; i < new_size; i++)
            {
                interp.Interpolate(ds_old->f, R[i], ds->f[i], dfdr, order);
                interp.Interpolate(ds_old->g, R[i], ds->g[i], dgdr, order);
                ds->dfdr[i] = dfdr;
                ds->dgdr[i] = dgdr;
            }
        }

        AddState(ds);
        it.Next();
    }
}

StateManager::~StateManager(void)
{
    Clear();
}

const StateManager& StateManager::operator=(const StateManager& other)
{
    lattice = other.lattice;
    
    ConstStateIterator it = other.GetConstStateIterator();
    it.First();
    while(!it.AtEnd())
    {
        pOrbitalConst old_orbital = it.GetState();
        pOrbital s = GetState(OrbitalInfo(old_orbital));

        if(s == NULL)
        {   s = pOrbital(new Orbital(*it.GetState()));
            AddState(s);
        }
        else
            *s = *old_orbital;

        it.Next();
    }

    return *this;
}

pOrbitalConst StateManager::GetState(const OrbitalInfo& info) const
{
    StateSet::const_iterator it;

    if((it = AllStates.find(info)) != AllStates.end())
        return it->second;
    else
        return pOrbitalConst();
}

pOrbital StateManager::GetState(const OrbitalInfo& info)
{
    StateSet::iterator it;

    if((it = AllStates.find(info)) != AllStates.end())
        return it->second;
    else
        return pOrbital();
}

void StateManager::AddState(pOrbital s)
{
    OrbitalInfo info(s);
    AllStates.insert(StateSet::value_type(info, s));
}

void StateManager::Clear()
{
    AllStates.clear();
}

double StateManager::TestOrthogonality() const
{
    double max_orth = 0.;

    ConstStateIterator it = GetConstStateIterator();
    ConstStateIterator jt = GetConstStateIterator();

    it.First();
    while(!it.AtEnd())
    {
        jt = it;
        jt.Next();
        while(!jt.AtEnd())
        {
            double orth = fabs(it.GetState()->Overlap(*jt.GetState(), lattice));
            if(orth > max_orth)
                max_orth = orth;

            jt.Next();
        }
        it.Next();
    }

    return max_orth;
}

StateIterator StateManager::GetStateIterator()
{
    StateIterator it(this);
    return it;
}

ConstStateIterator StateManager::GetConstStateIterator() const
{
    ConstStateIterator it(this);
    return it;
}

void StateManager::Write(FILE* fp) const
{
    unsigned int num_states = NumStates();

    fwrite(&num_states, sizeof(unsigned int), 1, fp);
    ConstStateIterator it = GetConstStateIterator();
    while(!it.AtEnd())
    {
        it.GetState()->Write(fp);
        it.Next();
    }
}

void StateManager::Read(FILE* fp)
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
            max_size = mmax(max_size, ds.Size());
        }
    }

    // Ensure lattice is large enough for new states
    lattice->R(max_size);
}
