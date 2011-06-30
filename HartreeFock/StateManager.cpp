#include "Include.h"
#include "StateManager.h"
#include "StateIterator.h"
#include "Universal/Interpolator.h"

StateManager::StateManager(Lattice* lat, unsigned int atomic_number, int ion_charge):
    lattice(lat), Z(atomic_number), Charge(ion_charge)
{}

StateManager::StateManager(const StateManager& other, Lattice* new_lattice):
    Z(other.Z), Charge(other.Charge)
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
    const double* dR = lattice->dR();

    ConstStateIterator it = other.GetConstStateIterator();
    it.First();
    while(!it.AtEnd())
    {
        const Orbital* ds_old = it.GetState();

        // Copy kappa, pqn, etc.
        Orbital* ds = new Orbital(*ds_old);

        if(interpolate)
        {
            unsigned int new_size = lattice->real_to_lattice(R_old[ds_old->Size() - 1]);
            double dfdr, dgdr;

            ds->ReSize(new_size);
            for(unsigned int i = 0; i < new_size; i++)
            {
                interp.Interpolate(ds_old->f, R[i], ds->f[i], dfdr, order);
                interp.Interpolate(ds_old->g, R[i], ds->g[i], dgdr, order);
                ds->df[i] = dfdr * dR[i];
                ds->dg[i] = dgdr * dR[i];
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

const Orbital* StateManager::GetState(const OrbitalInfo& info) const
{
    StateSet::const_iterator it;

    if((it = AllStates.find(info)) != AllStates.end())
        return it->second.GetState();
    else
        return NULL;
}

Orbital* StateManager::GetState(const OrbitalInfo& info)
{
    StateSet::iterator it;

    if((it = AllStates.find(info)) != AllStates.end())
        return it->second.GetState();
    else
        return NULL;
}

void StateManager::AddState(Orbital* s)
{
    OrbitalInfo info(s);
    StatePointer sp(s);

    AllStates.insert(StateSet::value_type(info, sp));
}

void StateManager::Clear()
{
    // Delete current states
    StateSet::iterator it = AllStates.begin();
    while(it != AllStates.end())
    {   it->second.DeleteState();
        it++;
    }
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
    Orbital ds;
    unsigned int max_size = 0;

    // Read states
    fread(&num_core, sizeof(unsigned int), 1, fp);
    for(i = 0; i<num_core; i++)
    {
        ds.Read(fp);
        Orbital* current = GetState(OrbitalInfo(&ds));
        if(current)
        {   *current = ds;
            max_size = mmax(max_size, ds.Size());
        }
    }

    // Ensure lattice is large enough for new states
    lattice->R(max_size);
}
