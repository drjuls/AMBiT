#include "Include.h"
#include "StateManager.h"
#include "StateIterator.h"
#include "Universal/Interpolator.h"

StateManager::StateManager(pLattice lat):
    lattice(lat)
{}

StateManager::StateManager(const StateManager& other)
{
    (*this) = other;
}

StateManager::~StateManager(void)
{
    Clear();
}

const StateManager& StateManager::operator=(const StateManager& other)
{
    lattice = other.lattice;
    AllStates = other.AllStates;

    return *this;
}

const StateManager& StateManager::Copy(const StateManager& other, pLattice new_lattice)
{
    bool interpolate = false;
    if(new_lattice && new_lattice != other.lattice)
    {   lattice = new_lattice;
        interpolate = true;
    }
    else
        lattice = other.lattice;

    AllStates.clear();

    ConstStateIterator it = other.GetConstStateIterator();

    if(interpolate)
    {
        Interpolator interp(other.lattice);
        unsigned int order = 6;

        const double* R_old = other.lattice->R();
        const double* R = lattice->R();

        it.First();
        while(!it.AtEnd())
        {
            pOrbitalConst old_orbital = it.GetState();

            // Copy kappa, pqn, etc.
            pOrbital s(new Orbital(old_orbital->Kappa(), old_orbital->GetPQN(), old_orbital->GetEnergy()));

            double real_orbital_size = R_old[old_orbital->Size() - 1];
            unsigned int new_size = lattice->real_to_lattice(real_orbital_size);
            s->ReSize(new_size);

            for(unsigned int i = 0; i < new_size; i++)
            {   interp.Interpolate(old_orbital->f, R[i], s->f[i], s->dfdr[i], order);
                interp.Interpolate(old_orbital->g, R[i], s->g[i], s->dgdr[i], order);
            }
            
            AddState(s);
            it.Next();
        }

    }
    else
    {   //Simply copy orbitals
        it.First();
        while(!it.AtEnd())
        {
            pOrbitalConst old_orbital = it.GetState();

            pOrbital s(new Orbital(*old_orbital));
            AddState(s);

            it.Next();
        }
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
