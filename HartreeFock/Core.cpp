#include "Include.h"
#include "Core.h"
#include "Universal/MathConstant.h"
#include "Universal/CoulombIntegrator.h"
#include "Universal/Interpolator.h"
#include "ConfigurationParser.h"

Core::Core(pLattice lat, const std::string& filling): StateManager(lat)
{
    if(!filling.empty())
        SetOccupancies(ConfigurationParser::ParseFractionalConfiguration(filling));
}

Core::Core(const Core& other):
    StateManager(other), occupancy(other.occupancy)
{}

const Core& Core::operator=(const Core& other)
{
    StateManager::operator=(other);
    occupancy = other.occupancy;
    return *this;
}

const Core& Core::Copy(const Core& other, pLattice new_lattice)
{
    StateManager::Copy(other, new_lattice);
    occupancy = other.occupancy;

    return *this;
}

Core Core::Copy(pLattice new_lattice) const
{
    Core ret(new_lattice);
    ret.Copy(*this, new_lattice);
    return ret;
}

void Core::Clear()
{
    StateManager::Clear();
    occupancy.clear();
}

double Core::GetOccupancy(OrbitalInfo info) const
{
    OccupationMap::const_iterator it = occupancy.find(info);
    if(it == occupancy.end())
        return 0.0;
    else
        return it->second;
}

double Core::GetOccupancy(pOrbital info) const
{
    return GetOccupancy(OrbitalInfo(info));
}

void Core::SetOccupancies(const OccupationMap& occupancies)
{
    occupancy = occupancies;

    // Remove orbitals with zero occupancy
    OrbitalMap::iterator it = AllStates.begin();
    while(it != AllStates.end())
    {
        if(!occupancy.count(it->first))
        {   OrbitalMap::iterator remove = it;
            it++;
            AllStates.erase(remove);
        }
        else
            it++;
    }

    // Add orbitals that now have occupancies
    OccupationMap::const_iterator occ = occupancy.begin();
    while(occ != occupancy.end())
    {
        if(!AllStates.count(occ->first))
        {   pOrbital s(new Orbital(occ->first.Kappa(), occ->first.PQN()));
            AddState(s);
        }
        occ++;
    }
}

const OccupationMap& Core::GetOccupancies() const
{
    return occupancy;
}

int Core::NumElectrons() const
{
    double ret = 0.0;
    for(OccupationMap::const_iterator it = occupancy.begin(); it != occupancy.end(); it++)
        ret += it->second;

    // Return int: add small value for rounding
    return int(ret + 0.01);
}

void Core::Write(FILE* fp) const
{
    // TODO: Output occupancies

    // Output core
    StateManager::Write(fp);
//
//    // Output open shell states
//    unsigned int num_open = OpenShellStorage.size();
//    fwrite(&num_open, sizeof(unsigned int), 1, fp);
//
//    OrbitalMap::const_iterator it = OpenShellStorage.begin();
//    while(it != OpenShellStorage.end())
//    {
//        it->second->Write(fp);
//        it++;
//    }
//
//    // Output original occupancies
//    num_open = OpenShellStates.size();
//    fwrite(&num_open, sizeof(unsigned int), 1, fp);
//
//    std::map<OrbitalInfo, double>::const_iterator info_it = OpenShellStates.begin();
//    while(info_it != OpenShellStates.end())
//    {
//        unsigned int pqn = info_it->first.PQN();
//        int kappa = info_it->first.Kappa();
//        double occ = info_it->second;
//
//        fwrite(&pqn, sizeof(unsigned int), 1, fp);
//        fwrite(&kappa, sizeof(int), 1, fp);
//        fwrite(&occ, sizeof(double), 1, fp);
//
//        info_it++;
//    }
}

void Core::Read(FILE* fp)
{
    Clear();    
    unsigned int num_core, i;
    unsigned int max_size = 0;

    // TODO: Read occupancies
    // Read original occupancies
//    fread(&num_open, sizeof(unsigned int), 1, fp);
//    for(i = 0; i<num_open; i++)
//    {
//        unsigned int pqn;
//        int kappa;
//        double occ;
//
//        fread(&pqn, sizeof(unsigned int), 1, fp);
//        fread(&kappa, sizeof(int), 1, fp);
//        fread(&occ, sizeof(double), 1, fp);
//
//        OrbitalInfo info(pqn, kappa);
//        OpenShellStates.insert(std::pair<OrbitalInfo, double>(info, occ));
//    }

    // Read core
    fread(&num_core, sizeof(unsigned int), 1, fp);
    for(i = 0; i<num_core; i++)
    {
        pOrbital ds(new Orbital(-1));
        ds->Read(fp);
        max_size = mmax(max_size, ds->Size());
        AddState(ds);
    }

    // Ensure lattice is large enough for all states
    lattice->R(max_size);
}
