#include "Include.h"
#include "Core.h"
#include "Universal/MathConstant.h"
#include "Universal/Interpolator.h"
#include "ConfigurationParser.h"

Core::Core(pLattice lat, const std::string& filling): OrbitalMap(lat)
{
    if(!filling.empty())
        SetOccupancies(ConfigurationParser::ParseFractionalConfiguration(filling));
}

Core::Core(const Core& other):
    OrbitalMap(other), occupancy(other.occupancy)
{}

Core::Core(Core&& other):
    OrbitalMap(other), occupancy(other.occupancy)
{}

Core::Core(const OrbitalMap& other):
    OrbitalMap(other)
{
    for(auto orbital: m_orbitals)
        occupancy[orbital.first] = 2 * abs(orbital.first.Kappa());
}

const Core& Core::operator=(const Core& other)
{
    OrbitalMap::operator=(other);
    occupancy = other.occupancy;
    return *this;
}

const Core& Core::Clone(const Core& other, pLattice new_lattice)
{
    OrbitalMap::Clone(other, new_lattice);
    occupancy = other.occupancy;

    return *this;
}

Core Core::Clone(pLattice new_lattice) const
{
    Core ret(new_lattice);
    ret.Clone(*this, new_lattice);
    return ret;
}

void Core::clear()
{
    OrbitalMap::clear();
    occupancy.clear();
}

auto Core::erase(const_iterator position) -> iterator
{
    iterator ret = OrbitalMap::erase(position);
    occupancy.erase(position->first);

    return ret;
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
    auto it = m_orbitals.begin();
    while(it != end())
    {
        if(!occupancy.count(it->first))
            it = m_orbitals.erase(it);
        else
            it++;
    }

    // Add orbitals that now have occupancies
    OccupationMap::const_iterator occ = occupancy.begin();
    while(occ != occupancy.end())
    {
        if(!count(occ->first))
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
    // Output core
    OrbitalMap::Write(fp);

    // Output occupancies
    unsigned int size = occupancy.size();
    fwrite(&size, sizeof(unsigned int), 1, fp);
    for(auto& pair: occupancy)
    {
        int kappa = pair.first.Kappa();
        int pqn = pair.first.PQN();
        double occ = pair.second;

        fwrite(&kappa, sizeof(int), 1, fp);
        fwrite(&pqn, sizeof(int), 1, fp);
        fwrite(&occ, sizeof(double), 1, fp);
    }
}

void Core::Read(FILE* fp)
{
    clear();
    occupancy.clear();
    unsigned int occ_size, i;

    // Read core
    OrbitalMap::Read(fp);

    // Read occupancies
    fread(&occ_size, sizeof(unsigned int), 1, fp);
    for(i = 0; i < occ_size; i++)
    {
        int kappa;
        int pqn;
        double occ;

        fread(&kappa, sizeof(int), 1, fp);
        fread(&pqn, sizeof(int), 1, fp);
        fread(&occ, sizeof(double), 1, fp);

        OrbitalInfo info(pqn, kappa);
        occupancy.insert(std::pair<OrbitalInfo, double>(info, occ));
    }
}
