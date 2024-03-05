#ifndef ORBITAL_MANAGER_H
#define ORBITAL_MANAGER_H

#include "HartreeFock/OrbitalMap.h"
#include "HartreeFock/Core.h"
#include "Configuration/IndexedIterator.h"
#include "Include.h"

namespace Ambit
{
/** Classification of orbitals matches maps OrbitalManager. */
enum class OrbitalClassification { deep, hole, particle, high, core, excited, valence, all };

typedef std::map<OrbitalInfo, unsigned int> OrbitalIndex;
typedef std::map<unsigned int, OrbitalInfo> ReverseOrbitalIndex;

/** An overgrown struct holding all OrbitalMaps including core and excited.
    Classification (each Orbital appears in exactly one of the following):
        deep        -- always occupied in the CI model space
        hole        -- sometimes occupied (below Fermi level)
        particle    -- sometimes occupied (above Fermi level)
        high        -- never occupied in the CI model space

    Furthermore each orbital will appear in some of the following:
        core        -- below Fermi level (includes deep + hole)
        excited     -- above Fermi level (includes particle + high)
        valence     -- sometimes occupied in CI model space (includes hole + particle)
        all         -- all orbitals, potentially including some never used in the calculation

    For convenience, all OrbitalMaps are public members, but they will usually be const.
 */
class OrbitalManager
{
public:
    pOrbitalMap deep;       //!< always occupied in the CI model space
    pOrbitalMap hole;       //!< sometimes occupied in the CI model space (below Fermi level)
    pOrbitalMap particle;   //!< sometimes occupied in the CI model space (above Fermi level)
    pOrbitalMap high;       //!< never occupied in the CI model space

    pOrbitalMap core;       //!< below Fermi level == deep + hole
    pOrbitalMap excited;    //!< above Fermi level == particle + high
    pOrbitalMap valence;    //!< sometimes occupied in CI model space == hole + particle

    pOrbitalMap all;        //!< all orbitals >= deep + hole + particle + high == core + excited

    OrbitalIndex state_index;            //!< Unique key for all orbitals
    ReverseOrbitalIndex reverse_state_index;    //!< Unique key for all orbitals

public:
    /** Create a new (empty) pOrbitalMap pointed to by all states. */
    OrbitalManager(pLattice lattice);

    /** Initialise simple frozen core, all-valence model where there are no holes or pure excited states. */
    OrbitalManager(pCore core, pOrbitalMap valence);

    /** Initialise using Read(). */
    OrbitalManager(const std::string& filename);

    pLattice GetLattice() const { return lattice; }

    /** Get total number of orbitals stored (for making keys). */
    unsigned int size() const { return all->size(); }

    pOrbitalMap GetOrbitalMap(OrbitalClassification type);
    pOrbitalMapConst GetOrbitalMap(OrbitalClassification type) const;

    /** Initialise state_index and reverse_state_index. */
    void MakeStateIndexes();

    void Read(const std::string& filename);
    void Write(const std::string& filename) const;

protected:
    pLattice lattice;

    void ReadInfo(FILE* fp, pOrbitalMap& orbitals);
    void WriteInfo(FILE* fp, const pOrbitalMap& orbitals) const;
};

typedef std::shared_ptr<OrbitalManager> pOrbitalManager;
typedef std::shared_ptr<const OrbitalManager> pOrbitalManagerConst;

/** Free function to read state indexes. */
inline void ReadOrbitalIndexes(OrbitalIndex& state_index, FILE* fp)
{
    unsigned int size;
    file_err_handler->fread(&size, sizeof(unsigned int), 1, fp);

    for(unsigned int i = 0; i < size; i++)
    {
        unsigned int index;
        int pqn, kappa;

        file_err_handler->fread(&index, sizeof(unsigned int), 1, fp);
        file_err_handler->fread(&pqn, sizeof(int), 1, fp);
        file_err_handler->fread(&kappa, sizeof(int), 1, fp);

        state_index[OrbitalInfo(pqn, kappa)] = index;
    }
}

/** Free function to write state indexes. */
inline void WriteOrbitalIndexes(const OrbitalIndex& state_index, FILE* fp)
{
    unsigned int size = state_index.size();
    file_err_handler->fwrite(&size, sizeof(unsigned int), 1, fp);

    for(auto& pair: state_index)
    {
        const int pqn = pair.first.PQN();
        const int kappa = pair.first.Kappa();

        file_err_handler->fwrite(&pair.second, sizeof(unsigned int), 1, fp);
        file_err_handler->fwrite(&pqn, sizeof(int), 1, fp);
        file_err_handler->fwrite(&kappa, sizeof(int), 1, fp);
    }
}

/** Reverse state index. */
inline ReverseOrbitalIndex GetReverseIndex(const OrbitalIndex& state_index)
{
    ReverseOrbitalIndex ret;
    for(auto& pair: state_index)
        ret.insert(std::make_pair(pair.second, pair.first));

    return ret;
}

}
#endif
