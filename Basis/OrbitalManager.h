#ifndef ORBITAL_MANAGER_H
#define ORBITAL_MANAGER_H

#include "HartreeFock/OrbitalMap.h"
#include "HartreeFock/Core.h"

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
        all         -- all orbitals
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

    pOrbitalMap all;        //!< all orbitals == deep + hole + particle + high == core + excited

    std::map<OrbitalInfo, unsigned int> state_index;
    std::map<unsigned int, OrbitalInfo> reverse_state_index;

public:
    OrbitalManager() {}

    unsigned int size() const { return all->size(); }
    void UpdateStateIndexes();

    void Read(const std::string& filename);
    void Write(const std::string& filename) const;

protected:
    void ReadInfo(FILE* fp, pOrbitalMap& orbitals);
    void WriteInfo(FILE* fp, const pOrbitalMap& orbitals) const;
};

typedef boost::shared_ptr<OrbitalManager> pOrbitalManager;
typedef boost::shared_ptr<const OrbitalManager> pOrbitalManagerConst;

enum class OrbitalClassification { deep, hole, particle, high, core, excited, valence, all };

#endif
