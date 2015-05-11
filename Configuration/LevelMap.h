#ifndef LEVEL_MAP_H
#define LEVEL_MAP_H

#include "Level.h"
#include "Include.h"
#include "NonRelConfiguration.h"
#include <vector>
#include <map>

/** All levels in LevelVector should come from the same Hamiltonian, although this is not enforced. */
typedef std::vector<pLevel> LevelVector;
void Print(const LevelVector& levels, double min_percentage = 1.0);
void Print(const LevelVector& levels, double min_percentage, double max_energy);

struct pHamiltonianIDComparator {
    bool operator() (const pHamiltonianID& lhs, const pHamiltonianID& rhs) const
    {   return *lhs < *rhs; }
};

/** Map from pHamiltonianID (which can be a subclass) to LevelVector. */
typedef std::map<pHamiltonianID, LevelVector, pHamiltonianIDComparator> LevelMap;
typedef std::shared_ptr<LevelMap> pLevelMap;
typedef std::shared_ptr<const LevelMap> pLevelMapConst;

/** Write LevelMap object to file.
    LevelMap maps HamiltonianID -> LevelVector, where the HamiltonianID defines
        Write(FILE*) and Read(FILE*) which define a unique HamiltonianID, as well as
        GetRelativisticConfigList().
 */
void WriteLevelMap(const LevelMap& level_map, const std::string& filename);

/** Return LevelMap object read from file.
    LevelMap maps pHamiltonianID -> LevelVector, where the HamiltonianID defines
        Write(FILE*) and Read(FILE*) which define a unique HamiltonianID, as well as
        GetRelativisticConfigList().
 */
pLevelMap ReadLevelMap(pHamiltonianIDConst hamiltonian_example, const std::string& filename, const std::string& angular_directory);

/** Print LevelVector to outstream, with all possible options for printing.
    All other print functions call this one.
 */
void Print(const LevelVector& levels, double min_percentage, bool use_max_energy, double max_energy);

#endif
