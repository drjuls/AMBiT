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

class symmetry_match_predicate
{
public:
    symmetry_match_predicate(const Symmetry& sym): m_sym(sym) {}
    symmetry_match_predicate(const symmetry_match_predicate& other): m_sym(other.m_sym) {}
    ~symmetry_match_predicate() {}

    symmetry_match_predicate operator=(const symmetry_match_predicate& other)
    {   m_sym = other.m_sym; return *this;
    }

    bool operator()(const std::pair<pHamiltonianID, LevelVector>& val) { return m_sym == val.first->GetSymmetry(); }
protected:
    Symmetry m_sym;
};

/** Symmetry filter iterator over LevelMap */
typedef boost::filter_iterator<symmetry_match_predicate, LevelMap::iterator> LevelMap_symmetry_iterator;
typedef boost::filter_iterator<symmetry_match_predicate, LevelMap::const_iterator> const_LevelMap_symmetry_iterator;

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
pLevelMap ReadLevelMap(pHamiltonianIDConst hamiltonian_example, const std::string& filename, pAngularDataLibrary angular_library);

/** Append object pointed to by LevelMap::iterator to file.
    Calls WriteLevelMap if file doesn't exist, otherwise appends to end and updates number of records stored.
 */
void AppendLevelMap(const LevelMap::const_iterator it, const std::string& filename);

/** Print LevelVector to outstream, with all possible options for printing.
    All other print functions call this one.
 */
void Print(const LevelVector& levels, double min_percentage, bool use_max_energy, double max_energy);

#endif
