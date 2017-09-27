#ifndef LEVEL_MAP_H
#define LEVEL_MAP_H

#include "Level.h"
#include "Include.h"
#include "NonRelConfiguration.h"
#include <vector>
#include <map>
#include <boost/filesystem.hpp>

/** LevelVector consists of Hamiltonian that levels come from,
    common relativistic config list, and vector of levels.
 */
class LevelVector
{
public:
    LevelVector(pHamiltonianID hID = nullptr): hID(hID) {}
    LevelVector(pHamiltonianID hID, pRelativisticConfigList configs, pLevel level):
        hID(hID), configs(configs), levels(1, level) {}
    LevelVector(pRelativisticConfigList configs, pLevel level): configs(configs), levels(1, level)
    {
        if(level)
            hID = level->GetHamiltonianID();
    }

    pHamiltonianID hID {nullptr};
    pRelativisticConfigList configs {nullptr};
    std::vector<pLevel> levels;

    /** Print levels to outstream. min_percentages outside of the range (0, 100] will switch off the printing of configurations. */
    void Print(double min_percentage = 1.0) const;
    void Print(double min_percentage, double max_energy) const;

    /** Print one line per level: J, P, index, E, g, Largest configuration, HamiltonianID. */
    void PrintInline(bool print_leading_configuration = true, bool print_gfactors = true, std::string separator = " ") const;
    void PrintInline(double max_energy, bool print_leading_configuration = true, bool print_gfactors = true, std::string separator = " ") const;

protected:
    /** Print LevelVector to outstream, with all possible options for printing.
        All other print functions call this one.
     */
    void Print(double min_percentage, bool use_max_energy, double max_energy) const;
    void PrintInline(bool use_max_energy, double max_energy, bool print_leading_configuration, bool print_gfactors, std::string separator = " ") const;
};

template <class Value>
struct DereferenceComparator
{
    bool operator()(const Value& first, const Value& second) const { return (*first < *second); }
};

/** LevelStore consists of
    "KeyLibrary", list of pHamiltonianID (which can be a subclass) that are required/desired,
    and storage for associated LevelVectors.
    Iterators act over key library.
 */
class LevelStore
{
protected:
    struct symmetry_match_predicate
    {
        Symmetry m_sym;
        symmetry_match_predicate(const Symmetry& sym): m_sym(sym) {}

        bool operator()(const pHamiltonianID& val) { return m_sym == val->GetSymmetry(); }
    };

public:
    LevelStore() = default;

    /** Set of all keys, whether associated with LevelVector or not. */
    std::set<pHamiltonianID, DereferenceComparator<pHamiltonianID>> keys;

    /** LevelStore has begin() and end() functions over keys. */
    auto begin() -> decltype(keys)::const_iterator
    {   return keys.begin();
    }

    auto end() -> decltype(keys)::const_iterator
    {   return keys.end();
    }

    /** Symmetry filter iterator over key library: matches all keys with correct symmetry. */
    using const_symmetry_iterator = boost::filter_iterator<symmetry_match_predicate, decltype(keys)::const_iterator>;

    const_symmetry_iterator begin(const Symmetry& sym) const
    {   return const_symmetry_iterator(symmetry_match_predicate(sym), keys.begin(), keys.end());
    }

    const_symmetry_iterator end(const Symmetry& sym) const
    {   return const_symmetry_iterator(symmetry_match_predicate(sym), keys.end(), keys.end());
    }

    /** Get LevelVector corresponding to key. */
    virtual LevelVector GetLevels(pHamiltonianID key) = 0;

    /** Store LevelVector (and generally write to file, although this is implementation dependent). */
    virtual void Store(pHamiltonianID key, const LevelVector& level_vector) = 0;
};

typedef std::shared_ptr<LevelStore> pLevelStore;

/** Implementation of LevelStore for when all levels can be stored simultaneously in memory. */
class LevelMap : public LevelStore
{
    typedef std::map<pHamiltonianID, LevelVector, DereferenceComparator<pHamiltonianID>> MapType;

public:
    /** Construct LevelMap with no read/write capacity. */
    LevelMap(pAngularDataLibrary lib);

    /** Construct LevelMap with no read. */
    LevelMap(const std::string& file_id, pAngularDataLibrary lib);

    /** Constructor attempts to read existing level map. */
    LevelMap(pHamiltonianIDConst hamiltonian_example, const std::string& file_id, pAngularDataLibrary lib);

    /** Get LevelVector corresponding to key. */
    virtual LevelVector GetLevels(pHamiltonianID key) override;
    
    /** Store LevelVector and write to file. */
    virtual void Store(pHamiltonianID key, const LevelVector& level_vector) override;

protected:
    /** Read entire LevelMap in single file. */
    void ReadLevelMap(pHamiltonianIDConst hamiltonian_example);

protected:
    std::string filename_prefix;
    pAngularDataLibrary angular_library;

    MapType m_map;
};

/** Implementation of LevelStore that writes everything to files and stores almost nothing.
    Files are separated by Key, in order to speed up retrieval.
 */
class FileSystemLevelStore : public LevelStore
{
public:
    FileSystemLevelStore(const std::string& file_prefix, pAngularDataLibrary lib);
    FileSystemLevelStore(const std::string& dir_name, const std::string& file_prefix, pAngularDataLibrary lib);

    /** Get LevelVector corresponding to key. */
    virtual LevelVector GetLevels(pHamiltonianID key) override;
    
    /** Store LevelVector and write to file. */
    virtual void Store(pHamiltonianID key, const LevelVector& level_vector) override;

protected:
    boost::filesystem::path directory;
    std::string filename_prefix;
    pAngularDataLibrary angular_library;
};

#endif
