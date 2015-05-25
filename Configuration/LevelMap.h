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
    bool operator() (const pHamiltonianIDConst& lhs, const pHamiltonianIDConst& rhs) const
    {   return *lhs < *rhs; }
};

/** LevelStore consists of
    "KeyLibrary", list of pHamiltonianID (which can be a subclass) that are required/desired,
    and storage for associated LevelVectors.
    Iterators act over key library.
 */
class LevelStore
{
protected:
    class symmetry_match_predicate;

public:
    typedef std::vector<pHamiltonianID> KeyLibrary;

    LevelStore() {}

    // Iterators over key library
    typedef KeyLibrary::const_iterator const_iterator;

    const_iterator begin() const { return m_lib.begin(); }
    unsigned int count(const pHamiltonianIDConst& key) const; //!< Test for existence of key in key library.
    bool empty() const { return m_lib.empty(); }
    const_iterator end() const { return m_lib.end(); }
    const_iterator find(const pHamiltonianIDConst& key) const;
    unsigned int size() const { return m_lib.size(); }


    /** Symmetry filter iterator over key library: matches all keys with correct symmetry. */
    typedef boost::filter_iterator<symmetry_match_predicate, const_iterator> const_symmetry_iterator;
    const_symmetry_iterator begin(const Symmetry& sym) const
    {   return const_symmetry_iterator(symmetry_match_predicate(sym), begin(), end());
    }
    const_symmetry_iterator end(const Symmetry& sym) const
    {   return const_symmetry_iterator(symmetry_match_predicate(sym), end(), end());
    }

    /** Insert into key library.
        Copy RelativisticConfigList from key into library if key already exists.
        Return iterator to key in library.
     */
    virtual const_iterator insert(pHamiltonianID key);

    /** Get LevelVector corresponding to key. */
    virtual LevelVector GetLevels(pHamiltonianID key) = 0;

    /** Store LevelVector (and generally write to file, although this is implementation dependent). */
    virtual void Store(pHamiltonianID key, const LevelVector& level_vector) = 0;

protected:
    KeyLibrary m_lib;

protected:
    class symmetry_match_predicate
    {
    public:
        symmetry_match_predicate(const Symmetry& sym): m_sym(sym) {}
        symmetry_match_predicate(const symmetry_match_predicate& other): m_sym(other.m_sym) {}
        ~symmetry_match_predicate() {}
        
        symmetry_match_predicate operator=(const symmetry_match_predicate& other)
        {   m_sym = other.m_sym; return *this;
        }
        
        bool operator()(const pHamiltonianID& val) { return m_sym == val->GetSymmetry(); }
    protected:
        Symmetry m_sym;
    };
};

typedef std::shared_ptr<LevelStore> pLevelStore;

/** Implementation of LevelStore for when all levels can be stored simultaneously in memory. */
class LevelMap : public LevelStore
{
    typedef std::map<pHamiltonianID, LevelVector, pHamiltonianIDComparator> MapType;

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
    Files are separated by Symmetry, in order to speed up retrieval.
 */
class FileSystemLevelStore : public LevelStore
{
public:
    FileSystemLevelStore(const std::string& file_prefix, pAngularDataLibrary lib);

    /** Get LevelVector corresponding to key. */
    virtual LevelVector GetLevels(pHamiltonianID key) override;
    
    /** Store LevelVector and write to file. */
    virtual void Store(pHamiltonianID key, const LevelVector& level_vector) override;

protected:
    std::string filename_prefix;
    pAngularDataLibrary angular_library;
};

/** Print LevelVector to outstream, with all possible options for printing.
    All other print functions call this one.
 */
void Print(const LevelVector& levels, double min_percentage, bool use_max_energy, double max_energy);

#endif
