#ifndef LEVEL_H
#define LEVEL_H

#include "Symmetry.h"
#include "RelativisticConfiguration.h"
#include "Universal/Enums.h"
#include <boost/iterator/filter_iterator.hpp>

/** LevelID provides a key to uniquely identify an eigenstate of the Hamiltonian
    using J, Parity, and an integer index drawn from the energy ordering of the eigenstates.
 */
class LevelID
{
public:
    LevelID(): m_twoJ(-1), m_parity(Parity::even), m_ID(0) {}
    LevelID(int two_J, Parity parity, unsigned int level_id): m_twoJ(two_J), m_parity(parity), m_ID(level_id) {}
    LevelID(const Symmetry& sym, unsigned int level_id): m_twoJ(sym.GetTwoJ()), m_parity(sym.GetParity()), m_ID(level_id) {}
    LevelID(const std::string& name);

    bool operator<(const LevelID& other) const;
    inline bool operator>(const LevelID& other) const
    {   return (other < *this);
    }

    inline bool operator==(const LevelID& other) const
    {   return ((m_twoJ == other.m_twoJ) && (m_parity == other.m_parity) && (m_ID == other.m_ID));
    }
    inline bool operator!=(const LevelID& other) const
    {   return ((m_twoJ != other.m_twoJ) || (m_parity != other.m_parity) || (m_ID != other.m_ID));
    }

    double GetJ() const { return double(m_twoJ)/2.0; }
    int GetTwoJ() const { return m_twoJ; }
    Parity GetParity() const { return m_parity; }
    unsigned int GetID() const { return m_ID; }
    Symmetry GetSymmetry() const { return Symmetry(m_twoJ, m_parity); }

    /** Name is twoJ, followed by parity (e/o), followed by id. */
    std::string Name() const;

protected:
    int m_twoJ;
    Parity m_parity;
    unsigned int m_ID;
};

/** A Level is an eigenstate of the Hamiltonian with eigenvector of length NumCSFs.
    The symmetry of the eigenstate is given by the stored RelativisticConfigList.
 */
class Level
{
    friend class LevelMap;
public:
    /** Initialise Level with energy eigenvalue, eigenvector, and pointer to relativistic configurations (CSFs).
        PRE: length(csf_eigenvector) == configlist->NumCSFs() == numCSFs (if supplied).
     */
    Level(const double& energy, const double* csf_eigenvector, pRelativisticConfigListConst configlist, unsigned int numCSFs = 0);
    Level(const Level& other);
    Level(Level&& other);
    ~Level();

    const Level& operator=(const Level& other);
    Level& operator=(Level&& other);

    double GetEnergy() const { return eigenvalue; }
    void SetEnergy(double energy) { eigenvalue = energy; }

    double GetgFactor() const { return gFactor; }
    void SetgFactor(double g_factor) { gFactor = g_factor; }

    const double* GetEigenvector() const { return eigenvector; }
    unsigned int GetEigenvectorLength() const { return N; }     //!< Eigenvector length = NumCSFs

    pRelativisticConfigListConst GetRelativisticConfigList() const { return configs; }

//    inline ConfigurationSet* GetConfigurationSet()
//    {   return mConfigurationSet;
//    }
//    inline Configuration GetLeadingConfiguration()
//    {   return GetConfigurationSet()->GetLargestConfiguration();
//    }
//    inline TransitionSet* GetTransitionSet()
//    {   return mTransitionSet;
//    }

protected:
    Level(): N(0), eigenvector(nullptr) {}  //!< For use by LevelMap::Read()

    double eigenvalue;
    unsigned int N;     // length of eigenvector = configs->NumCSFs()
    double* eigenvector;
    pRelativisticConfigListConst configs;
    double gFactor;
//    ConfigurationSet* mConfigurationSet;
//    TransitionSet* mTransitionSet;
};

typedef std::shared_ptr<Level> pLevel;
typedef std::shared_ptr<const Level> pLevelConst;

/** Map LevelID to pLevel. For each level, M = J.
    As well as iterators over all levels, LevelMap provides a filter_iterator to select only
    levels with a given Symmetry (J, Parity).
 */
class LevelMap : private std::map<LevelID, pLevel>
{
private:
    typedef std::map<LevelID, pLevel> Parent;

protected:
    /** Predicate for use in boost::filter_iterator to match iterators based on symmetry. */
    class symmetry_match_predicate
    {
    public:
        symmetry_match_predicate(const Symmetry& sym): m_sym(sym) {}
        symmetry_match_predicate(const symmetry_match_predicate& other): m_sym(other.m_sym) {}
        ~symmetry_match_predicate() {}

        symmetry_match_predicate operator=(const symmetry_match_predicate& other)
        {   m_sym = other.m_sym; return *this;
        }

        bool operator()(Parent::value_type x) { return m_sym == x.first.GetSymmetry(); }
    protected:
        Symmetry m_sym;
    };

public:
    LevelMap(): Parent() {}

    typedef Parent::iterator iterator;
    typedef Parent::const_iterator const_iterator;

    using Parent::operator=;
    using Parent::operator[];
    using Parent::begin;
    using Parent::clear;
    using Parent::empty;
    using Parent::end;
    using Parent::find;
    using Parent::insert;
    using Parent::size;

    /** symmetry_iterators are boost::filter_iterators that include only elements of the map
        that have the correct symmetry.
     */
    typedef boost::filter_iterator<symmetry_match_predicate, iterator> symmetry_iterator;
    typedef boost::filter_iterator<symmetry_match_predicate, const_iterator> const_symmetry_iterator;

    /** Return a symmetry_iterator that selects on symmetry sym. */
    symmetry_iterator begin(const Symmetry& sym)
    {   return symmetry_iterator(symmetry_match_predicate(sym), begin(), end());
    }

    const_symmetry_iterator begin(const Symmetry& sym) const
    {   return const_symmetry_iterator(symmetry_match_predicate(sym), begin(), end());
    }

    const_symmetry_iterator end(const Symmetry& sym) const
    {   return const_symmetry_iterator(symmetry_match_predicate(sym), end(), end());
    }

    /** Number of stored levels with symmetry sym. */
    unsigned int size(const Symmetry& sym) const;

    /** Return set of all symmetries for which size(symmetry) > 0. */
    std::set<Symmetry> GetSymmetries() const;

    /** Print all levels with given symmetry.
        Include parentage of all configurations which contribute more than min_percentage.
        PRE: Assumes all levels with same symmetry point to same RelativisticConfigList.
     */
    inline void Print(const Symmetry& sym, double min_percentage = 1.0) const
    {   Print(sym, min_percentage, false, 0.0);
    }

    /** Print all levels with given symmetry and energy < max_energy.
        Include parentage of all configurations which contribute more than min_percentage.
        PRE: Assumes all levels with same symmetry point to same RelativisticConfigList.
     */
    inline void Print(const Symmetry& sym, double min_percentage, double max_energy) const
    {   Print(sym, min_percentage, true, max_energy);
    }

    bool Read(const std::string& filename, const std::string& angular_directory);
    void Write(const std::string& filename) const;

protected:
    /** Hidden version with all possible options for printing. All other print functions call this one. */
    void Print(const Symmetry& sym, double min_percentage, bool use_max_energy, double max_energy) const;
};

typedef std::shared_ptr<LevelMap> pLevelMap;
typedef std::shared_ptr<const LevelMap> pLevelMapConst;

#endif
