#ifndef ANGULAR_DATA_H
#define ANGULAR_DATA_H

#include "HartreeFock/OrbitalInfo.h"
#include "Projection.h"
#include <list>
#include <boost/shared_ptr.hpp>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <boost/iterator/iterator_facade.hpp>

class RelativisticConfiguration;

class AngularData
{
    /** Store projections and configuration state functions (CSF) corresponding to
        the angular part of a RelativisticConfiguration.
        Access is only provided by const iterators.
    */
public:
    AngularData(): CSFs(nullptr) {}
    AngularData(const RelativisticConfiguration& config, int two_m); //!< Generate projections but not CSFs.
    AngularData(const RelativisticConfiguration& config, int two_m, int two_j); //!< Generate projections and CSFs.
    ~AngularData();

    typedef std::vector< std::pair<int, int> > ConfigKeyType;   // pair(kappa, pqn) for all orbitals in RelativisticConfiguration
    typedef std::list< std::vector<int> >::const_iterator const_projection_iterator;
    typedef const double* const_CSF_iterator;

    /** Bidirectional list iterator over list of "projections", i.e. list of std::vector<int> */
    const_projection_iterator projection_begin() const { return projections.begin(); }
    const_projection_iterator projection_end() const { return projections.end(); }
    unsigned int projection_size() const { return projections.size(); }

    /** Random access iterator over CSFs corresponding to projection i.
        PRE: 0 <= i <= projection_size()
     */
    const_CSF_iterator CSF_begin(int i) const { return CSFs + i * num_CSFs; }
    const_CSF_iterator CSF_end(int i) const { return CSFs + (i+1) * num_CSFs; }

    int GetTwoM() const { return two_m; }
    int GetTwoJ() const { return two_j; }
    unsigned int NumCSFs() const { return num_CSFs; }

    /** Generate CSFs by diagonalising projections over J^2. */
    int GenerateCSFs(const RelativisticConfiguration& config, int two_j);
    int GenerateCSFs(const ConfigKeyType& key, int two_j);

protected:
    int GenerateProjections(const RelativisticConfiguration& config, int two_m);
    double GetJSquared(const Projection& first, const Projection& second) const;

    static bool ProjectionCompare(const std::vector<int>& first, const std::vector<int>& second);

    /** List of "projections": in this context, vectors of two_Ms. */
    std::list< std::vector<int> > projections;
    int two_m;

    /** CSF coefficients for a given J. Usually one requires all coefficients for a given projection,
        so the projection comes first: CSFs[proj * N + csf] where N is NumCSFs().
     */
    double* CSFs;
    unsigned int num_CSFs;
    int two_j;
};

typedef boost::shared_ptr<AngularData> pAngularData;
typedef boost::shared_ptr<const AngularData> pAngularDataConst;

class AngularDataLibrary
{
    /** Collection of AngularData elements, indexed by RelativisticConfiguration.
        The collection is stored on disk in
            AMBiT/AngularData/particle_number.two_m.two_j.angular
     */
public:
    AngularDataLibrary(int particle_number, int two_m, int two_j);
    ~AngularDataLibrary() {}

    /** Retrieve or create an AngularData object for the given configuration, two_m, and two_j.
        If there are no projections possible with our two_m, it returns a null pointer.
        PRE: config is sorted correctly.
     */
    pAngularData operator[](const RelativisticConfiguration& config);

    void GenerateCSFs();

    void Read();
    void Write() const;

protected:
    typedef AngularData::ConfigKeyType KeyType;

    int two_m, two_j;
    KeyType GenerateKey(const RelativisticConfiguration& config) const;

    std::string filename;
    std::unordered_map<KeyType, pAngularData, boost::hash<KeyType> > library;
};

typedef boost::shared_ptr<AngularDataLibrary> pAngularDataLibrary;

#endif
